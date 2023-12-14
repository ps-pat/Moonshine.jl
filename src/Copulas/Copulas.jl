import StatsAPI: loglikelihood, fit!

using Combinatorics: combinations

import Distributions

using Random: AbstractRNG, GLOBAL_RNG

using JuMP

using NLopt: NLopt

using FunctionWrappers: FunctionWrapper

###################
# Packaged Alphas #
###################

include("Alphas.jl")

export AbstractΦCopula
"""
    AbstractΦCopula{P, A, N}

Abstract type for copulas joining two phenotypes.

# Type Parameters

  - `P<:Distribution`: type of the phenotypes
  - `A<:AbstractAlpha`: type of the α
"""
abstract type AbstractΦCopula{P<:Distribution,A<:AbstractAlpha} end

#############
# Interface #
#############

"""
    pdf(copula, φ1[, φ2, d, αpars...])
    logpdf(copula, φ1[, φ2, d, αpars...])
    conditional_pdf(copula, φ1, φ2, d[, αpars...])
    logconditional_pdf(copula, φ1, φ2, d[, αpars...])

(Logarithm of) the (conditional) probability density function of the copula
given a genetic distance and a set of parameters for the associated α. If
`αpars...` is not provided, use the values stored in the α.

The one phenotype version of `pdf` and `logpdf` is the marginal (log)density
of that penotype. As such, these methods do not need nor accept `d` or `αpars`
arguments.

# Implementation

Only one of `logpdf` or `pdf` and one of `logconditional_pdf` or
`conditional_pdf` need a custom implementation. Two methods need implementation
for `pdf`:

  - `pdf(copula, φ)`
  - `pdf(copula, φ1, φ2, d, αpars)`
"""
function pdf_marginal end,
function logpdf_marginal end,
function pdf_joint end,
function logpdf_joint end,
function pdf_conditional end,
function logpdf_conditional end

for pdftype ∈ ["marginal", "joint", "conditional"]
    fun_str = "pdf_" * pdftype
    fun = Symbol(fun_str)
    logfun = Symbol("log" * fun_str)

    @eval begin
        export $fun, $logfun

        $fun(copula::AbstractΦCopula) = exp ∘ $logfun(copula)

        $logfun(copula::AbstractΦCopula) = log ∘ $fun(copula)
    end
end

## Marginal pdfs

function pdf_marginal(copula::AbstractΦCopula{<:Bernoulli})
    p = succprob(marginal(copula))

    φ -> φ ? p : 1 - p
end

export alpha
"""
    alpha(copula)

α function associated with the copula.

# Implementation

If the α is stored as a field named `α`, no custom implementation is needed.
"""
alpha(copula::AbstractΦCopula) = getfield(copula, :α)

export marginal
"""
    marginal(copula)

Marginal distribution associated with the copula.

# Implementation

If the distribution is stored as a field named `marginal_distribution`, no
custom impementation is needed.
"""
marginal(copula::AbstractΦCopula) = getfield(copula, :marginal_distribution)

#############
# Inference #
#############

export loglikelihood
"""
    loglikelihood(copula, Φ, genealogy)
    loglikelihood([rng, ]copula, Φ, H, G; idx = 1, n = 1000, genpars...)

Log-likelihood of a copula.

Both forms return functions that take one positional argument by parameter in
the order given by `(parameters ∘ alpha)(copula)`.

# Arguments

  - `rng::AbstractRNG`: random number generator
  - `copula::AbstractΦCopula`: copula for the phenotypes
  - `Φ`: iterable of phenotypes
  - `H`: iterable of genetic sequences
  - `G`: type of genealogy
  - `idx`: index of the marker for which to build genealogies
  - `n`: number of genealogy to sample
  - `genpars...`: genetic parameters passed as keyword arguments to the
    constructor of `G`
"""
function loglikelihood end

function loglikelihood(copula::AbstractΦCopula, Φ, genealogy::AbstractGenealogy)
    function (pars...)
        sum(combinations(leaves(genealogy), 2), init = zero(Float64)) do (i, j)
            logpdf(copula, Φ[i], Φ[j], distance(genealogy, i, j), pars...)
        end
    end
end

function compute_distances(tree::Tree, Φ::AbstractVector{Bool})
    n = nleaves(tree)

    ## 0-1 matrix of ancestors for leaves. Entry i, j is 1 iff leaf j
    ## has vertex n + i in its ancestry. Since every leaf has vertex
    ## 2n - 1 (the root) as an ancestor, corresponding row is not
    ## included in the matrix.
    ancestry = falses(n - 2, n)

    @inbounds for ivertex ∈ Iterators.drop(reverse(ivertices(tree)), 1)
        descendant_leaves = filter(<=(nleaves(tree)), descendants(tree, ivertex))
        ancestry[ivertex - n, descendant_leaves] .= true
    end

    ## mults stores the multiplicity of a distance for a given phenotypes pair.
    mults = zeros(Int, nivertices(tree), 3)
    dists = 2 .* latitudes(tree)

    @inbounds for (i, j) ∈ combinations(leaves(tree), 2)
        line = something(findfirst(ancestry[:,i] .& ancestry[:,j]), n - 1)
        col = Φ[i] + Φ[j] + 1
        mults[line, col] += 1
    end

    mults, dists
end

function loglikelihood(rng::AbstractRNG, copula::AbstractΦCopula{<:Bernoulli},
                       Φ, H, G; idx = 1, n = 1000, genpars...)
    npars = (length ∘ parameters)(alpha(copula))
    fs = Vector{FunctionWrapper{Real, NTuple{npars, Real}}}(undef, n)

    logpdf = logpdf_joint(copula)

    for k ∈ eachindex(fs)
        genealogy = G(H; genpars...)
        build!(rng, genealogy, idx)

        mults, dists = compute_distances(genealogy, Φ)

        fs[k] = function (pars::T...) where T<:Real
            sum((enumerate ∘ eachcol)(mults), init = zero(promote_type(T, Float64))) do (k, mults_col)
                ## Phenotypes are Gray-encoded in k - 1.
                k -= 1
                φ1, φ2 = .!iszero.((k ⊻ (k >> 1)) .& (1, 2))

                sum(zip(mults_col, dists), init = zero(promote_type(T, Float64))) do (mult, dist)
                    iszero(mult) && return zero(dist)

                    mult * logpdf(φ1, φ2, dist, pars...)
                end
            end
        end
    end

    function (pars::T...) where T<:Real
        sum(f -> f(pars...)::promote_type(T, Float64), fs,
            init = zero(promote_type(T, Float64)))
    end
end

function loglikelihood(copula::AbstractΦCopula, Φ, H, G; kwargs...)
    loglikelihood(GLOBAL_RNG, copula, Φ, H, G; kwargs...)
end

export fit!
"""
    fit!(copula, Φ, H, G; global_attrs, local_attrs, genpars...)

Fit a copula. The domain of the parameters is assumed to be unbounded unless
stated otherwise by `bounds(alpha(colupa))`

# Implementation

The default implementation fits parameters by maximizing the loglikelihood of
the copula. The optimization is done in two passes:

 1. global optimizations using ESCH (an evolutionary algorithm);
 2. local optimization using LBFGS.
    Both implementation are from the [NLopt](https://nlopt.readthedocs.io/en/latest/)
    library. Since only finite domain is supported by ESCH, a call to
    `bound(::T)` must return a bounded domain for each parameter of the likelihood
    in order to use the default `fit!(::AbstractΦCopula, ...)` method.

## Arguments

The default implementation accepts two keyword arguments `global_attrs` and
`local_attrs` which must be iterable of attribute/value pairs. They are used to
set the attributes of the global/local optimizer respectively. See
[NLopt.jl documentation](https://github.com/JuliaOpt/NLopt.jl). Their
default values are:

  - `global_attrs = ("algorithm" => :GN_ESCH, "maxtime" => 5 * (length ∘ parameters)(alpha(copula)), "maxeval" => 2000)`
  - `local_attrs = ("algorithm" => :LN_NELDERMEAD, "maxtime" => 5 * (length ∘ parameters)(alpha(copula)))`
"""
function fit!(rng, copula::AbstractΦCopula, Φ, H, G;
              global_attrs = ("algorithm" => :GN_ESCH,
                              "maxeval" => 1000),
              local_attrs = ("algorithm" => :LN_SBPLX,
                             "maxtime" => ceil(Int, length(H) / 10) * (length ∘ parameters)(alpha(copula))),
              genpars...)
    α = alpha(copula)

    ## Models
    global_model, local_model = Model(NLopt.Optimizer), Model(NLopt.Optimizer)

    ## Add variables to model.
    _bounds = bounds(α)

    ## It is mandatory for variables to be added to the model in the order
    ## they appear in the signature of the objective function.
    for arg ∈ parameters(α)
        if arg ∈ keys(_bounds)
            @variable(global_model, base_name = string(arg),
                      lower_bound = first(_bounds[arg]),
                      upper_bound = last(_bounds[arg]))
            @variable(local_model, base_name = string(arg),
                      lower_bound = first(_bounds[arg]),
                      upper_bound = last(_bounds[arg]))
        else
            @variable(global_model, base_name = string(arg))
            @variable(local_model, base_name = string(arg))
        end

        set_start_value(variable_by_name(global_model, string(arg)),
                        getparameter(α, arg))
    end
    global_vars, local_vars = all_variables(global_model), all_variables(local_model)

    ## Objective Function
    objective = loglikelihood(rng, copula, Φ, H, G; genpars...)
    nparameters = length(parameters(α))
    register(global_model, :objective,
             nparameters, objective,
             autodiff = true)
    register(local_model, :objective,
             nparameters, objective,
             autodiff = true)
    @NLobjective(global_model, Max, objective(global_vars...))
    @NLobjective(local_model, Max, objective(local_vars...))

    ## Global Optimization
    set_attributes(global_model, global_attrs...)
    optimize!(global_model)

    termination_status(global_model) ∈
    (LOCALLY_SOLVED, OPTIMAL, TIME_LIMIT, ITERATION_LIMIT) ||
        @warn "Global optimization did not converge"

    ## Local Optimization
    for (global_var, local_var) ∈ zip(global_vars, local_vars)
        set_start_value(local_var, value(global_var))
    end

    set_attributes(local_model, local_attrs...)
    optimize!(local_model)

    termination_status(local_model) ∈ (LOCALLY_SOLVED, OPTIMAL) ||
        @warn "Local optimization did not converge"

    ## Store estimated values.
    for (var, val) ∈ zip(Symbol.(local_vars), value.(local_vars))
        setparameter!(α, var, val)
    end

    (global_model = global_model, local_model = local_model)
end

function fit!(copula::AbstractΦCopula, Φ, H, G; kwargs...)
    fit!(GLOBAL_RNG, copula, Φ, H, G; kwargs...)
end

####################
# Packaged Copulas #
####################

macro copula_struct(copula)
    ## Generate name
    copula_string = string(copula)
    if length(copula_string) < 6 || copula_string[1:6] != "Copula"
        copula = Symbol("Copula" * copula_string)
    end

    quote
        export $copula
        struct $copula{P,A} <: AbstractΦCopula{P,A}
            marginal_distribution::P
            α::A

            #$copula(P, α) = new{typeof(P),typeof(α)}(P, α)
        end
    end
end

include("Frechet.jl")
include("CuadrasAuge.jl")
