import StatsAPI: fit!

using Combinatorics: combinations

import Distributions

using Random: AbstractRNG, GLOBAL_RNG

using JuMP, Ipopt

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
    pdf_marginal(copula,)
    logpdf_marginal(copula)
    pdf_joint(copula)
    logpdf_joint(copula)
    pdf_conditional(copula)
    logpdf_conditional(copula)

(Logarithm of) the marginal/joint/conditional probability density
function of the copula.

# Implementation

- Only one of `logpdf_joint` or `pdf_joint` need to be implemented.
- Only one of `logpdf_conditional` or `pdf_conditional` need to be implemented.
- Implemented methods should return a function taking 3 positional arguments and
  a vararg argument. The first 2 arguments are the phenotypes. The third
  argument is the genetic distance between the individuals to which the
  phenotypes are associated. The remaining arguments are to be passed to the
  α function associated with the copula.
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

## Gradient

export ∇logpdf_joint
export ∇²logpdf_joint
"""
    ∇logpdf_joint(copula)

Logarithmic derivatives of the joint pdf.

# Implementation

If `∇logpdf_joint` or `∇²logpdf_joint` are not implemented, automatic
differentiation will be used by [`fit!`](@ref) if needed.
"""
function ∇logpdf_joint end,
function ∇²logpdf_joint end

∇logpdf_joint(copula::AbstractΦCopula) = nothing
∇²logpdf_joint(copula::AbstractΦCopula) = nothing

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
        @views line = something(findfirst(ancestry[:,i] .& ancestry[:,j]), n - 1)
        col = Φ[i] + Φ[j] + 1
        mults[line, col] += 1
    end

    mults, dists
end

function decode_phenos(k)
    k -= 1
    .!iszero.((k ⊻ (k >> 1)) .& (1, 2))
end

fun_bernoulli(mults, dists, f, ::Val{N}) where N =
    isone(N) ? fun_bernoulli(mults, dists, f) :
    fun_bernoulli_vector(mults, dists, f, Val(N))

fun_bernoulli(mults, dists, f) = function (pars::T...) where T<:Real
    z0 = zero(promote_type(T, Float64))

    sum((enumerate ∘ eachcol)(mults), init = z0) do (k, mults_col)
        ## Phenotypes are Gray-encoded in k - 1.
        φ1, φ2 = decode_phenos(k)

        sum(zip(mults_col, dists), init = z0) do (mult, dist)
            iszero(mult) && return z0

            mult .* f(φ1, φ2, dist, pars...)
        end
    end
end

fun_bernoulli_vector(mults, dists, f, npars) = function (pars::T...) where T<:Real
    z = fill(z0, npars)

    sum((enumerate ∘ eachcol)(mults), init = z) do (k, mults_col)
        ## Phenotypes are Gray-encoded in k - 1.
        φ1, φ2 = decode_phenos(k)

        sum(zip(mults_col, dists), init = z) do (mult, dist)
            iszero(mult) && return z

            mult .* f(φ1, φ2, dist, pars...)
        end
    end
end

export fit!
"""
    fit!([rng, ]copula, Φ, H, G; n = 10, linsolver = "mumps",
         global_attrs, local_attrs, genpars...)

Fit a copula. The domain of the parameters is assumed to be unbounded unless
stated otherwise by `bounds(alpha(colupa))`. Optimization is done in two passes:

 1. global optimizations using DIRECT;
 2. local optimization using Ipopt.

The algorithms are implemented by
[NLopt](https://nlopt.readthedocs.io/en/latest/) and
[Ipopt](https://github.com/coin-or/Ipopt) respectively.

# Arguments

 - `rng`: Random number generator
 - `copula`: Copula to fit
 - `Φ`: Phenotypes
 - `H`: Haplotypes
 - `G`: Type of genealogy
 - `n`: Number of genealogies to sample
 - `linsolver`: Linear solver used by Ipopt
 - `global_attrs`: attributes for the global optimizer
 - `local_attrs`: attributes for the local optimizer
 - `genpars`: parameters passed directly to the constructor of the genealogies.

# Linear Solver

Two linear solvers for Ipopt are available out of the box for Julia ≥
1.9: MUMPS and SPRAL. If you want to use the latter, both
OMP_CANCELLATION and OMP_PROC_BIND environment have to be set to TRUE.

# Default attributes

Overriden by the `global_attrs` and `local_attrs` arguments.

## Global solver:

 - `"maxtime" => 1`

A complete list is available on
[NLopt.jl](https://github.com/JuliaOpt/NLopt.jl).

## Local solver:
 - `"check_derivatives_for_naninf" => "yes"`

A complete list is available on
[Ipopt.jl](https://github.com/jump-dev/Ipopt.jl).
"""
function fit! end

function fit!(rng, copula::AbstractΦCopula, Φ, H, G;
              n = 100,
              linsolver = "mumps",
              global_attrs = (), local_attrs = (),
              genpars...)
    α = alpha(copula)

    global_model, local_model = Model(NLopt.Optimizer), Model(Ipopt.Optimizer)

    ## Add variables to the model.
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

        set_start_value(something(variable_by_name(global_model, string(arg))),
                        getparameter(α, arg))
    end

    ## Objective function.

    ## Create n JuMP operators and assign them to both models. Each
    ## operator is the log-likelihood of a randomly sampled genealogy.
    let logpdf = logpdf_joint(copula),
        ∇logpdf = ∇logpdf_joint(copula),
        ∇²logpdf = ∇²logpdf_joint(copula),
        npars = (length ∘ parameters)(α)

        global_model[:logliks] = Vector{NonlinearOperator}(undef, n)
        local_model[:logliks] = Vector{NonlinearOperator}(undef, n)
        legacy_∇fs = Vector{Function}(undef, n)
        legacy_∇²fs = Vector{Function}(undef, n)

        for k ∈ 1:n
            genealogy = G(H; genpars...)
            build!(rng, genealogy)
            mults, dists = compute_distances(genealogy, Φ)

            opname = Symbol("l" * string(k))
            f = fun_bernoulli(mults, dists, logpdf)
            ∇f = fun_bernoulli(mults, dists, ∇logpdf, Val(npars))
            legacy_∇fs[k] = ∇f
            ∇²f = fun_bernoulli(mults, dists, ∇²logpdf, Val(npars))
            legacy_∇²fs[k] = ∇²f

            global_model[:logliks][k] =
                add_nonlinear_operator(local_model, npars, f, ∇f, ∇²f, name = opname)
            local_model[:logliks][k] =
                add_nonlinear_operator(local_model, npars, f, ∇f, ∇²f, name = opname)
        end

        ## Support for legacy interface :(
        legacy_f = x -> sum(f -> f(x), global_model[:logliks])
        legacy_∇f = x -> sum(f -> f(x), legacy_∇fs)
        legacy_∇²f = x -> sum(f -> f(x), legacy_∇²fs)

        register(global_model, :f, npars, legacy_f, legacy_∇f, legacy_∇²f)
    end

    global_vars, local_vars = all_variables(global_model), all_variables(local_model)
    @NLobjective(global_model, Max, f(global_vars...))
    @objective(local_model, Max, sum(f -> f(local_vars...), local_model[:logliks]))


    ## Global optimization.
    ## TODO: Switch to STOGO
    set_attributes(global_model,
                   "algorithm" => :GN_DIRECT,
                   "maxtime" => 1,
                   global_attrs...)
    optimize!(global_model)

    ## Local optimization. Warm start from the global optimization step.
    for (global_var, local_var) ∈ zip(global_vars, local_vars)
        set_start_value(local_var, value(global_var))
    end

    set_attributes(local_model,
                   "linear_solver" => linsolver,
                   "check_derivatives_for_naninf" => "yes",
                   local_attrs...)
    optimize!(local_model)

    ## Store estimations.
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
