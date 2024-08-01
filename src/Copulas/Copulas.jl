import StatsAPI: fit!

using Combinatorics: combinations

import Distributions

using Random: AbstractRNG, GLOBAL_RNG

using Metaheuristics: ECA, Options, optimize, minimizer

using NLPModels, NLPModelsIpopt

using StrideArrays

abstract type AbstractAlpha end

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

export pdf_marginal, pdf_joint, pdf_conditional
"""
    pdf_marginal(copula, φ; logpdf = false)
    pdf_joint(copula, φ, ψ, t, αpars; logscale = false)
    pdf_conditional(copula, φ, ψ, t, logscale = false)

Marginal/joint/conditional probability density
function of the copula.
"""
function pdf_marginal end,
function pdf_joint end,
function pdf_conditional end

## Marginal pdfs

function pdf_marginal(copula::AbstractΦCopula{<:Bernoulli}, φ; logscale = false)
    p = succprob(marginal(copula))
    ret = φ ? p : 1 - p

    logscale ? log(ret) : ret
end

## Gradient & Hessian

export ∇scale, ∇²scale
"""
    ∇scale
    ∇²scale

Functions evaluating to the scaling factors used for the application
of the chain rule. Specifically:

  - `∇scale`: f'(α) / f(α)
  - `∇scale`: f''(α) / f(α)
"""
function ∇scale end,
function ∇²scale end

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

function compute_mults(tree::Tree, Φ::AbstractVector{Bool})
    n = nleaves(tree)

    ## 0-1 matrix of ancestors for leaves. Entry i, j is 1 iff leaf j
    ## has vertex n + i in its ancestry. Since every leaf has vertex
    ## 2n - 1 (the root) as an ancestor, corresponding row is not
    ## included in the matrix. First dimension is padded so that its
    ## length is a multiple of 64.
    nrow = n - 2 + 64 - (n - 2) % 64
    nchunks_bycol = nrow ÷ 64
    ancestry = falses(nrow, n)

    @inbounds for ivertex ∈ Iterators.drop(reverse(ivertices(tree)), 1)
        for descendant ∈ descendants(tree, ivertex)
            isleaf(tree, descendant) || continue
            ancestry[ivertex - n, descendant] = true
        end
    end

    ## mults stores the multiplicity of a distance for a given phenotypes pair.
    mults = StrideArray{Int}(undef, nivertices(tree), StaticInt(3))
    fill!(mults, 0)

    @inbounds for i ∈ range(0, length = n - 1)
        firsti = div(i * nrow + 1, 64, RoundUp)
        for j ∈ range(i + 1, n - 1)
            firstj = div(j * nrow + 1, 64, RoundUp)
            line = n - 1
            col = Φ[i + 1] + Φ[j + 1] + 1

            for k ∈ range(0, length = nchunks_bycol)
                p = trailing_zeros(
                    ancestry.chunks[firsti + k] &
                        ancestry.chunks[firstj + k]) + 1
                if p < 64
                    line = 64k + p
                    break
                end
            end
            mults[line, col] += 1
        end
    end

    mults
end

export AlphaOptimization
mutable struct AlphaOptimization{C<:AbstractΦCopula{<:Bernoulli}} <: AbstractNLPModel{Float64, Vector{Float64}}
    meta::NLPModelMeta{Float64, Vector{Float64}}
    counters::Counters

    const multiplicities::Vector{Matrix{Int}}
    const distances::Vector{Vector{Float64}}
    const copula::C
end

## Packaged Alphas
include("Alphas.jl")

function AlphaOptimization(rng, copula::AbstractΦCopula{<:Bernoulli},
                           Φ, H, ::Type{Tree};
                           n = 10,
                           treepars = (;), genpars...)
    nrow = length(H) - 1
    multiplicities = Vector{Matrix{Int}}(undef, n)
    distances = Vector{Vector{Float64}}(undef, n)

    for k ∈ 1:n
        genealogy = Tree(H; genpars...)
        build!(rng, genealogy, treepars...)

        distances[k] = 2 * latitudes(genealogy)
        multiplicities[k] = compute_mults(genealogy, Φ)
    end

    AlphaOptimization(nlpmeta(alpha(copula)), Counters(),
                      multiplicities, distances, copula)
end

function decode_phenos(k)
    k -= 1
    .!iszero.((k ⊻ (k >> 1)) .& (1, 2))
end

function fun_bernoulli!(ret, multsvec, distsvec, f)
    @inbounds for (mults, dists) ∈ zip(multsvec, distsvec)
        for j ∈ axes(mults, 2)
            φ, ψ = decode_phenos(j)

            for (i, dist) ∈ enumerate(dists)
                mult = mults[i, j]

                iszero(mult) && continue
                ret += mult * f(φ, ψ, dist)
            end
        end
    end

    ret
end

function NLPModels.obj(nlp::AlphaOptimization{<:AbstractΦCopula{<:Bernoulli}},
                       αpars::AbstractVector{T}) where T
    ret = zero(T)
    f = (φ, ψ, t) -> pdf_joint(nlp.copula, φ, ψ, t, αpars, logscale = true)
    fun_bernoulli!(ret, nlp.multiplicities, nlp.distances, f)
end

function fun_bernoulli!!(multsvec, distsvec, f!)
    @inbounds for (mults, dists) ∈ zip(multsvec, distsvec)
        for j ∈ axes(mults, 2)
            φ, ψ = decode_phenos(j)

            for (i, dist) ∈ enumerate(dists)
                mult = mults[i, j]

                iszero(mult) && continue
                f!(φ, ψ, dist, mult)
            end
        end
    end
end

function NLPModels.grad!(nlp::AlphaOptimization{<:AbstractΦCopula{<:Bernoulli}},
                         αpars, g)
    fill!(g, 0)
    α = alpha(nlp.copula)

    f! = function (φ, ψ, t, mult)
        σ = ∇scale(nlp.copula)(φ, ψ, t, αpars)
        ∇!(α, g, t, αpars, mult * σ)
    end

    fun_bernoulli!!(nlp.multiplicities, nlp.distances, f!)
    g
end

function NLPModels.hess_coord!(nlp::AlphaOptimization{<:AbstractΦCopula{<:Bernoulli}},
                               αpars::AbstractVector, vals::AbstractVector;
                               obj_weight = 1)
    fill!(vals, 0)
    α = alpha(nlp.copula)

    f! = function(φ, ψ, t, mult)
        σ = ∇scale(nlp.copula)(φ, ψ, t, αpars)
        σ2 = ∇²scale(nlp.copula)(φ, ψ, t, αpars)
        ∇²!(α, vals, t, αpars, mult, σ, σ2)
    end

    fun_bernoulli!!(nlp.multiplicities, nlp.distances, f!)
    vals .*= obj_weight

    vals
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

function fit!(rng, copula::AbstractΦCopula{<:Bernoulli},
              Φ, H, ::Type{Tree};
              n = 10,
              linsolver = "mumps",
              libpardisopath = "",
              global_attrs = (;), local_attrs = (;),
              genpars...)
    α = alpha(copula)
    npars = (length ∘ parameters)(α)

    model = AlphaOptimization(rng, copula, Φ, H, Tree; n = n, genpars...)

    #######################
    # Global Optimization #
    #######################

    global_algorithm = ECA(; global_attrs...,
                           options = Options(rng = rng))

    obj_global = function (αpars)
        f = -obj(model, αpars)
        cs = cons(model, αpars)
        g = [cs - model.meta.ucon; model.meta.lcon - cs]

        f, g, Float64[]
    end

    global_bounds = hcat(model.meta.lvar, model.meta.uvar)'
    res_global = optimize(obj_global, global_bounds, global_algorithm)

    println(res_global)

    ######################
    # Local Optimization #
    ######################

    res_local = ipopt(model;
                      x0 = minimizer(res_global),
                      check_derivatives_for_naninf = "yes",
                      honor_original_bounds = "yes",
                      max_iter = 1000,
                      nlppars(α)...,
                      local_attrs...)

    for (k, parameter) ∈ enumerate(parameters(α))
        setparameter!(α, parameter, res_local.solution[k])
    end

    (global_model = res_global, local_model = res_local)
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

## Packaged Copulas
include("Frechet.jl")
include("CuadrasAuge.jl")
