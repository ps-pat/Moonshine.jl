using LinearAlgebra: Symmetric

using StaticArrays

using NLPModels

export AbstractAlpha
"""
    AbstractAlpha

Abstract type for αs.
"""
abstract type AbstractAlpha end

#############
# Interface #
#############

export bounds
"""
    bounds(α)

Bounds for the parameters of an α.

# Implementation

Must return a named tuple. Keys must be symbols representing parameters. The
corresponding value must be a container of length 2 where the first and last
entries are lower and upper bounds for the parameter respectively.
"""
function bounds end

export parameters
"""
    parameters(α)

Parameters to be optimized as an iterable of symbols.
"""
function parameters end

export getparameter
"""
    getparameter(α, parameter)

Get the value of a parameter.
"""
getparameter(α::AbstractAlpha, parameter) = getfield(α, parameter)::Float64

export setparameter!
"""
    setparameter!(α, parameter, value)

Set the value of a parameter.
"""
setparameter!(α::AbstractAlpha, parameter, value) =
    setfield!(α, parameter, value)

export nlpmeta
"""
    nlpmeta(α)

NLPMeta object.
[Parameters here](https://github.com/JuliaSmoothOptimizers/NLPModels.jl/blob/main/src/nlp_types.jl#L28)
"""
function nlpmeta end

export nlppars
"""
    nlppars(α)

Parameters passed to ipopt.
"""
nlppars(::AbstractAlpha) = (; mu_strategy = "adaptive",
                            adaptive_mu_globalization = "kkt-error")

export ∇, ∇!, ∇², ∇²!

###############
# Exponential #
###############

export AlphaExponential

mutable struct AlphaExponential <: AbstractAlpha
    λ::Float64
    λbounds::NTuple{2, Float64}
end

## Constructors.
AlphaExponential(λ) = AlphaExponential(λ, (1e-6, 10))

AlphaExponential() = AlphaExponential(1)

## Evaluation
function(::AlphaExponential)(t::T, λvec) where T
    λ = first(λvec)
    λ ≤ 0 && return zero(T)

    -expm1(-first(λvec) * t)
end

(α::AlphaExponential)(t) = α(t, [α.λ])

## Gradient
function ∇!(::AlphaExponential, g, t, λvec, σ = 1.)
    g[1] += t * exp(-first(λvec) * t) * σ
    g
end

function ∇(α::AlphaExponential, t, λvec, σ = 1.)
    g = zeros(MVector{1})
    ∇!(α, g, t, λvec, σ)
end

## Hessian
function ∇²!(α::AlphaExponential, vals, t, λvec, mult = 1., σ = 1., σ2 = 1.)
    λ = first(λvec)

    vals[1] += (-t^2 * exp(-λ * t) * σ + (σ2 - σ^2) * t^2 * exp(-2λ * t)) *
        mult

    vals
end

function ∇²(α::AlphaExponential, t, αpars, mult = 1., σ = 1., σ2 = 1.)
    vals = zeros(MVector{1})
    ∇²!(α, t, αpars, mult, σ, σ2)
end

@generated parameters(::AlphaExponential) = (:λ,)

bounds(α::AlphaExponential) = (; λ = α.λbounds)

nlpmeta(α::AlphaExponential) = NLPModelMeta(
    1;
    x0 = [α.λ],
    lvar = [first(α.λbounds)], uvar = [last(α.λbounds)],
    minimize = false
)

const AOExponential = AlphaOptimization{<:AbstractΦCopula{<:Any, <:AlphaExponential}}

function NLPModels.hess_structure!(::AOExponential, rows, cols)
    rows[1] = 1
    cols[1] = 1

    rows, cols
end

############
# Gompertz #
############

export AlphaGompertz

mutable struct AlphaGompertz <: AbstractAlpha
    η::Float64
    b::Float64
    ηbounds::NTuple{2, Float64}
    bbounds::NTuple{2, Float64}
end

## Constructors.
AlphaGompertz(η, b) = AlphaGompertz(η, b, (1e-6, 10), (1e-6, 10))

AlphaGompertz() = AlphaGompertz(1., 1.)

## Evaluation
function (::AlphaGompertz)(t::T, αpars) where T
    η, b = αpars
    η ≤ 0 || b ≤ 0 && return zero(T)

    -expm1(-η * expm1(b * t))
end

(α::AlphaGompertz)(t) = α(t, (α.η, α.b))

## Gradient
function ∇!(::AlphaGompertz, g, t, αpars, σ = 1.)
    η, b = αpars

    p = b * t
    q = expm1(p)
    s = exp(-η * q)

    g[1] += q * s * σ
    g[2] += η * t * exp(p) * s * σ
    g
end

function ∇(α::AlphaGompertz, t, αpars, σ = 1.)
    g = zeros(MVector{2})
    ∇!(α, g, t, αpars, σ)
end

## Hessian
function ∇²!(α::AlphaGompertz, vals, t, αpars, mult = 1., σ = 1., σ2 = 1.)
    η, b = αpars

    ∇α = ∇(α, t, αpars)
    p = b * t
    q = expm1(p)
    s = exp(-η * q)

    vals[1] += (-q^2 * s *
        σ + ∇α[1]^2 * (σ2 - σ^2)) * mult
    vals[2] += (t * exp(p) * s * (1 - η * q) *
        σ + ∇α[1] * ∇α[2] * (σ2 - σ^2)) * mult
    vals[3] += (η * t^2 * exp(p) * s * (1 - η * exp(b * t)) *
        σ + ∇α[2]^2 * (σ2 - σ^2)) * mult

    vals
end

function ∇²(α::AlphaGompertz, t, αpars, mult = 1., σ = 1., σ2 = 1.)
    vals = zeros(MVector{3})
    ∇²!(α, t, αpars, mult, σ, σ2)
end

@generated parameters(::AlphaGompertz) = (:η, :b)

bounds(α::AlphaGompertz) = (; η = α.ηbounds, b = α.bbounds)

nlpmeta(α::AlphaGompertz) = NLPModelMeta(
    2;
    x0 = [α.η, α.b],
    lvar = [first(α.ηbounds), first(α.bbounds)],
    uvar = [last(α.ηbounds), last(α.bbounds)],
    minimize = false
)

const AOGompertz = AlphaOptimization{<:AbstractΦCopula{<:Any, <:AlphaGompertz}}

function NLPModels.hess_structure!(::AOGompertz, rows, cols)
    rows .= (1, 2, 2)
    cols .= (1, 1, 2)

    rows, cols
end

################
# Gumponential #
################

export AlphaGE

mutable struct AlphaGE <: AbstractAlpha
    p::Float64
    const αE::AlphaExponential
    const αG::AlphaGompertz
end

## Constructors.
AlphaGE(p, λ, η, b) = AlphaGE(p, AlphaExponential(λ), AlphaGompertz(η, b))

AlphaGE() = AlphaGE(0.5, AlphaExponential(), AlphaGompertz())

## Evaluation
function (α::AlphaGE)(t, αpars)
    p, λ, η, b = αpars
    p * α.αE(t, (λ,)) + (1 - p) * α.αG(t, (η, b))
end

function (α::AlphaGE)(t)
    p = α.p
    p * α.αE(t) + (1 - p) * α.αG(t)
end

## Gradient
function ∇!(α::AlphaGE, g, t, αpars, σ = 1.)
    αE = α.αE
    αG = α.αG
    p = α.p

    g[1] = (αE(t, αpars) - αG(t, αpars)) * σ
    ∇!(αE, view(g, 2), t, view(αpars, 2), σ * p)
    ∇!(αG, view(g, 3:4), t, view(αpars, 3:4), σ * (1 - p))

    g
end

function ∇(α::AlphaGE, t, αpars, σ = 1.)
    g = zeros(MVector{4})
    ∇!(α, g, t, αpars, σ)
end

## Hessian
function ∇²!(α::AlphaGE, vals, t, αpars, mult = 1., σ = 1., σ2 = 1.)
    αE = α.αE
    αG = α.αG
    p = α.p

    ∇α = ∇(α, t, αpars)
    ∇αE = ∇(αE, t, view(αpars, 2))
    ∇αG = ∇(αG, t, view(αpars, 3:4))

    vals[1] = (∇αE[1] * σ + ∇α[1] * ∇α[2] * (σ2 - σ^2)) * mult
    vals[2] = (-∇αG[1] * σ + ∇α[1] * ∇α[3] * (σ2 - σ^2)) * mult
    vals[3] = (-∇αG[2] * σ + ∇α[1] * ∇α[4] * (σ2 - σ^2)) * mult
    ∇²!(αE, view(vals, 4), t, view(αpars, 2), mult, σ * p, σ2 * p^2)
    ∇²!(αG, view(vals, 5:7), t, view(αpars, 3:4), mult, σ * p, σ2 * p^2)

    vals
end

function ∇²(α::AlphaGE, t, αpars, mult = 1., σ = 1., σ2 = 1.)
    vals = zeros(MVector{7})
    ∇²!(α, vals, t, αpars, mult, σ, σ2)
end

@generated parameters(::AlphaGE) = (:p, :λ, :η, :b)

bounds(α::AlphaGE) = (p = (1e-3, 1 - 1e-3), bounds(α.αE)..., bounds(α.αG)...)

function nlpmeta(α::AlphaGE)
    nlpE = nlpmeta(α.αE)
    nlpG = nlpmeta(α.αG)

    NLPModelMeta(4;
                 x0 = [α.p; nlpE.x0; nlpG.x0],
                 lvar = [first(bs) for bs ∈ bounds(α)],
                 uvar = [last(bs) for bs ∈ bounds(α)],
                 minimize = false,
                 nnzh = 7)
end

nlppars(α::AlphaGE) = (
    max_iter = 100,
    mu_strategy = "adaptive",
    #adaptive_mu_globalization = "kkt-error",
    nlp_scaling_max_gradient = 1e-5,
    nlp_scaling_min_value = 1e-6,
    tol = 1e-5,
    dual_inf_tol = 1e4
    #compl_inf_tol = 1e1,
)

const AOGE = AlphaOptimization{<:AbstractΦCopula{<:Any, <:AlphaGE}}

function NLPModels.hess_structure!(::AOGE, rows, cols)
    rows .= (2, 3, 4, 2, 3, 4, 4)
    cols .= (1, 1, 1, 2, 3, 3, 4)

    rows, cols
end
