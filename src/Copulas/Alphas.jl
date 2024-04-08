using Random: AbstractRNG, GLOBAL_RNG

using SpecialFunctions: erf, erfc

using Distributions: BetaPrime, Chi, Chisq, cdf

using StatsFuns

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

Must return a named tuple. The key of an entry must be the name of an
argument of `loglikelihood(::T, ...)`. The corresponding value must be a
container of length 2 where the first and last entries are the lower and
upper bounds respectively. The bounds
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
function getparameter end

export setparameter!
"""
    setparameter!(α, parameter, value)

Set the value of a parameter.
"""
function setparameter! end

#################
# Alpha Factory #
#################

let
    zeroinf = (1e-4, 1e2)
    zeroone = (1e-4, 1 - 1e-4)

    ## Maxwell-Boltzmann
    maxwellboltzmann_cdf = function (t, a)
        erf(t * invsqrt2 / a) - inv(sqrthalfπ) * (t / a) * exp(-t^2 * 0.5 / a^2)
    end

    ## Folded Normal
    foldednormal_cdf = function (t, μ, σ)
        c = invsqrt2 / σ

        0.5 * (erf((t + μ) * c) + erf((t - μ) * c))
    end

    ## Gompertz
    ∇gompertz = function (t, b, η)
        p = b * t
        q = expm1(p)
        s = exp(-η * q)

        [η * t * s * exp(p), q * s]
    end

    #! format: off
    As = (
        Exponential =
            ((t, λ) -> -expm1(-λ * t),
             (t, λ) -> t * exp(-λ * t),
             (t, λ) -> -t^2 * exp(-λ * t),
             (; λ = ("rate", 1, zeroinf)),
             "CDF of an exponential random variable."),)
        # Gompertz =
        #     ((t, b, η) -> -expm1(-η * expm1(b * t)), ∇gompertz,
        #      (b = ("scale", 1, zeroinf), η = ("shape", 1, zeroinf)),
        #      "CDF of a shifted Gompertz distributed random variable."))

#     As = (
#         MaxwellBoltzmann =
#             (maxwellboltzmann_cdf,
#              (; a = ("scale", 1, zeroinf)),
#              "CDF of a Maxwell-Boltzmann distributed random variable."),
#         BetaPrime =
#             ((t, a, b) -> cdf(BetaPrime(a, b), t),
#              (a = ("shape", 1, zeroinf), b = ("shape", 1, zeroinf)),
#              "CDF of a beta prime random variable."),
#         Dagum =
#             ((t, p, a, b) -> inv(1 + (b / t)^a)^p,
#              (p = ("shape", 1, zeroinf),
#               a = ("shape", 1, zeroinf),
#               b = ("scale", 1, zeroinf)),
#              """
# CDF of a Dagum distributed random variable.

# See also [`AlphaLogLogistic`](@ref) for the p = 1 constrained case.
#     """),
#         LogLogistic =
#             ((t, a, b) -> inv(1 + b * inv(t)^a),
#              (a = ("shape", 1, (1e-1, 1e1)),
#               b = ("scale", 1, (1e-1, 1e1))),
#              """
# CDF of a log-logistic distributed random variable.

# See also [`AlphaDagum`](@ref) for a 3 parameters parametrization.
#     """),
#         Pareto =
#             ((t, σ, ξ) -> 1 - (1 + ξ * t / σ)^(-inv(ξ)),
#              (σ = ("scale", 1, zeroinf), ξ = ("scale", 1, zeroinf)),
#              """
# CDF of a generalized Pareto distributed random variable with μ = 0.

# See also [`AlphaLomax`](@ref).
#     """),
#         FoldedNormal =
#             (foldednormal_cdf,
#              (μ = ("location", 0, (0, last(zeroinf))),
#               σ = ("scale", 1, zeroinf)),
#              "CDF of a folded-normal distributed random variable"),
#         Lomax =
#             ((t, a, λ) -> 1 - (1 + t / λ)^(-a),
#              (a = ("shape", 1, zeroinf), λ = ("scale", 1, zeroinf)),
#              """
# CDF of a Lomax distributed random variable. Special case of the generalized Pareto distribution.

# See also [`AlphaPareto`](@ref).
#     """),
#         Rayleigh =
#             ((t, σ) -> 1 - exp(-t^2 / (2 * σ^2)),
#              (;σ = ("scale", 1, zeroinf)),
#              "CDF of a Rayleigh distributed random variable."),
#         Chi =
#             ((t, k) -> cdf(Chi(k), t),
#              (;k = ("degrees of freedom", 1, zeroinf)),
#              "CDF of a chi distributed random variable."),
#         Chisq =
#             ((t, k) -> cdf(Chisq(k), t),
#              (;k = ("degrees of freedom", 1, zeroinf)),
#              "CDF of a chi-squared distributed random variable."),
#         ExpLog =
#             ((t, p, β) -> 1 - log(1 - (1 - p) * exp(-β * t)) / log(p),
#              (p = ("probability", 0.5, zeroone),
#               β = ("rate", 1, zeroinf)),
#              "CDF of an exponential-logarithm distributed random variable."),
#         LogCauchy =
#             ((t, μ, σ) -> inv(π) * atan((log(t) - μ) / σ) + 0.5,
#              (μ = ("location", 0, (-Inf, Inf)),
#               σ = ("scale", 1, zeroinf)),
#              "CDF of a log-Cauchy distributed random variable."),
#         Levy =
#             ((t, c) -> erfc(sqrt(c / (2 * t))),
#              (; c = ("scale", 1, zeroinf)),
#              "CDF of a Levy distributed random variable."),
#         Gompertz =
#             ((t, b, η) -> -expm1(-η * expm1(b * t)),
#              (t, b, η) -> b * η * exp(η + b * t - η * exp(b * t)),
#              (b = ("scale", 1, zeroinf), η = ("shape", 1, zeroinf)),
#              "CDF of a shifted Gompertz distributed random variable."),
#         GompertzShifted =
#             ((t, b, η) -> -expm1(-b * t) * exp(-η * exp(-b * t)),
#              (b = ("scale", 1, zeroinf), η = ("shape", 1, zeroinf)),
#              "CDF of a shifted Gompertz distributed random variable."))
    #! format: on
    export grad, hessian

    for (A, pars) ∈ pairs(As)
        fun, grad, hessian, parameters, description = pars

        ## Generate name
        Astring = string(A)
        if length(Astring) < 5 || Astring[1:5] != "Alpha"
            A = Symbol("Alpha" * Astring)
        end

        ## Generate fields
        fields = [quote
                      $parameter::Float64
                      $(Symbol(string(parameter) * "bounds"))::NTuple{2,Float64}
                  end
                  for parameter ∈ keys(parameters)]

        ## Generate docstring.
        docstring = """
$A

$description

# Parameters
    """

        for (key, val) ∈ pairs(parameters)
            docstring *= "\n- `$(string(key))::Float64`: $(first(val))"
        end

        ## Generate methods signature
        args = [:($(parameter)) for parameter ∈ keys(parameters)]
        default_bounds = [:($(last(p))) for p ∈ values(parameters)]
        complete_args = mapreduce(vcat, vcat, args, default_bounds)

        @eval begin
            export $A
            @doc $docstring
            mutable struct $A <: AbstractAlpha
                $(fields...)
            end

            # Constructors
            $A() = $A(mapreduce(x -> getindex(x, 2), vcat, $parameters)...)

            $A($(args...)) = $A($(complete_args...))

            ## Evaluation
            (α::$A)(t, $(args...)) = $fun(t, $(args...))

            function (α::$A)(t)
                parameter_values = map(par -> getfield(α, par), parameters(α))
                α(t, parameter_values...)
            end

            ## Gradient
            function grad(α::$A)
                parameter_values = map(par -> getfield(α, par), parameters(α))
                (t, $(args...)) -> $grad(t, $(args...))
            end

            ## Hessian
            function hessian(α::$A)
                parameter_values = map(par -> getfield(α, par), parameters(α))
                (t, $(args...)) -> $hessian(t, $(args...))
            end

            ## AbstractAlpha interface
            function bounds(α::$A)
                bs = map(p -> getfield(α, Symbol(string(p) * "bounds")),
                         parameters(α))
                NamedTuple{parameters(α)}(Tuple(bs))
            end

            @generated parameters(::$A) = Tuple($args)

            ## TODO: Skip type assertion
            getparameter(α::$A, parameter) = getfield(α, parameter)::Float64

            function setparameter!(α::$A, parameter, value)
                setfield!(α, parameter, convert(Float64, value))
            end
        end
    end
end

## TODO: Implement phase-type distribution.
## Phase-type distribution
mutable struct AlphaPhaseType <: AbstractAlpha
    α::Vector{Float64}
    S::Matrix{Float64}
end

function AlphaPhaseType(n)
    α = ones(n) / n

    S = zeros(n, n)
    for k ∈ 1:n
        S[k, k] = -1
    end

    AlphaPhaseType(α, S)
end

## TODO: Implement metalog distribution.
