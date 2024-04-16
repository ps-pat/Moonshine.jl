using Random: AbstractRNG, GLOBAL_RNG

using SpecialFunctions: erf, erfc

using Distributions: BetaPrime, Chi, Chisq, cdf

using StatsFuns

using LinearAlgebra: Symmetric

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

export grad, grad!, hessian, hessian!

export @Alpha
"""
    @Alpha(A, f, ∇f, ∇²f, description, parameters...)

Generate types for simple alpha functions.

# Arguments

  - `A`: name of the function; "Alpha" will be prepended to the name of the type
    automatically
  - `f`: alpha function; its signature must be `(t, αpars...)`
  - `∇f`, `∇²f`: gradient and Hessian of `f`; details below
  - `description`: short description
  - `parameters...`: optimizable parameters; details below

# Gradient & Hessian
## One parameter

Signature of both functions must be the same as for `f`. Functions must return
the first and second derivatives of `f`.

## Multiple parameters
### Gradient

Signature must be `(g, scale, t, αpars...)`. `∇f` must mutate vector `g` by
adding the gradient evaluated at point `(t, αpars...)` multiplied by scalar
`scale` to it.

### Hessian

Signature must be `(g, scale, scale2, mult, t, αpars...)`. `∇²f` must mutate
the lower triangle of matrix `H` by adding v_ij to `H[i, j]` where

vij = `scale * f_ij(t, αpars...) + f_i(t, αpars...) * f_i(t, αpars...) * (scale2 - scale^2)`.

`f_i`, `f_j` and `f_ij` are partial derivatives of `f`.

#### JuMP compatibility warning!

`∇²f` must mutate **only the lower triangle of the matrix**. Mutating
`H[i, j]` where j > i may cause undefined behaviour.

# Parameters
Each of those arguments must have the structure

`arg = ("type", default_value, (lower_bound, upper_bound))`.

For instance, the rate parameter of the exponential alpha is provided by

`λ = ("rate", 1, (1e-6, 10))`.
"""
macro Alpha(A, f, ∇f, ∇²f, description, parameters...)
    ## Generate name
    A_string = string(A)
    if length(A_string) < 5 || A_string[1:5] != "Alpha"
        A = Symbol("Alpha" * A_string)
    end

    ## Arguments
    args = [:($(first(parameter.args))) for parameter ∈ parameters]
    nargs = length(args)
    default_bounds = [:($(last(parameter.args).args[3])) for parameter ∈ parameters]
    complete_args = mapreduce(vcat, vcat, args, default_bounds)

    ## Docstring.
    docstring_parameters = ""
    for parameter ∈ parameters
        arg = string(first(parameter.args))
        type = string(last(parameter.args).args[1])
        docstring_parameters *= "\n- `$(arg)::Float64`: $type"
    end

    ## Fields.
    fields = [quote
                  $(parameter.args[1])::Float64
                  $(Symbol(string(parameter.args[1]) * "bounds"))::NTuple{2, Float64}
              end
              for parameter ∈ parameters]

    ## Constructors
    constructors = [quote
                        function $A()
                            arguments = mapreduce(vcat, $parameters) do parameter
                                values = last(parameter.args).args
                                init = values[2]
                                bounds = Tuple(values[3].args)

                                init, bounds
                            end

                            mod = $nargs > 1 ? Iterators.flatten : identity
                            $A(mod(arguments)...)
                        end

                        $A($(args...)) = $A($(complete_args...))
                    end
                    for parameter ∈ parameters]

    esc(quote
            @doc """
    $(string($A))

$(string($description))

# Parameters
""" * $docstring_parameters
            mutable struct $A <: AbstractAlpha
                $(fields...)
            end

            ## Constructors
            $(constructors...)

            ## Evaluation
            (α::$A)(t, $(args...)) = $f(t, $(args...))

            ## TODO: make such function for gradient and hessian too.
            function (α::$A)(t)
                parameter_values = map(par -> getfield(α, par), parameters(α))
                α(t, parameter_values...)
            end

            ## Derivatives
            ## TODO: do that more elegantly, probably outside of quote.
            if $nargs > 1
                grad!(α::$A) = (g, scale, t, $(args...)) ->
                    $∇f(g, scale, t, $(args...))

                grad(α::$A) = function (t::T, $(args...)) where T
                    g = zeros(T, $nargs)
                    grad!(α)(g, 1, t, $(args...))
                    g
                end

                hessian!(α::$A) = (H, scale, scale2, mult, t, $(args...)) ->
                    $∇²f(H, scale, scale2, mult, t, $(args...))

                hessian(α::$A) = function (t::T, $(args...)) where T
                    H = zeros(T, $nargs, $nargs)
                    hessian!(α)(H, 1, 1, 1, t, $(args...))
                    Symmetric(H, :L)
                end
            else
                grad(α::$A) = (t, $(args...)) -> $∇f(t, $(args...))

                hessian(α::$A) = (t, $(args...)) -> $∇²f(t, $(args...))
            end

            ## AbstractAlpha interface
            function bounds(α::$A)
                bs = map(p -> getfield(α, Symbol(string(p) * "bounds"))::NTuple{2, Float64},
                         parameters(α))
                NamedTuple{parameters(α)}(Tuple(bs))
            end

            @generated parameters(::$A) = Tuple($args)

            getparameter(α::$A, parameter) = getfield(α, parameter)::Float64

            function setparameter!(α::$A, parameter, value)
                setfield!(α, parameter, convert(Float64, value))
            end
        end)
end

## Exponential
export AlphaExponential

const exponential = (t, λ) -> -expm1(-λ * t)
const ∇exponential = (t, λ) -> t * exp(-λ * t)
const ∇²exponential = (t, λ) -> -t^2 * exp(-λ * t)
const exponential_description = "CDF of an exponential random variable."
@Alpha(Exponential, exponential, ∇exponential, ∇²exponential,
       exponential_description,
       λ = ("rate", 1, (1e-6, 10)))

## Gompertz
export AlphaGompertz

const gompertz = (t, η, b) -> -expm1(-η * expm1(b * t))

const ∇gompertz = function (g, scale, t, η, b)
    p = b * t
    q = expm1(p)
    s = exp(-η * q)

    g[1] += q * s * scale
    g[2] += η * t * exp(p) * s * scale
    g
end

const ∇²gompertz = function (H, scale, scale2, mult, t::T, η, b) where T
    ∇α = zeros(T, 2)
    ∇gompertz(∇α, 1, t, η, b)

    p = b * t
    q = expm1(p)
    s = exp(-η * q)

    H[1, 1] += (-q^2 * s *
        scale + ∇α[1]^2 * (scale2 - scale^2)) * mult
    H[2, 2] += (η * t^2 * exp(p) * s * (1 - η * exp(b * t)) *
        scale + ∇α[2]^2 * (scale2 - scale^2)) * mult
    H[2, 1] += (t * exp(p) * s * (1 - η * q) *
        scale + ∇α[1] * ∇α[2] * (scale2 - scale^2)) * mult

    H
end

const gompertz_description = "CDF of a Gompertz distributed random variable."

@Alpha(Gompertz, gompertz, ∇gompertz, ∇²gompertz,
       gompertz_description,
       η = ("shape", 1, (1e-6, 10)), b = ("scale", 1, (1e-6, 10)))

## Maxwell-Boltzmann
export AlphaMaxwellBoltzmann

const mb = function(t, a)
    p = t / a
    erf(invsqrt2 * p) - inv(sqrthalfπ) * p * exp(-0.5 * p^2)
end

const ∇mb = (t, a) -> -inv(sqrthalfπ) * (t^3 / a^4) * exp(-0.5 * (t / a)^2)

const ∇²mb = (t, a) ->
    inv(sqrthalfπ) * (t^3 / a^7) * (2a - t) * (2a + t) * exp(-0.5 * (t / a)^2)

const mb_description = "CDF of a Maxwell-Boltzmann distributed random variable."

@Alpha(MaxwellBoltzmann, mb, ∇mb, ∇²mb,
       mb_description,
       a = ("scale", 1, (1e-6, 10)))

## Fréchet
export AlphaFrechet

const frechet = (t, a, s) -> exp(-(s / t)^a)

const ∇frechet = function (g, scale, t, a, s)
    r = s / t
    c = r^a * exp(-r^a)

    g[1] -= log(r) * c * scale
    g[2] -= (a / s) * c * scale

    g
end

const ∇²frechet = function (H, scale, scale2, mult, t::T, a, s) where T
    ∇α = zeros(T, 2)
    ∇frechet(∇α, 1, t, a, s)

    r = s / t
    c = r^a * exp(-r^a)
    r1 = r^a - 1
    l = log(r)
    q = r1 - t

    H[1, 1] += (l^2 * r1 * c *
        scale + ∇α[1]^2 * (scale2 - scale^2)) * mult
    H[2, 2] += ((a / s^2) * (a * r1 + 1) * c *
        scale + ∇α[2]^2 * (scale2 - scale^2)) * mult
    H[2, 1] += ((a * r1 * l - 1) / s * c *
        scale + ∇α[1] * ∇α[2] * (scale2 - scale^2)) * mult
    H
end

const frechet_description = "CDF of a Frechet distributed random variable."

@Alpha(Frechet, frechet, ∇frechet, ∇²frechet,
       frechet_description,
       a = ("shape", 1, (1e-6, 10)), s = ("scale", 1, (1e-6, 10)))

## Lomax

export AlphaLomax

const lomax = (t, a, λ) -> 1 - inv(1 + t / λ)^a

const ∇lomax = function (g, scale, t, a, λ)
    q = 1 + t / λ
    qa = inv(q)^a

    g[1] += qa * log(q) * scale
    g[2] -= a * t * qa / (λ * (λ + t)) * scale
    g
end

const ∇²lomax = function(H, scale, scale2, mult, t::T, a, λ) where T
    ∇α = zeros(T, 2)
    ∇lomax(∇α, 1, t, a, λ)

    q = 1 + t / λ
    logq = log(q)
    qa = inv(q)^a
    λt = λ * (λ + t)

    H[1, 1] += (-qa * logq^2 *
        scale + ∇α[1]^2 * (scale2 - scale^2)) * mult
    H[2, 2] += (qa * a * t * (2λ + (1 - a)t) / λt^2 *
        scale + ∇α[2]^2 * (scale2 - scale^2)) * mult
    H[2, 1] += (qa * t * (logq * a - 1) / λt *
        scale + ∇α[1] * ∇α[2] * (scale2 - scale^2)) * mult
    H
end

const lomax_description = "CDF of a Lomax distributed random variable."

@Alpha(Lomax, lomax, ∇lomax, ∇²lomax,
       lomax_description,
       a = ("shape", 1, (1e-6, 1e6)), λ = ("scale", 1, (1e-6, 1e6)))
