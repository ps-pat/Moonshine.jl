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

include("basic_alphas.jl")
