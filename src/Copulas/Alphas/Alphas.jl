using Random: AbstractRNG, GLOBAL_RNG

## TODO: Generalize to ARGs.

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

###################
# Packaged Alphas #
###################

include("Exponential.jl")
