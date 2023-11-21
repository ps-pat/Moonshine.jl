using Random: GLOBAL_RNG

export AlphaExponential
mutable struct AlphaExponential<:AbstractAlpha
    λ::Float64
    const bounds::NTuple{2, Float64}
end

AlphaExponential() = AlphaExponential(1)

AlphaExponential(λ) = AlphaExponential(λ, eps(Float64), 100)

AlphaExponential(λ, lbound, ubound) = AlphaExponential(λ, (lbound, ubound))

(α::AlphaExponential)(t) = α(t, α.λ)
(α::AlphaExponential)(t, λ) = -expm1(-λ * t)

###########################
# AbstractAlpha Interface #
###########################

bounds(α::AlphaExponential) = (;λ = α.bounds)

@generated parameters(::AlphaExponential) = (:λ,)

getparameter(α::AlphaExponential, parameter) = getfield(α, parameter)

setparameter!(α::AlphaExponential, parameter, value) =
    setfield!(α, parameter, value)
