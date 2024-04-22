__precompile__(false)

###############
# Exponential #
###############

export AlphaExponential

const exponential = (t, λ) -> -expm1(-λ * t)
const ∇exponential = (t, λ) -> t * exp(-λ * t)
const ∇²exponential = (t, λ) -> -t^2 * exp(-λ * t)
const exponential_description = "CDF of an exponential random variable."
@Alpha(Exponential, exponential, ∇exponential, ∇²exponential,
       exponential_description,
       λ = ("rate", 1, (1e-6, 10)))

#####################
# Maxwell-Boltzmann #
#####################

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

############
# Gompertz #
############

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

###########
# Fréchet #
###########

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

#########
# Lomax #
#########

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
