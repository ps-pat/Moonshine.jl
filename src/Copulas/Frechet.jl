export FrechetCopula

struct FrechetCopula{P, A, N} <: AbstractΦCopula{P, A, N}
    α::A
    parameters::NTuple{N, Float64}

    FrechetCopula(P, α, pars...) = new{P, typeof(α), npars(P)}(α, pars)
end

####################
# Binary Phenotype #
####################

## Utilities

bernoulli_cdf(x, p) = (0 ≤ x < 1) * (1 - p) + (x ≥ 1)
bernoulli_cdf(φ::Bool, p) = φ ? 1 : 1 - p

## AbstractΦCopula Interface

function cdf(copula::FrechetCopula{PhenotypeBinary}, φ1, φ2, d, αpars...)
    α = copula.α(d, αpars...)
    p = first(copula.parameters)

    F1 = bernoulli_cdf(φ1, p)
    F2 = bernoulli_cdf(φ2, p)

    α * F1 * F2 + (1 - α) * min(F1, F2)
end

function logpdf(copula::FrechetCopula{PhenotypeBinary}, φ1, φ2, d, αpars...)
    α = copula.α(d, αpars...)
    p = first(copula.parameters)

    if φ2
        t1 = log(p)
        s2 = (1 - p) * α
    else
        t1 = log(1 - p)
        s2 = p * α
    end

    t2 = log((φ1 ⊻ φ2) ? s2 : 1 - s2)

    t1 + t2
end

function conditional_pdf(copula::FrechetCopula{PhenotypeBinary}, φ1, φ2, d, αpars...)
    α = copula.α(d, αpars...)
    αtilds = α * (2 - α)
    p = first(copula.parameters)

    if φ2
        s2 = (1 - p) * αtilde
    else
        s2 = p * αtilde
    end

    φ1 ⊻ φ2 ? s2 : 1 - s2
end
