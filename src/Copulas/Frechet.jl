@copula_struct Frechet

####################
# Binary Phenotype #
####################

function pdf_joint(copula::CopulaFrechet{<:Bernoulli}, φ, ψ, t, αpars;
                   logscale = false)
    α = alpha(copula)
    p = succprob(marginal(copula))

    αt = α(t, αpars)

    if isone(φ)
        t1 = log(p)
        s2 = (1 - p) * αt
    else
        t1 = log(1 - p)
        s2 = p * αt
    end

    t2 = log((φ ⊻ ψ) ? s2 : 1 - s2)

    logscale ? t1 + t2 : exp(t1) * exp(t2)
end

function pdf_conditional(copula::CopulaFrechet{<:Bernoulli}, φ, ψ, t;
                         logscale = false)
    α = alpha(copula)
    p = succprob(marginal(copula))

    αt = α(t)

    if isone(ψ)
        s2 = (1 - p) * α(t)
    else
        s2 = p * α(t)
    end

    t2 = φ ⊻ ψ ? s2 : 1 - s2
    logscale ? log(t2) : t2
end

∇scale(copula::CopulaFrechet{<:Bernoulli}) = function (φ, ψ, t::T, αpars) where T
    α = alpha(copula)(t, αpars)

    φ ⊻ ψ && return inv(α)

    pq = (φ ? failprob : succprob)(marginal(copula))
    inv(α - inv(pq))
end

∇²scale(copula::CopulaFrechet{<:Bernoulli}) = function (φ, ψ, t::T, αpars...) where T
    zero(T)
end
