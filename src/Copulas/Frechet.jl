@copula_struct Frechet

####################
# Binary Phenotype #
####################

function logpdf_joint(copula::CopulaFrechet{<:Bernoulli})
    α = alpha(copula)
    p = succprob(marginal(copula))

    function (φ, ψ, t, αpars...)
        αt = α(t, αpars...)

        if isone(ψ)
            t1 = log(p)
            s2 = (1 - p) * αt
        else
            t1 = log(1 - p)
            s2 = p * αt
        end

        t2 = log((φ ⊻ ψ) ? s2 : 1 - s2)

        t1 + t2
    end
end

function pdf_conditional(copula::CopulaFrechet{<:Bernoulli})
    α = alpha(copula)
    p = succprob(marginal(copula))

    function (φ, ψ, t, αpars...)
        αt = α(t, αpars...)

        if isone(ψ)
            s2 = (1 - p) * αt
        else
            s2 = p * αt
        end

        φ ⊻ ψ ? s2 : 1 - s2
    end
end

∇scale(copula::CopulaFrechet{<:Bernoulli}) = function (φ, ψ, t::T, αpars...) where T
    α = alpha(copula)(t, αpars...)

    φ ⊻ ψ && return inv(α)

    pq = (φ ? failprob : succprob)(marginal(copula))
    inv(α - inv(pq))
end

∇²scale(copula::CopulaFrechet{<:Bernoulli}) = function (φ, ψ, t::T, αpars...) where T
    zero(T)
end
