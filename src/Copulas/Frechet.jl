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

function deriv_frechet(copula::CopulaFrechet{<:Bernoulli}, deriv_fun)
    α = alpha(copula)
    ∇α = deriv_fun(α)

    function (φ, ψ, t, αpars...)
        num = ∇α(t, αpars...)
        denum = α(t, αpars...)
        φ ⊻ ψ && return num ./ denum

        if φ
            mult = failprob(marginal(copula))
        else
            mult = succprob(marginal(copula))
        end

        (mult .* num) ./ ( mult .* denum .- 1)
    end
end

∇logpdf_joint(copula::CopulaFrechet{<:Bernoulli}) = deriv_frechet(copula, grad)

∇²logpdf_joint(copula::CopulaFrechet{<:Bernoulli}) = function (φ, ψ, t, αpars...)
    deriv_frechet(copula, hessian)(φ, ψ, t, αpars...) .-
        (∇logpdf_joint(copula)(φ, ψ, t, αpars...)).^2
end
