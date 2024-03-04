@copula_struct CuadrasAuge

####################
# Binary Phenotype #
####################

function pdf_joint(copula::CopulaCuadrasAuge{<:Bernoulli})
    α = alpha(copula)
    q = failprob(marginal(copula))

    function (φ, ψ, t, αpars...)
        αt = α(t, αpars...)

        ret = q^αt

        s = φ + ψ

        if !iszero(s)
            ret = s - ret
        end

        ret *= q

        if s == 2
            ret = 1 - ret
        end

        ret

    end
end

function pdf_conditional(copula::CopulaCuadrasAuge{<:Bernoulli})
    α = alpha(copula)
    p = succprob(marginal(copula))

    function (φ, ψ, t, αpars...)
        αt = α(t, αpars...)

        ret = (1 - p)^αt

        s = φ + ψ

        if !iszero(s)
            ret = s - ret
        end

        if isone(ψ)
            ret *= (1 - p) / p
        end

        if s == 2
            ret = 1 / p - ret

        end

        ret
    end
end

function ∇logpdf_joint(copula::CopulaCuadrasAuge{<:Bernoulli})
    α = alpha(copula)
    ∇α = grad(α)
    q = failprob(marginal(copula))

    function (φ, ψ, t, αpars...)
        f00 = log(q) .* ∇α(t, αpars...)

        !(φ || ψ) && return f00

        qα = inv(q^(α(t, αpars...)))
        if φ ⊻ ψ
            d = 1 - qα
        else
            d = qα * (inv(q) - 2) + 1
        end

        f00 ./ d
    end
end
