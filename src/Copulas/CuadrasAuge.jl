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
