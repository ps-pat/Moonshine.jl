@copula_struct CuadrasAuge

####################
# Binary Phenotype #
####################

function pdf_joint(copula::CopulaCuadrasAuge{<:Bernoulli}, φ, ψ, t, αpars;
                   logscale = false)
    α = alpha(copula)(t, αpars)
    q = failprob(marginal(copula))

    ret = q^α

    s = φ + ψ

    if !iszero(s)
        ret *= -1
        ret += s
    end

    ret *= q

    if s > 1
        ret *= -1
        ret += 1
    end

    logscale ? log(ret) : ret
end

function pdf_conditional(copula::CopulaCuadrasAuge{<:Bernoulli}, φ, ψ, t;
                         logscale = false)
    α = alpha(copula)(t)
    q = failprob(marginal(copula))

    s = φ + ψ

    ret = q^α - s

    if isone(s)
        if φ
            ret *= -1
        else
            ret /= 1 - inv(q)
        end
    elseif s > 1
        ret *= q
        ret += 1
        ret *= inv(1 - q)
    end

    logscale ? log(ret) : ret
end

∇scale(copula::CopulaCuadrasAuge{<:Bernoulli}) = function (φ, ψ, t, αpars)
    α = alpha(copula)(t, αpars)
    q = failprob(marginal(copula))
    qα = q^α

    s = φ + ψ

    div = q * (qα - s)

    if s > 1
        div += 1
        div /= 3
    end

    α * qα / div
end

∇²scale(copula::CopulaCuadrasAuge{<:Bernoulli}) = function (φ, ψ, t, αpars)
    α = alpha(copula)(t, αpars)
    q = failprob(marginal(copula))
    qα = q^α

    s = φ + ψ

    div = q

    if !iszero(s)
        div *= qα - s
        div /= qα
    end

    if s > 1
        div += inv(qα)
        div /= 3
    end

    (1 + α^2 / q) / div
end
