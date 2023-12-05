export CopulaCuadrasAuge

struct CopulaCuadrasAuge{P,A,N} <: AbstractΦCopula{P,A,N}
    α::A
    parameters::NTuple{N,Float64}

    CopulaCuadrasAuge(P, α, pars...) = new{P,typeof(α),npars(P)}(α, pars)
end

####################
# Binary Phenotype #
####################

## AbstractΦCopula Interface

function pdf(copula::CopulaCuadrasAuge{PhenotypeBinary}, φ, ψ, t, αpars...)
    α = copula.α(t, αpars...)
    q = 1 - first(copula.parameters)
    s = φ + ψ

    ret = q^α

    if !iszero(s)
        ret = s - ret
    end

    ret *= q

    if s == 2
        ret = 1 - ret
    end

    ret
end

function conditional_pdf(copula::CopulaCuadrasAuge{PhenotypeBinary}, φ, ψ, t, αpars...)
α = copula.α(t, αpars...)
    q = 1 - first(copula.parameters)
    s = φ + ψ

    ret = q^α

    if !iszero(s)
        ret = s - ret
    end

    if isone(ψ)
        ret *= q / (1 - q)
    end

    if s == 2
        ret = 1 / q - ret
    end

    ret
end
