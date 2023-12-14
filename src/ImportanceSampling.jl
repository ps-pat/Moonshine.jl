using Statistics: mean

import Base:
    ## Iteration interface.
    iterate,
    eltype,
    length,
    size,
    ## Indexing interface.
    firstindex,
    lastindex,
    getindex,
    setindex!

using Base.Threads: nthreads, @threads, threadid
using RandomNumbers.PCG: PCGStateSetseq

#########################
# Definition of IsChain #
#########################

export IsChain
struct IsChain{G, H, P, T}
    fH::H
    fΦ::P

    H::Vector{Sequence{T}}

    genealogies::Vector{G}

    seed::UInt128
end

function IsChain{G}(fH, fΦ, haplotypes::AbstractVector{Sequence{T}}, n, seed) where {G,T}
    H, P = typeof(fH), typeof(fΦ)
    IsChain{G, H,P,T}(fH, fΦ, haplotypes, Vector{G}(undef, n), seed)
end

IsChain{G}(fH, fΦ, haplotypes::AbstractVector{Sequence{T}}, n) where {G,T} =
    IsChain{G}(fH, fΦ, haplotypes::AbstractVector{Sequence{T}}, n, rand(UInt128))

#######################
# Iteration interface #
#######################

function iterate(chain::IsChain, state = 1)
    # state > lastindex(chain) ? nothing : (chain[state], state + 1)
    while state <= lastindex(chain)
        isassigned(chain.genealogies, state) && return (chain[state], state + 1)
        state += 1
    end

    nothing
end

@generated eltype(::Type{IsChain{G}}) where G = G

@generated eltype(chain::IsChain{G}) where G = eltype(chain)

length(chain::IsChain) = length(chain.genealogies)

size(chain::IsChain) = (length(chain),)

## Indexing interface.
lastindex(chain::IsChain) = length(chain)

firstindex(chain::IsChain) = lastindex(chain) > 0 ? 1 : 0

getindex(chain::IsChain, i) = getindex(chain.genealogies, i)

## No setindex!

###########
# Methods #
###########

phenotypes(ic) = ic.fΦ.Φ

export simulate!
"""
    simulate!(chain[, idx])

Simulate genealogies for a chain. Indices of the genealogies to be simulated
cand be specified via the optional argument `idx`.
"""
function simulate! end

function simulate!(chain::IsChain{G}, idx) where G
    Threads.@threads for k ∈ idx
        rng = PCGStateSetseq((chain.seed, k))
        chain.genealogies[k] = G(chain.H; genpars(chain.fH)...)
        build!(rng, chain.genealogies[k])
    end

    chain
end

simulate!(chain) = simulate!(chain, eachindex(chain.genealogies))

export weights
"""
    weights(chain; logscale = false)

Importance samping weights of the chain.
"""
function weights(chain; logscale = false)
    fG = chain.fH

    GetLogPrs(f) = Map(g -> f(g, logscale = true))

    chain |>
        Zip(GetLogPrs(prob), GetLogPrs(fG)) ⨟ MapSplat(-) ⨟ Map(logscale ? identity : exp)
end

export ess
"""
    ess(chain)

Effective sample size of the chain.
"""
function ess(chain)
    ws = weights(chain)
    sum(ws)^2 / sum(w -> w^2, ws)
end

macro crazy_dispatch(fun, dist)
    :($fun(ic::IsChain))
end

##########################
# Density reconstruction #
##########################

let dists = (:Bernoulli,)
    for dist ∈ dists
        name = Symbol("Phenotype" * string(dist))
        @eval begin
            const $name = PhenotypeDensity{T,C} where {D<:$dist, T, C<:AbstractΦCopula{D}}
        end
    end
end

function reduce_dist(part, nmissing)
    res = zeros(BigFloat, Int(exp2(nmissing)))

    for (log_fΦ, log_w) ∈ part
        res .+= exp.(log_fΦ .+ log_w)
    end

    res
end

function raw_sample_pdf(chain::IsChain{G,H,P}) where {G,H,P<:PhenotypeBernoulli}
    fΦ = chain.fΦ
    fG = chain.fH
    nmissing = sum(ismissing.(phenotypes(chain)))

    GetLogPrs(f) = Map(g -> f(g, logscale = true))
    log_fΦs = Map(genealogy -> log.(fΦ(genealogy)))
    log_ws = Zip(GetLogPrs(prob), GetLogPrs(fG)) ⨟ MapSplat(.-)


    ## TODO: parallelize.
    res = zeros(BigFloat, Int(exp2(nmissing)))

    for (log_fΦ, log_w) ∈ chain |> Zip(log_fΦs, log_ws)
        res .+= exp.(log_fΦ .+ log_w)
    end

    res ./ sum(res, dims = 1)
end

export sample_pdf
function sample_pdf(chain::IsChain{G,H,P}) where {G,H,P<:PhenotypeBernoulli}
    p = raw_sample_pdf(chain)

    length(p) > 1 ? BernoulliMulti(p) : Bernoulli(p)
end
