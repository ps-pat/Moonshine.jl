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
struct IsChain{G, H, P}
    fH::H
    fΦ::P

    H::Vector{Sequence}

    genealogies::Vector{G}

    seed::UInt128
end

function IsChain{G}(fH, fΦ, haplotypes::AbstractVector{Sequence}, n, seed) where G
    H, P = typeof(fH), typeof(fΦ)
    IsChain{G,H,P}(fH, fΦ, haplotypes, Vector{G}(undef, n), seed)
end

IsChain{G}(fH, fΦ, haplotypes::AbstractVector{Sequence}, n) where G =
    IsChain{G}(fH, fΦ, haplotypes::AbstractVector{Sequence}, n, rand(UInt128))

#######################
# Iteration interface #
#######################

function iterate(chain::IsChain, state = 1)
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
weights(chain, logscale = false) = Iterators.map(chain) do genealogy
    fG = chain.fH
    log_weight = fG(genealogy, logscale = true) - prob(genealogy, logscale = true)
    logscale ? log_weight : exp(log_weight)
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

export sample_dist
function sample_dist(chain::IsChain{G,H,P}) where {G,H,P<:PhenotypeBernoulli}
    fΦ = chain.fΦ
    fG = chain.fH
    nmissing = sum(ismissing.(phenotypes(chain)))

    res = zeros(BigFloat, Int(exp2(nmissing)))

    Threads.@threads for genealogy ∈ chain
        log_fΦ = fΦ(genealogy)
        log_w =  fG(genealogy, logscale = true) - prob(genealogy, logscale = true)

        res .+= exp.(log_fΦ .+ log_w)
    end

    BernoulliMulti(res ./ sum(res, dims = 1))
end
