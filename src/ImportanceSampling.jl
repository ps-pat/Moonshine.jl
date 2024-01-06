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
    setindex!,
    view

using Base.Threads
using RandomNumbers.PCG: PCGStateSetseq

using ChunkSplitters
import ChunkSplitters: chunks

#########################
# Definition of IsChain #
#########################

export IsChain
struct IsChain{G, H}
    fG::H

    H::Vector{Sequence}

    genealogies::Vector{G}

    seed::UInt128
end

function IsChain{G}(fG, haplotypes::AbstractVector{Sequence}, n, seed) where G
    H = typeof(fG)
    IsChain{G,H}(fG, haplotypes, Vector{G}(undef, n), seed)
end

IsChain{G}(fG, haplotypes::AbstractVector{Sequence}, n) where G =
    IsChain{G}(fG, haplotypes::AbstractVector{Sequence}, n, rand(UInt128))

## ChunkSplitters
chunks(x::IsChain, nchunks::Int, type = :batch) =
    chunks(x.genealogies, nchunks, type)

view(chain::IsChain, idx) = view(chain.genealogies, idx)

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
        chain.genealogies[k] = G(chain.H; genpars(chain.fG)...)
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
    fG = chain.fG
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
function sample_dist(chain::IsChain{G,H}, fΦ::P) where {G,H,P<:PhenotypeBernoulli}
    fG = chain.fG
    nmissing = sum(ismissing.(phenotypes(fΦ)))

    res_parts = map(chunks(chain, Threads.nthreads())) do (rg, _)
        @spawn begin
            local acc = zeros(BigFloat, Int(exp2(nmissing)))

            for genealogy ∈ view(chain, rg)
                local log_fΦ = log.(fΦ(genealogy))
                local log_w =  fG(genealogy, logscale = true) - prob(genealogy, logscale = true)

                acc .+= exp.(log_fΦ .+ log_w)
            end

            acc
        end
    end

    res = sum(fetch, res_parts)

    BernoulliMulti(res ./ sum(res, dims = 1))
end
