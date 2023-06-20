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

## Ugly circular reference hack, sorry!

export IsChain
struct IsChain{H<:AbstractGraphDensity, P<:AbstractGraphDensity, T, R}
    dens_haps::H
    dens_phenos::P

    haplotypes::Vector{Sequence{T}}

    records::Vector{R}
end

IsChain(dens_haps, dens_phenos, haplotypes, n) =
    IsChain(dens_haps, dens_phenos, haplotypes, Vector{IsRecord}(undef, n))

## Iteration interface.
iterate(chain::IsChain, state = 1) =
    state > lastindex(chain) ? nothing : (chain[state], state + 1)

@generated eltype(::Type{IsChain{H, P, T, R}}) where{H, P, T, R} = R

@generated eltype(::IsChain) = eltype(IsChain)

length(chain::IsChain) = length(chain.records)

size(chain::IsChain) = (length(chain),)

## Indexing interface.
lastindex(chain::IsChain) = length(chain)

firstindex(chain::IsChain) = lastindex(chain) > 0 ? 1 : 0

getindex(chain::IsChain, i) = getindex(chain.records, i)

setindex!(chain::IsChain, x, i) = setindex!(chain.records, x, i)

export joint
"""
    joint(chain, phenotypes = nothing; logscale = false)

Compute the joint density of a pair of haplotypes/phenotypes vector. If
`phenotypes` is provided, it must be the same length as the number of `missing`
phenotypes in `chain`. The returned value corresponds to the "complete data"
joint density. Otherwise, only the phenotypes in `chain` are accounted for;
this corresponds to the "incomplete data" joint density.
"""
function joint(chain, phenotypes = nothing; logscale = false)
    dens_tree = chain.dens_haps

    ret = mean(chain) do record
        dens_φ = record.dens_pheno
        arg = record.arg

        weight_log = dens_tree(arg, logscale = true) - arg.logprob
        pdf_log = dens_φ(phenotypes, logscale = true)

        exp(weight_log + pdf_log)
    end

    logscale ? log(ret) : ret
end

export likelihood
"""
    likelihood(chain, phenotypes; logscale = false)

Compute the complete/incomplete data joint densities ratio. This corresponds to
the likelihood of `phenotypes`.
"""
function likelihood(chain, phenotypes; logscale = false)
    ret_log = joint(chain, phenotypes, logscale = true) -
        joint(chain, logscale = true)

    logscale ? ret_log : exp(ret_log)
end

##########################
# Definition of IsRecord #
##########################

struct IsRecord
    parent_chain::IsChain{H, P, T, IsRecord} where {H, P, T}

    arg::Arg
    dens_pheno::Function
end

"""
    compute_joint([rng], haplotypes, phenotypes, dens_haps, dens_phenos;
                  n_is = 1000, n_mcmc = 1000)

Compute the joint density of a vector of haplotypes and
phenotypes. Returns an `IsChain`.

# Arguments
- `haplotypes::Vector`: vector of haplotypes
- `dens_haps::AbstractGraphDensity`: density for `haplotypes`
- `dens_phenos::AbstractGraphDensity`: density for phenotypes conditional on
  `haplotypes`
- `n_is`: number of is samples to draw (for integration wrt arg)
- `n_mcmc`: number of MCMC samples to draw (for integration wrt phenotypes)
"""
function compute_joint end
export compute_joint

function compute_joint(rng, haplotypes, dens_haps, dens_phenos;
                       n_is = 1000, n_mcmc = 1000)
    chain = IsChain(dens_haps, dens_phenos, haplotypes, n_is)

    ## Initialize the rng of each worker.
    seed = rand(rng, Int)

    @threads for k ∈ 1:n_is
        rng_local = PCGStateSetseq((seed, k))

        tree = buildtree!(rng_local, Arg(haplotypes))
        dens_φ = chain.dens_phenos(rng_local, tree, M = n_mcmc)

        chain[k] = IsRecord(chain, tree, dens_φ)
    end

    chain
end
