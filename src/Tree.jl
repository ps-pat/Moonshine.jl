using Graphs

import Graphs: add_vertex!, add_edge!, rem_edge!

using Random

using Distributions

using StatsBase: AbstractWeights, FrequencyWeights, ProbabilityWeights, sample

import Base: iterate, eltype, length, size

using RandomNumbers.PCG: PCGStateOneseq

#######################
# TreeCore Definition #
#######################

struct TreeCore
    graph::SimpleDiGraph{VertexType}
    latitudes::Vector{Float64}
    sequences::Vector{Sequence}
    nleaves::Int

    positions::Vector{Float64}
    seq_length::Float64
    Ne::Float64
    μ_loc::Float64
end

function TreeCore(leaves::AbstractVector{Sequence};
                  positions = (collect ∘ range)(0, 1, length = length(leaves)),
                  seq_length = one(Float64),
                  Ne = one(Float64),
                  μ_loc = zero(Float64))
    n = length(leaves)

    sequences = similar(leaves, 2n - 1)
    for (k, leaf) ∈ enumerate(leaves)
        sequences[k] = leaf
    end

    positions = validate_positions(positions, (length ∘ first)(leaves))

    TreeCore(SimpleDiGraph(n), zeros(Float64, n - 1), sequences, n,
                positions, seq_length, Ne, μ_loc)
end

function TreeCore(rng::AbstractRNG,
                  nmin::Integer, minlength::Integer,
                  nmax::Integer = 0, maxlength::Integer = 0;
                  genpars...)
    n = iszero(nmax) ? nmin : rang(rng, nmin:nmax)
    nmarkers = iszero(maxlength) ? minlength : rang(rng, minlength:maxlength)

    TreeCore([Sequence(rng, nmarkers) for _ ∈ 1:n]; genpars...)
end

function TreeCore(nmin::Integer, minlength::Integer,
                  nmax::Integer = 0, maxlength::Integer = 0;
                  genpars...)
    TreeCore(GLOBAL_RNG, nmin, minlength, nmax, maxlength; genpars...)
end

###################
# Tree Definition #
###################

export Tree
mutable struct Tree <: AbstractGenealogy
    core::TreeCore
    logprob::BigFloat
    nextvertex::Int
end

function Tree(leaves::AbstractVector{Sequence}; genpars...)
    Tree(TreeCore(leaves; genpars...), zero(BigFloat), one(Int))
end

function Tree(rng::AbstractRNG,
              nmin::Integer, minlength::Integer,
              nmax::Integer = 0, maxlength::Integer = 0;
              genpars...)
    Tree(TreeCore(rng, nmin, minlength, nmax, maxlength; genpars...),
         zero(BigFloat), one(Int))
end

function Tree(nmin::Integer, minlength::Integer,
              nmax::Integer = 0, maxlength::Integer = 0;
              genpars...)
    Tree(GLOBAL_RNG, nmin, minlength, nmax, maxlength; genpars...)
end

##################
# Simple Methods #
##################

for field ∈ [:(:nleaves), :(:sequences), :(:latitudes), :(:graph),
             :(:seq_length), :(:Ne), :(:positions)]
    fun_name = eval(field)
    @eval begin
        export $fun_name

        function $fun_name(tree::Tree)
            getfield(tree.core, $field)
        end
    end
end

###############################
# AbstractGenealogy Interface #
###############################

describe(::Tree, long = true) = long ? "Coalescent Tree" : "Tree"

latitudes(tree::Tree, ivs) = getindex(latitudes(tree), ivs)

function latitude(tree::Tree, v)
    isleaf(tree, v) ? zero(Float64) : latitudes(tree)[v - nleaves(tree)]
end

leaves(tree::Tree) = Base.OneTo(nleaves(tree))

ivertices(tree::Tree) = range(nleaves(tree) + 1, nv(tree))

function mrca(tree::Tree)
    any(iszero.(latitudes(tree))) ?
    zero(VertexType) : argmax(latitudes(tree)) + nleaves(tree)
end

function distance(tree::Tree, v1, v2)
    2tmrca(tree, (v1, v2)) - latitude(tree, v1) - latitude(tree, v2)
end

function prob(tree::Tree; logscale = false)
    ret = tree.logprob

    logscale ? ret : exp(ret)
end

###########################
# AbstractGraph Interface #
###########################

function add_vertex!(tree::Tree, seq = Sequence(), lat = ∞)
    latitudes(tree)[tree.nextvertex] = lat
    sequences(tree)[nleaves(tree) + tree.nextvertex] = seq

    tree.nextvertex += 1

    add_vertex!(graph(tree))
end

add_edge!(tree::Tree, e) = add_edge!(graph(tree), e)

rem_edge!(tree::Tree, e) = rem_edge!(graph(tree), e)

############
# Plotting #
############

function graphplot(tree::Tree, int::Ω;
                   wild_color = :blue,
                   derived_color = :red,
                   attributes...)
    mask = fill(ancestral_mask(tree, int), nv(tree))
    node_color = ifelse.(any.(sequences(tree) .& mask),
                         derived_color, wild_color)

    invoke(graphplot, Tuple{AbstractGenealogy}, tree;
           node_color = node_color)
end

graphplot(tree::Tree; attributes...) = graphplot(tree, Ω(0); attributes...)

###########
# Methods #
###########

export dad, sibling
"""
    dad(tree, v)
    sibling(tree, v)

Parent/sibling of a vertex.
"""
function dad end,
function sibling end

for fun ∈ [:dad, :sibling]
    fun_gen = Symbol(string(fun) * 's')

    @eval function $fun(tree, v)
        (first ∘ $fun_gen)(tree, v)
    end
end

export mut_rate
mut_rate(tree::Tree, scaled = true) = tree.core.μ_loc * (scaled ? 4 * Ne(tree) : 1)

#################
# Tree Building #
#################

## TODO: get rid of idx_weights.
struct MutationSampler{T, W<:AbstractWeights}
    types::Vector{T}
    nlive::W
    transition_matrix::Matrix{Float64}

    idx::Vector{Vector{Int}}
    idx_weights::Vector{W}

    n::Int
    event_number::Base.RefValue{Int}
    logprob::Base.RefValue{BigFloat}
end

function MutationSampler(data, Dist = Hamming, args = (), kwargs = ())
    types = compute_types(data)
    ntypes = length(types)

    idx = [findall(==(type), data) for type ∈ types]
    idx_weights = [FrequencyWeights(ones(Int, length(v)), length(v)) for v ∈ idx]

    nlive = FrequencyWeights(sum.(idx_weights), length(data))

    transition_matrix = coalescence_matrix(Dist, types, args...; kwargs...)

    MutationSampler(types, nlive,
                    transition_matrix, idx, idx_weights,
                    length(data), Ref(0), Ref(big(0.0)))
end

function sample_pair!(rng, ms)
    nlive = ms.nlive
    P = ms.transition_matrix

    ## Sample a type of item.
    type_idx = sample(rng, eachindex(ms.types), nlive)
    ms.logprob[] += log(nlive[type_idx]) - log(nlive.sum)

    ## Check if there is any type greather than the one sampled with
    ## only one live item remaining which can only coalesce with
    ## sampled type. If so, make one of those coalesce with sampled
    ## type.
    can_only_coalesce_with_type = type_idx .+
        findall(range(type_idx + 1, length(ms.types))) do k
            ## Check if a coalescence is possible. For that to be the
            ## case, there has to be live items and the probability of a
            ## coalescence must be > 0.
            iszero(nlive[k] * P[type_idx, k]) && return false

            ## Check if the coalescence with type_idx is the only
            ## possibility. For that to be the case, there has to be live
            ## items other than type_idx with coalescence probability > 0.
            any(>(0),
                nlive[1:end .≠ type_idx] .* P[1:end .≠ type_idx, k]) &&
                    return false

            true
        end

    ## At that point, if there is such an item, we cancel the
    ## coalescence event and simulate a new one.
    ## TODO: make sampled type_idx taboo for next call.
    isempty(can_only_coalesce_with_type) || return sample_pair!(rng, ms)

    ## Otherwise, we're good to go!
    nlive[type_idx] -= 1
    ms.event_number[] += 1

    ## Conditionally sample a pair of indices.
    if iszero(nlive[type_idx])
        ## If there is only one item of the sampled type with a positive
        ## weight, we mutate it.

        ## Get first item.
        η1 = ms.idx[type_idx][findfirst(!iszero, ms.idx_weights[type_idx])]

        ## Sample mutation.
        type_probs = ProbabilityWeights(P[:,type_idx] .* nlive)
        newtype_idx = sample(rng, eachindex(ms.types), type_probs)
        ms.logprob[] += log(ms.transition_matrix[newtype_idx, type_idx])

        ## Get second item.
        j = sample(rng, eachindex(ms.types), ms.idx_weights[newtype_idx])
        η2 = ms.idx[newtype_idx][j]
        n = count(!iszero, ms.idx_weights[newtype_idx])
        ms.logprob[] += -log(n)

        ## Update indices.
        ms.idx[newtype_idx][j] = ms.n + ms.event_number[]
    else
        ## Sample items.
        i, j = sample(rng, eachindex(ms.idx[type_idx]), ms.idx_weights[type_idx], 2,
                      replace = false)
        η1, η2 = ms.idx[type_idx][[i, j]]

        ## Compute probability.
        n = count(!iszero, ms.idx_weights[type_idx])
        ms.logprob[] += logtwo - log(n) - log(n - 1)

        ## Update weights & indices
        ms.idx_weights[type_idx][[i, j]] .= [0, 1]
        ms.idx[type_idx][j] = ms.n + ms.event_number[]
    end

    η1, η2
end

function tree_coalesce!(rng, tree, vertices, nlive)
    ## Sample coalescing vertices. A + log(2) is missing, don't forget
    ## to account for it!
    n = length(vertices)
    tree.logprob -= log(n) + log(n - 1)
    shuffle!(rng, vertices)
    child1, child2 = pop!(vertices), pop!(vertices)

    ## Sample a latitude for the coalescence event.
    Δdist = Exponential(inv(nlive - 1))
    Δ = rand(rng, Δdist)
    tree.logprob += logpdf(Δdist, Δ)
    newlat = latitude(tree, nv(tree)) + Δ

    ## Perform the coalescence.
    _dad = coalesce!(tree, child1, child2, newlat)

    ## Update live vertices.
    push!(vertices, _dad)

    @debug("$child1 and $child2 coalesced into $_dad at $(latitude(tree, dad))")

    _dad
end

## The Coalessor
mutable struct Coalessor
    const idx::Vector{Int}
    rng::PCGStateOneseq{UInt128, Val{:XSH_RS}, UInt64}
end

Coalessor(rng, n) = Coalessor(collect(1:n), PCGStateOneseq(rand(rng, UInt)))

function iterate(iter::Coalessor, state = 1)
    n = length(iter)
    state == n + 1 && return nothing

    nlive = n - state + 2
    idx = iter.idx
    rng = iter.rng

    i = rand(rng, 1:nlive)
    x = idx[i]
    idx[i] = idx[nlive]

    j = rand(rng, 1:(nlive - 1))
    y = idx[j]
    idx[j] = n + state + 1

    (x, y), state + 1
end

eltype(::Type{Coalessor}) = NTuple{2, Int}

length(c::Coalessor) = length(c.idx) - 1

size(c::Coalessor) = (length(c),)

export build!
"""
    build!(rng, tree, idx)

Build a coalescent tree consistent with a given marker.

`tree` must only contain leaves.
"""
function build!(rng, tree::Tree, idx = 1)
    nv(tree) > nleaves(tree) && error("Tree contains non-leaf vertices")

    _sequences = view(sequences(tree), 1:nleaves(tree))

    derived = findall(σ -> getindex(σ, idx), _sequences)
    wild = setdiff(range(1, nleaves(tree)), derived)

    nlive = nleaves(tree)
    nlive_derived, nlive_wild = length(derived), length(wild)

    ## While there are more than 1 live derived vertex, we build two
    ## distinct subtrees.
    while nlive_derived > 1
        ## Determine the type of coalescence.
        wild_weight = nlive_wild * (nlive_wild - 1)
        derived_weight = nlive_derived * (nlive_derived - 1)
        derived_prob = derived_weight / (wild_weight + derived_weight)

        coalescence_dist = Bernoulli(derived_prob)
        derived_coalescence = rand(rng, coalescence_dist)
        cvertices = derived_coalescence ? derived : wild
        tree.logprob += logpdf(coalescence_dist, derived_coalescence)

        ## Perform the coalescence.
        tree_coalesce!(rng, tree, cvertices, nlive)

        ## Update the number of live nodes.
        nlive -= 1
        nlive_derived -= derived_coalescence
        nlive_wild -= !derived_coalescence
    end

    cvertices = [derived; wild]
    while nlive > 1
        tree_coalesce!(rng, tree, cvertices, nlive)
        nlive -= 1
    end

    ## Accounts for the the constant in the coalescence pprobabilities.
    tree.logprob += (nleaves(tree) - 1) * log(2)

    tree
end

build!(tree::Tree, idx = 1) = build!(GLOBAL_RNG, tree, idx)
