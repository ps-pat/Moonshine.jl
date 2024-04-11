using Graphs

import Graphs: add_vertex!, add_edge!, rem_edge!, ne

using Random

using Distributions

using StatsBase: AbstractWeights, FrequencyWeights, ProbabilityWeights, sample

import Base: iterate, length, size, isempty

using RandomNumbers.PCG: PCGStateOneseq

#######################
# TreeCore Definition #
#######################

struct TreeCore
    graph::SimpleDiGraph{VertexType}
    latitudes::Vector{Float64}
    sequences::Vector{Sequence}

    positions::Vector{Float64}
    seq_length::Float64
    Ne::Float64
    μ_loc::Float64
end

TreeCore() = TreeCore(SimpleDiGraph{VertexType}(), [], [], [],
                      typemin(Float64), typemin(Float64), typemin(Float64))

function TreeCore(leaves::AbstractVector{Sequence};
                  positions = nothing,
                  seq_length = one(Float64),
                  Ne = one(Float64),
                  μ_loc = zero(Float64))
    n = length(leaves)

    sequences = similar(leaves, 2n - 1)
    for (k, leaf) ∈ enumerate(leaves)
        sequences[k] = leaf
    end

    if isnothing(positions)
        positions =
            isone(n) ? zeros(1) : (collect ∘ range)(0, 1, length = length(leaves))
    else
        positions = validate_positions(positions, (length ∘ first)(leaves))
    end

    TreeCore(SimpleDiGraph(n), zeros(Float64, n - 1), sequences,
                positions, seq_length, Ne, μ_loc)
end

function TreeCore(rng::AbstractRNG,
                  nmin::Integer, minlength::Integer,
                  nmax::Integer = 0, maxlength::Integer = 0;
                  genpars...)
    n = iszero(nmax) ? nmin : rand(rng, nmin:nmax)
    nmarkers = iszero(maxlength) ? minlength : rand(rng, minlength:maxlength)

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

Tree() = Tree(TreeCore(), typemin(BigFloat), typemin(Int))

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

for field ∈ [:(:sequences), :(:latitudes), :(:graph),
             :(:seq_length), :(:Ne), :(:positions)]
    fun_name = eval(field)
    @eval begin
        export $fun_name

        function $fun_name(tree::Tree)
            getfield(tree.core, $field)
        end
    end
end

export nleaves
nleaves(tree::Tree) = (length(sequences(tree)) + 1) ÷ 2

isempty(treecore::TreeCore) = isempty(treecore.sequences)

isempty(tree::Tree) = isempty(tree.core)

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
    any(iszero, latitudes(tree)) && return zero(VertexType)
    isone(nv(tree)) && return one(VertexType)
    argmax(latitudes(tree)) + nleaves(tree)
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

function add_vertex!(tree::Tree, seq, lat)
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

export depth
"""
    depth(tree, v)

Depth of a vertex.
"""
function depth(tree, v)
    d = zero(Int)
    _mrca = mrca(tree)

    while v ≠ _mrca
        v = dad(tree, v)
        d += 1
    end

    d
end

export leaves_permutation
"""
    leaves_permutation(tree)
    leaves_permutation(tree, vs)

Permutation of (a subset of ) 1:n induced by a tree.
"""
function leaves_permutation(tree, vs)
    sort!(vs)
    adepths = Dict{Int, Vector{eltype(tree)}}()

    for v ∈ vs
        d = depth(tree, v)
        if !haskey(adepths, d)
            adepths[d] = Vector{eltype(tree)}()
        end
        push!(adepths[d], v)
    end

    depths = (sort ∘ collect)(keys(adepths))

    ret = similar(vs)
    retptr = 1
    @inbounds for d ∈ depths
        for v ∈ adepths[d]
            ret[retptr] = v
            retptr += 1
        end
    end

    ret
end

leaves_permutation(tree) = leaves_permutation(tree, leaves(tree))

#################
# Tree Building #
#################

struct MutationSampler{T, W<:AbstractWeights}
    types::Vector{T}
    nlive::W
    transition_matrix::Matrix{Float64}

    idx::Vector{Vector{Int}}

    n::Int
    event_number::Base.RefValue{Int}
    logprob::Base.RefValue{BigFloat}
end

function MutationSampler(data, Dist = Hamming{Float64}, args = (), kwargs = ())
    types = compute_types(data)
    ntypes = length(types)

    idx = [findall(==(type), data) for type ∈ types]

    nlive = FrequencyWeights(length.(idx), length(data))

    transition_matrix = coalescence_matrix(Dist, types, args...; kwargs...)

    MutationSampler(types, nlive,
                    transition_matrix, idx,
                    length(data), Ref(0), Ref(big(0.0)))
end

function sample_pair!(rng, ms)
    ms.event_number[] += 1

    nlive = ms.nlive
    P = ms.transition_matrix

    ## Sample a type of item.
    type_idx = something(sample(rng, eachindex(ms.types), nlive))
    ms.logprob[] += log(nlive[type_idx]) - log(nlive.sum)

    ## Check if there is any type greater than the one sampled with
    ## only one live item remaining which can only coalesce with
    ## sampled type. If any is found, switch type_idx for the greatest
    ## element with at least 1 live item.
    for i ∈ range(type_idx + 1, length(ms.types))
        ## Check if a coalescence is possible. For that to be the
        ## case, there has to be live items and the probability of a
        ## coalescence must be > 0.
        iszero(nlive[i]) || iszero(P[type_idx, i]) && continue

        ## Check if the coalescence with type_idx is the only
        ## possibility. For that to be the case, there has to be live
        ## items other than type_idx with coalescence probability > 0.
        for j ∈ eachindex(nlive)
            j == type_idx && continue
            !iszero(nlive[j]) && !iszero(P[j, i]) && @goto bypass
        end

        ## If we reached that point (and didn't jump straight to the
        ## bypass), there is such an item.
        type_idx = something(findlast(!iszero, nlive))
        break

        @label bypass
    end

    ## Conditionally sample a pair of indices.
    if isone(nlive[type_idx])
        ## If there is only one item of the sampled type with a positive
        ## weight, we mutate it.

        ## Get first item.
        η1 = first(ms.idx[type_idx])
        nlive[type_idx] -= 1

        ## Sample mutation.
        type_probs = ProbabilityWeights(P[:,type_idx] .* nlive.values)
        newtype_idx = sample(rng, eachindex(ms.types), type_probs)
        ms.logprob[] += log(ms.transition_matrix[newtype_idx, type_idx])

        ## Sample second item.
        j = rand(rng, 1:nlive[newtype_idx])
        η2 = ms.idx[newtype_idx][j]

        ## Update indices.
        ms.idx[newtype_idx][j] = ms.n + ms.event_number[]

        ## Update probability.
        ms.logprob[] -= log(nlive[newtype_idx])
    else
        idx = ms.idx[type_idx]
        n = nlive[type_idx]

        ## Sample items. i and j are the indices of the first and
        ## second sampled items respectively. η1 and η2 are the actual
        ## items. After this iteration:
        ##   1. nlive[type_idx] will be decreased by 1;
        ##   2. idx[i] will contain the item that was previously at index
        ##      nlive[type_idx];
        ##   3. idx[j] will contain the item created by the coalescence of
        ##      η1 with η2.

        i = rand(rng, 1:n)
        η1 = idx[i]
        idx[i] = idx[n]
        nlive[type_idx] -= 1

        j = rand(rng, 1:(n - 1))
        η2 = idx[j]
        idx[j] = ms.n + ms.event_number[]

        ## Update probability.
        ms.logprob[] += logtwo - log(n) - log(n - 1)
    end

    η1, η2
end

export build!
"""
    build!(rng, tree)

Build a coalescent tree consistent with a given marker.

`tree` must only contain leaves.
"""
function build! end

function build!(rng, tree::Tree)
    n = nleaves(tree)
    nv(tree) > nleaves(tree) && error("Tree contains non-leaf vertices")

    ## Create internal vertices.
    add_vertices!(graph(tree), n - 1)

    ## Sample coalescence vertices latitudes.
    Δlatitudes = randexp!(rng, latitudes(tree)) ./ range(n - 1, 1, step = -1)
    cumsum!(latitudes(tree), Δlatitudes)
    tree.logprob += sum((enumerate ∘ reverse)(big.(Δlatitudes)),
                        init = zero(BigFloat)) do (λ, t)
        log(-expm1(-λ * t))
    end

    ## Sample coalescence events
    ms = MutationSampler((collect ∘ sequences)(tree, leaves(tree)))
    for parent ∈ range(n + 1, length = n - 1)
        child1, child2 = sample_pair!(rng, ms)
        add_edge!(tree, parent, child1)
        add_edge!(tree, parent, child2)

        tree.core.sequences[parent] =
            tree.core.sequences[child1] & tree.core.sequences[child2]
    end

    tree.logprob += ms.logprob[]

    tree
end

build!(tree::Tree) = build!(GLOBAL_RNG, tree)

function isvalid(tree::Tree)
    n = nleaves(tree)

    ## Check that the latitudes of the internal vertices is increasing.
    let (last_latitude, latitudes_iterator) = Iterators.peel(latitudes(tree))
        for (k, latitude) ∈ enumerate(latitudes_iterator)
            v = n + k + 1
            if latitude < last_latitude
                @info "Non increasing latitude" vertex = v
                return false
            end
            last_latitude = latitude
        end
    end

    ## Check that the children of every internal vertices are lesser
    ## than itself.
    for v ∈ ivertices(tree)
        _children = children(tree, v)
        all(<(v), _children) && continue

        @info "Invalid ordering" vertex = v children = _children
        return false
    end

    ## Check the number of children and number of parents of each
    ## vertex.

    ## The root must have 2 children and no parent.
    let r = mrca(tree)
        if !isempty(dads(tree, r)) || length(children(tree, r)) ≠ 2
            @info "Invalid tree mrca" mrca = r parents = dads(tree, r) children = children(tree, r)
            return false
        end
    end

    ## Non root internal vertices must have 2 children and 1 parent.
    for v ∈ setdiff(ivertices(tree), mrca(tree))
        isone(length(dads(tree, v))) && length(children(tree, v)) == 2 && continue

        @info "Invalid internal vertex" vertex = v parents = dads(tree, v) children = children(tree, v)
        return false
    end

    ## Leaves must have no children and 1 parent
    for v ∈ leaves(tree)
        isone(length(dads(tree, v))) && isempty(children(tree, v)) && continue

        @info "Invalid leaf" vertex = v parents = dads(tree, v) children = children(tree, v)
        return false
    end

    true
end
