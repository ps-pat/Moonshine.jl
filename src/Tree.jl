using Graphs

import Graphs: add_edge!

using Random

using Distributions

using StatsBase: sample, FrequencyWeights

import Base: iterate, length, size, isempty

using SpecialFunctions: loggamma

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
    μloc::Float64
end

TreeCore() = TreeCore(SimpleDiGraph{VertexType}(), [], [], [],
                      typemin(Float64), typemin(Float64), typemin(Float64))

function TreeCore(leaves::AbstractVector{Sequence};
                  positions = nothing,
                  seq_length = one(Float64),
                  Ne = one(Float64),
                  μloc = 1e-7)
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

    TreeCore(SimpleDiGraph(2n - 1), zeros(Float64, n - 1), sequences,
             positions, seq_length, Ne, μloc)
end

###################
# Tree Definition #
###################

export Tree
mutable struct Tree <: AbstractGenealogy
    core::TreeCore
    logprob::Float64x2
    nextvertex::Int
end

Tree() = Tree(TreeCore(), zero(BigFloat), typemin(Int))

Tree(core::TreeCore, logprob) = Tree(core, logprob, one(Int))

function Tree(leaves::AbstractVector{Sequence}; genpars...)
    Tree(TreeCore(leaves; genpars...), zero(BigFloat), one(Int))
end

###############################
# AbstractGenealogy Interface #
###############################

nleaves(tree::Tree) = (length(sequences(tree)) + 1) ÷ 2

describe(::Tree, long = true) = long ? "Coalescent Tree" : "Tree"

function mrca(tree::Tree)
    any(iszero, latitudes(tree)) && return zero(VertexType)
    isone(nv(tree)) && return one(VertexType)
    argmax(latitudes(tree)) + nleaves(tree)
end

###########################
# AbstractGraph Interface #
###########################

add_edge!(tree::Tree, e) = add_edge!(graph(tree), e)
add_edge!(tree::Tree, v1, v2) = add_edge!(graph(tree), v1, v2)

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

export build!
"""
    build!(rng, tree, distance; bias0 = ∞, toilet_prop = 1)

Build a coalescent tree

`tree` must contain 2n-1 vertices (where n is the number of leaves) and no edge.
"""
function build! end

function _sample_toilet(rng, xs, potential, threshold_prop)
    threshold = floor(Int, length(xs) * threshold_prop)
    iter = (Iterators.Stateful ∘ enumerate)(xs)
    z = ∞
    μ = zero(Float64)
    ret = zero(Int)
    for (k, x) ∈ iter
        z = potential(x)
        isinf(z) && continue
        μ += z
        z -= log(randexp(rng))
        ret = k
        break
    end

    for (k, x) ∈ iter
        newz = potential(x)
        isinf(newz) && continue
        tmin, tmax = minmax(newz, μ)
        μ = tmax + (log1p ∘ exp)(tmin - tmax)

        newz -= log(randexp(rng))
        newz <= z && continue

        z = newz
        ret = k

        k <= threshold || break
    end

    ret, (z, μ)
end

function build!(rng, tree::Tree;
                Dist::Distance = Hamming{Int}(), bias0 = ∞, toilet_prop = 1)
    n = nleaves(tree)
    nv(tree) ≠ 2n - 1 || !iszero(ne(tree)) && error("Invalid tree")
    μ = mut_rate(tree, true)
    η0 = zeros(Sequence, nmarkers(tree))

    ## Compute the norm of the sequences. They are stored as
    ## FrequencyWeights since they will be used for sampling sequences.
    norms = FrequencyWeights([distance(Dist, η0, η) + 1 for η ∈ sequences(tree, 1:n)])

    live = collect(leaves(tree))
    for nlive ∈ range(length(live), 2, step = -1)
        ## Sample first sequence
        v1_idxs = sample(rng, 1:nlive, norms[1:nlive], 2, replace = false)
        v1_idx = norms[first(v1_idxs)] > norms[last(v1_idxs)] ?
            first(v1_idxs) : last(v1_idxs)
        tree.logprob += log(norms[v1_idx]) - log(norms.sum)
        v1 = live[v1_idx]
        live[v1_idx] = live[nlive]
        norms[v1_idx] = norms[nlive]
        norms[nlive] = 0
        nlive -= 1

        ## Sample second sequence
        potential2 =
            let η1 = sequence(tree, v1),
                η1_d0 = distance(Dist, η0, η1)
                function(η)
                    d = distance(Dist, η1, η)
                    if distance(Dist, η0, η) > η1_d0
                        ## ∞ - loggamma(∞) incorrectly returns NaN
                        # isinf(bias0) && return -∞
                        d += bias0
                    end

                    iszero(d) && return zero(Float64)

                    ## Use Stirling approximation.
                    logd = log(d)
                    ret = d * (log(μ) - logd + 1) - 0.5 * (log2π + logd)
                    ret
                end
            end

        v2_idx, (gumbel_x, gumbel_μ) =
            _sample_toilet(rng,
                           sequences(tree, live[1:nlive]),
                           potential2,
                           toilet_prop)
        v2 = live[v2_idx]
        v = 2n - nlive
        live[v2_idx] = v
        if !isinf(gumbel_x)
            tree.logprob += logpdf(Gumbel(gumbel_μ), gumbel_x)
        end

        ## Sample coalescence lattitude
        Δcoal_dist = Exponential(inv(nlive))
        Δ = rand(rng, Δcoal_dist)
        tree.logprob += logpdf(Δcoal_dist, Δ)

        d12 = distance(Dist, sequence(tree, v1), sequence(tree, v2))
        if !iszero(d12)
            Δmut_dist = Gamma(d12, inv(μ))
            Δmut = rand(rng, Δmut_dist)
            tree.logprob += logpdf(Δmut_dist, Δmut)
            Δ += Δmut
        end

        ## Add coalescence event to tree
        add_edge!(tree, v, v1)
        add_edge!(tree, v, v2)

        sequences(tree)[v] = sequence(tree, v1) & sequence(tree, v2)
        norms[v2_idx] = distance(Dist, η0, sequence(tree, v)) + 1

        latitudes(tree)[n - nlive] = latitude(tree, v - 1) + Δ
    end

    tree
end

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
