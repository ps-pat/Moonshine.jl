using Graphs

import Graphs: add_edge!

using Random

using Distributions

using StatsBase: sample, FrequencyWeights

import Base: iterate, length, size, isempty

using SpecialFunctions: loggamma

###################
# Tree Definition #
###################

export Tree
struct Tree <: AbstractGenealogy
    graph::SimpleDiGraph{VertexType}
    latitudes::Vector{Float64}
    sequences::Vector{Sequence}
    sample::Sample
    logprob::Base.RefValue{Float64x2}
end

function Tree(sample::Sample)
    n = length(sample)

    sequences = Vector{Sequence}(undef, 2n - 1)
    sequences[1:n] .= sample.H

    Tree(SimpleDiGraph(2n - 1),
         zeros(Float64, n - 1),
         sequences,
         sample, Ref(zero(Float64x2)))
end

###############################
# AbstractGenealogy Interface #
###############################

nleaves(tree::Tree) = (length(sequences(tree)) + 1) ÷ 2

describe(::Tree, long = true) = long ? "Coalescent Tree" : "Tree"

mrca(tree::Tree, vs, ::Any) = mrca(tree, vs)

###########################
# AbstractGraph Interface #
###########################

add_edge!(tree::Tree, e) = add_edge!(graph(tree), e)
add_edge!(tree::Tree, v1, v2) = add_edge!(graph(tree), v1, v2)

###########
# Methods #
###########

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

## TODO: implement `mutation_edges!` for Tree.

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
    μ = mut_rate(tree, true) * sam(tree).sequence_length
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
        tree.logprob[] += log(norms[v1_idx]) - log(norms.sum)
        v1 = live[v1_idx]
        live[v1_idx] = live[nlive]
        norms[v1_idx] = norms[nlive]
        norms[nlive] = 0
        nlive -= 1

        ## Sample second sequence
        potential2 =
            let η1 = sequence(tree, v1),
                η1_d0 = distance(Dist, η0, η1)

                function (η)
                    d = distance(Dist, η1, η)
                    if distance(Dist, η0, η) > η1_d0
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
        tree.logprob[] += logpdf(Gumbel(gumbel_μ), gumbel_x)

        ## Sample coalescence latitude
        Δcoal_dist = Exponential(inv(nlive))
        Δ = rand(rng, Δcoal_dist)
        tree.logprob[] += logpdf(Δcoal_dist, Δ)

        ## Add some latitude to account for mutations.
        d12 = distance(Dist, sequence(tree, v1), sequence(tree, v2))
        if !iszero(d12)
            Δmut = randexp(rng) / μ
            for _ ∈ 2:d12
                Δmut_new = randexp(rng) / μ
                Δmut_new ≤ Δmut && continue
                Δmut = Δmut_new
            end
            tree.logprob[] += log(μ) + log(d12) +
                              (d12 - 1) * log(1 - exp(-μ * Δmut)) - μ * Δmut
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

function validate(tree::Tree)
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
