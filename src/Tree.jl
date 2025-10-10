using Graphs

using Random

using StatsBase: sample, FrequencyWeights

import Base: iterate, length, size, isempty

using SpecialFunctions: loggamma

###################
# Tree Definition #
###################

export Tree
"""
    $(TYPEDEF)

Coalescent tree.

See also [`Arg`](@ref).

# Fields
$(TYPEDFIELDS)

# Constructors
!!! info
    Random constructor calls [`Sample`](@ref)'s random constructor.

!!! warning
    These **do not** actually build the tree. For that, see
    [`build!(rng, tree)`](@ref).

$(METHODLIST)

# Arguments
Arguments are the same as for [`Sample`](@ref).
"""
struct Tree <: AbstractGenealogy
    "Tree's topology"
    graph::ThreeTree{VertexType}
    "Vertices' latitudes"
    latitudes::Vector{Float64}
    "Vertices' haplotypes"
    sequences::Vector{Sequence}
    "Associated [`Sample`](@ref)"
    sample::Sample
    "Log-value of the associated pdf"
    logdensity::Base.RefValue{Double64}
end

function Tree(sample::Sample)
    n = length(sample)

    sequences = Vector{Sequence}(undef, 2n - 1)
    sequences[1:n] .= sample.H

    Tree(ThreeTree(VertexType(n)),
         zeros(Float64, n - 1),
         sequences,
         sample, Ref(zero(Double64)))
end

function Tree(rng::AbstractRNG, n, μ, ρ, Ne, sequence_length)
    sample = Sample(rng, n, μ, ρ, Ne, sequence_length)
    Tree(sample)
end

###############################
# AbstractGenealogy Interface #
###############################

nleaves(tree::Tree) = (length(sequences(tree)) + 1) ÷ 2

describe(::Tree, long = true) = long ? "Coalescent Tree" : "Tree"

mrca(tree::Tree) = nv(tree)

mrca(tree::Tree, ::Any) = mrca(tree)

@generated maxdads(::Type{Tree}) = 1

@generated maxchildren(::Type{Tree}) = 2

iscoalescence(tree::Tree, v) = nleaves(tree) < v <= nv(tree)

###########
# Methods #
###########

export depth
"""
    $(SIGNATURES)

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

function descendants_leaves!(vertices, tree::Tree, v; buffer = default_buffer())
    n = nleaves(tree)

    vertices_len = 0

    @no_escape buffer begin
        store = @alloc(VertexType, n)
        stack = CheapStack(store)
        push!(stack, v)

        while !isempty(stack)
            x = pop!(stack)
            if isleaf(tree, x)
                vertices_len += 1
                vertices[vertices_len] = x
                continue
            end

            for child ∈ children(tree, x)
                push!(stack, child)
            end
        end
    end

    resize!(vertices, vertices_len)
end

descendants_leaves(tree::Tree, v; buffer = default_buffer()) =
    descendants_leaves!(Vector{VertexType}(undef, nleaves(tree)), tree, v,
                                           buffer = buffer)

#################
# Tree Building #
#################

function _sample_toilet(rng, xs, potential, threshold_prop)
    threshold = floor(Int, length(xs) * threshold_prop)
    z = ∞
    logπ = zero(z)
    logΣπ = zero(logπ)
    idx = zero(Int)

    k = firstindex(xs)
    while k <= lastindex(xs)
        x = xs[k]
        logπ = potential(x)

        if isfinite(logπ)
            logΣπ = logπ
            z = logπ - log(randexp(rng))
            idx = k

            k += 1
            break
        end

        k += 1
    end

    while k <= lastindex(xs)
        k <= threshold || break

        x = xs[k]
        logπ_new = potential(x)

        if isfinite(logπ_new)
            πmin, πmax = minmax(logπ_new, logΣπ)
            logΣπ = πmax + (log1p ∘ exp)(πmin - πmax)

            newz = logπ_new - log(randexp(rng))
            if newz >= z
                logπ = logπ_new
                z = newz
                idx = k
            end

        end

        k += 1
    end

    if iszero(idx) # All potentials were infinite
        idx = sample(rng, eachindex(xs))
        logπ = 0
        logΣπ = (log ∘ length)(xs)
    end

    idx, logπ - logΣπ
end

function build!(rng, tree::Tree;
                Dist::Distance = Hamming{Int}(), bias0 = 1, threshold_prop = 1)
    n = nleaves(tree)
    μ = mut_rate(tree, false)

    live = collect(leaves(tree))
    for nlive ∈ range(length(live), 2, step = -1)
        ## Sample first sequence
        v1_idx = sample(rng, 1:nlive)
        add_logdensity!(tree, -log(nlive))
        v1 = live[v1_idx]
        live[v1_idx] = live[nlive]
        nlive -= 1

        ## Sample second sequence
        potential2 =
            let η1 = sequence(tree, v1)
                function (η)
                    d = distance(Dist, η1, η)
                    iszero(d) && return zero(Float64)
                    d += d * bias0

                    ## Poisson potential with Stirling's approximation
                    d * (log(μ) - log(d) + 1)
                end
            end

        v2_idx, logprob = _sample_toilet(rng,
                                         sequences(tree, view(live, 1:nlive)),
                                         potential2,
                                         threshold_prop)
        v2 = live[v2_idx]
        v = 2n - nlive
        live[v2_idx] = v
        add_logdensity!(tree, logprob)

        ## Sample coalescence latitude
        Δcoal_dist = Exponential(inv(nlive))
        Δ = rand(rng, Δcoal_dist)
        add_logdensity!(tree, Δcoal_dist, Δ)

        ## Add coalescence event to tree
        add_coalescence_vertex!(tree.graph, v1, v2)

        sequences(tree)[v] = sequence(tree, v1) & sequence(tree, v2)

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
