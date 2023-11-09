using Graphs: SimpleDiGraph

import Graphs: add_vertex!, add_edge!, rem_edge!

using Random: AbstractRNG, GLOBAL_RNG

#######################
# TreeCore Definition #
#######################

struct TreeCore{T}
    graph::SimpleDiGraph{VertexType}
    latitudes::Vector{Float64}
    sequences::Vector{Sequence{T}}
    nleaves::Int

    positions::Vector{Float64}
    seq_length::Float64
    eff_popsize::Float64
    μ_loc::Float64
end

function TreeCore{T}(leaves::AbstractVector{Sequence{T}};
                     positions = (collect ∘ range)(0, 1, length = length(leaves)),
                     seq_length = one(Float64),
                     effective_popsize = one(Float64),
                     μ_loc = zero(Float64)) where T
    n = length(leaves)

    sequences = similar(leaves, 2n - 1)
    for (k, leaf) ∈ enumerate(leaves)
        sequences[k] = leaf
    end

    positions = validate_positions(positions, (length ∘ first)(leaves))

    TreeCore{T}(SimpleDiGraph(n), zeros(Float64, n - 1), sequences, n,
                positions, seq_length, effective_popsize, μ_loc)
end

function TreeCore{T}(rng::AbstractRNG,
                     nmin::Integer, minlength::Integer,
                     nmax::Integer = 0, maxlength::Integer = 0;
                     genpars...) where T
    n = iszero(nmax) ? nmin : rang(rng, nmin:nmax)
    nmarkers = iszero(maxlength) ? minlength : rang(rng, minlength:maxlength)

    TreeCore{T}([Sequence{T}(rng, nmarkers) for _ ∈ 1:n]; genpars...)
end

function TreeCore{T}(nmin::Integer, minlength::Integer,
                     nmax::Integer = 0, maxlength::Integer = 0;
                     genpars...) where T
    TreeCore{T}(GLOBAL_RNG, nmin, minlength, nmax, maxlength; genpars...)
end


###################
# Tree Definition #
###################

export Tree
mutable struct Tree{T} <: AbstractGenealogy
    core::TreeCore{T}
    logprob::BigFloat
    nextvertex::Int
end

Tree(leaves::AbstractVector{Sequence{T}}, genpars...) where T =
    Tree(TreeCore{T}(leaves; genpars...), zero(BigFloat), one(Int))

function Tree{T}(rng::AbstractRNG,
                 nmin::Integer, minlength::Integer,
                 nmax::Integer = 0, maxlength::Integer = 0;
                 genpars...) where T
    Tree(TreeCore{T}(rng, nmin, minlength, nmax, maxlength; genpars...),
         zero(BigFloat), one(Int))
end

function Tree(rng::AbstractRNG,
              nmin::Integer, minlength::Integer,
              nmax::Integer = 0, maxlength::Integer = 0;
              genpars...)
    Tree{UInt}(rng, nmin, minlength, nmax, maxlength; genpars...)
end

function Tree{T}(nmin::Integer, minlength::Integer,
                 nmax::Integer = 0, maxlength::Integer = 0;
                 genpars...) where T
    Tree{T}(GLOBAL_RNG, nmin, minlength, nmax, maxlength; genpars...)
end

function Tree(nmin::Integer, minlength::Integer,
              nmax::Integer = 0, maxlength::Integer = 0;
              genpars...)
    Tree{UInt}(nmin, minlength, nmax, maxlength; genpars...)
end

##################
# Simple Methods #
##################

for field ∈ [:(:nleaves), :(:sequences), :(:latitudes), :(:graph),
             :(:seq_length), :(:eff_popsize), :(:positions)]
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

latitude(tree::Tree, v) =
    isleaf(tree, v) ? zero(Float64) : latitudes(tree)[v - nleaves(tree)]

leaves(tree::Tree) = Base.OneTo(nleaves(tree))

ivertices(tree::Tree) = range(nleaves(tree) + 1, nv(tree))

mrca(tree::Tree) = any(iszero.(latitudes(tree))) ?
    zero(VertexType) : argmax(latitudes(tree)) + nleaves(tree)

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
mut_rate(tree::Tree, scaled = true) =
    tree.core.μ_loc * (scaled ? 4 * eff_popsize(tree) : 1)

function tree_coalesce!(rng, tree, vertices, nlive)
    ## Sample coalescing vertices.
    shuffle!(rng, vertices)
    child1, child2 = pop!(vertices), pop!(vertices)
    tree.logprob += (log ∘ inv ∘ binomial)(length(vertices), 2)

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

    tree
end

build!(tree::Tree, idx = 1) = build!(GLOBAL_RNG, tree, idx)
