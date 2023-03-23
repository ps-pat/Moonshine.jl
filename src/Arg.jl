using Graphs:
    SimpleDiGraph,
    SimpleEdge,
    AbstractSimpleGraph,
    degree,
    src, dst

import Graphs:
    edges,
    edgetype,
    has_edge,
    has_vertex,
    inneighbors,
    ne,
    nv,
    outneighbors,
    vertices,
    is_directed,
    add_vertex!,
    add_edge!,
    rem_edge!

import Base:
    eltype,
    show,
    union, intersect

using Random:
    AbstractRNG,
    GLOBAL_RNG,
    shuffle!

using GraphMakie:
    graphplot,
    register_interaction!,
    deregister_interaction!,
    NodeDrag,
    NodeHoverHighlight

using Distributions:
    Bernoulli,
    Exponential,
    Uniform,
    DiscreteUniform,
    Categorical,
    truncated,
    logccdf,
    logpdf

using IntervalSets

const VertexType = Int
const EdgeType = SimpleEdge{VertexType}
const Ω = Interval{:closed, :open, Float64}
const ∞ = Inf

######################
# ArgCore definition #
######################

struct ArgCore{T}
    graph::SimpleDiGraph{VertexType}
    latitudes::Vector{Float64}
    sequences::Vector{Sequence{T}}
    nleaves::Int
    ancestral_interval::Dict{EdgeType,
                             Set{Interval{:closed, :open, Float64}}}

    positions::Vector{Float64}
    seq_length::Float64
    eff_popsize::Float64
    μ_loc::Float64
    ρ_loc::Float64
end

## TODO: Check performance impact of sizehint!.
function ArgCore(leaves::AbstractArray{Sequence{T}};
                 positions = [],
                 seq_length = 1.0,
                 effective_popsize = 1.0,
                 μ_loc = 0.0,
                 ρ_loc = 0.0) where T
    n = length(leaves)
    latitudes = sizehint!(Float64[], 10n)
    sequences = sizehint!(leaves, 10n)

    nmarkers = (length ∘ first)(leaves)
    if length(positions) != nmarkers || !issorted(positions)
        @info((isempty(positions) ? "A" : "Invalid positions: a") *
            "ssuming equally spaced markers")
        positions = (collect ∘ range)(0, 1, length = nmarkers)
    end
    if minimum(positions) < 0
        @info("First position is < 0: shifting positions right")
        positions .-= minimum(positions)
    end
    if maximum(positions) > 1
        @info("Last position is > 0: scaling positions")
        positions ./= maximum(positions)
    end

    ArgCore{T}(SimpleDiGraph(n), latitudes, sequences, n, Dict(),
               positions, seq_length, effective_popsize, μ_loc, ρ_loc)
end

function ArgCore{T}(rng::AbstractRNG,
                    nmin::Integer, minlength::Integer,
                    nmax::Integer = 0, maxlength::Integer = 0;
                    genpars...) where T
    n = iszero(nmax) ? nmin : rand(rng, nmin:nmax)
    nmarkers = iszero(maxlength) ? minlength : rand(rng, minlength:maxlength)

    ArgCore([Sequence{T}(rng, nmarkers) for _ ∈ 1:n]; genpars...)
end

function ArgCore{T}(nmin::Integer, minlength::Integer,
                    nmax::Integer = 0, maxlength::Integer = 0;
                    genpars...) where T
    ArgCore{T}(GLOBAL_RNG, nmin, minlength, nmax, maxlength; genpars...)
end

function union(x::T, xs::Set{T}) where T
    Set([x]) ∪ xs
end

intersect(x::T, xs::Set{T}) where T =
    (Set ∘ broadcast)(Fix1(intersect, x), xs)

for fun ∈ [:union, :intersect]
    @eval $fun(xs::Set{T}, x::T) where T = $fun(x, xs)
end

## TODO: simplify sets of intervals?

#############################
# End of ArgCore definition #
#############################

##################
# Arg definition #
##################

mutable struct Arg{T} <: AbstractSimpleGraph{VertexType}
    core::ArgCore{T}
    nrecombinations::Int
    logprob::BigFloat
end
export Arg

Arg(leaves::AbstractArray{Sequence{T}}, genpars...) where T =
    Arg(ArgCore{T}(leaves; genpars...), 0, zero(BigFloat))

function Arg{T}(rng::AbstractRNG,
                nmin::Integer, minlength::Integer,
                nmax::Integer = 0, maxlength::Integer = 0;
                genpars...) where T
    Arg(ArgCore{T}(rng, nmin, minlength, nmax, maxlength; genpars...),
        0, zero(BigFloat))
end
function Arg(rng::AbstractRNG,
             nmin::Integer, minlength::Integer,
             nmax::Integer = 0, maxlength::Integer = 0;
             genpars...)
    Arg{UInt}(rng, nmin, minlength, nmax, maxlength; genpars...)
end

function Arg{T}(nmin::Integer, minlength::Integer,
                nmax::Integer = 0, maxlength::Integer = 0;
                genpars...) where T
    Arg{T}(GLOBAL_RNG, nmin, minlength, nmax, maxlength; genpars...)
end

function Arg(nmin::Integer, minlength::Integer,
             nmax::Integer = 0, maxlength::Integer = 0;
             genpars...)
    Arg{UInt}(nmin, minlength, nmax, maxlength; genpars...)
end

## AbstractGraph Interface.
## See https://juliagraphs.org/Graphs.jl/stable/ecosystem/interface/
for fun ∈ [:edges, :vertices, :ne, :nv, :add_vertex!]
    @eval function $fun(arg::Arg)
        $fun(arg.core.graph)
    end
end

for (fun, ret) ∈ Dict(:eltype => VertexType,
                      :edgetype => EdgeType,
                      :is_directed => :true)
    @eval begin
        @generated $fun(::Type{Arg}) = $ret
        @generated $fun(::Arg) = $fun(Arg)
    end
end

let arg_edge = Expr(:(::), :e, EdgeType),
    arg_vertex = Expr(:(::), :v, VertexType)
    for (fun, a) ∈ Dict(:has_edge => arg_edge,
                        :add_edge! => arg_edge,
                        :rem_edge! => arg_edge,
                        :has_vertex => arg_vertex,
                        :inneighbors => arg_vertex,
                        :outneighbors => arg_vertex)
        varname = first(a.args)
        @eval function $fun(arg::Arg, $a)
            $fun(arg.core.graph, $varname)
        end
    end
end

## Simple methods.
for field ∈ [:(:nleaves), :(:sequences), :(:latitudes),
             :(:seq_length), :(:eff_popsize), :(:positions)]
    fun_name = eval(field)
    @eval begin
        export $fun_name

        function $fun_name(arg)
            getfield(arg.core, $field)
        end
    end
end

nrecombinations(arg) = arg.nrecombinations

leaves(arg) = Base.OneTo(arg.core.nleaves)
isleaf(arg, v) = v <= nleaves(arg)

function ivertices(arg)
    n = nleaves(arg)
    range(n + 1, length = n - 1)
end
nivertices(arg) = nleaves(arg) - 1

nmarkers(arg) = (length ∘ first)(arg.core.sequences)

sequences(arg, v::VertexType) = sequences(arg)[v]
sequences(arg, e::EdgeType) = [sequences(arg, src(e)), sequences(arg, dst(e))]

export latitude
latitude(arg, v) =
    isleaf(arg, v) ? zero(Float64) : arg.core.latitudes[v - nleaves(arg)]

export mrca, tmrca
mrca(arg) = isempty(arg.core.latitudes) ?
    zero(Int) : argmax(arg.core.latitudes) + nleaves(arg)
tmrca(arg) = isempty(arg.core.latitudes) ?
    zero(Float64) : maximum(arg.core.latitudes)

"""
    children(arg, v, int = 0 .. ∞)
    parents(arg, v, int = 0 .. ∞)

Return the children/parents of a vertex. Optionally, only return those bearing
ancestral material for a given interval.
"""
function children end,
function parents end
export children, parents

children(arg, v) = outneighbors(arg, v)

parents(arg, v) = inneighbors(arg, v)

for fun ∈ [:children, :parents]
    @eval function $fun(arg, v, int)
        filter($fun(arg, v)) do x
            ancestral_intervals(arg, x) ∩ int |> !isempty
        end
    end
end

"""
    idxtopos(arg, idx)

Return the position of the marker given its index.
"""
idxtopos(arg, idx) = arg.core.positions[idx]

"""
    postoidx(arg, pos)

Return the largest marker's index that is at a position lesser than the one
given.
"""
function postoidx(arg, pos)
    @inbounds for (k, p) ∈ enumerate(arg.core.positions)
        p > pos && return k - 1
    end

    nmarkers(arg)
end

"""
    ancestral_intervals(arg, e)
    ancestral_intervals(arg, σ, δ)
    ancestral_intervals(arg, v)

Compute the interval for which an edge or a vertex is ancestral.
"""
function ancestral_intervals end
export ancestral_intervals

ancestral_intervals(arg, e::EdgeType) = arg.core.ancestral_interval[e]
ancestral_intervals(arg, σ, δ) = ancestral_intervals(arg, EdgeType(σ, δ))

function ancestral_intervals(arg, v::VertexType)
    isleaf(arg, v) && return IntervalSet(Ω(0, ∞))
    reduce(outneighbors(arg, v)) do x, y
        ancestral_intervals(arg, v, x) ∪ ancestral_intervals(arg, v, y)
    end
end

@generated blocksize(arg::Arg{T}) where T = blocksize(Sequence{T})

export nmutations
nmutations(arg, e) = (count_ones ∘ xor)(sequences(arg, e)...)
nmutations(arg) = mapreduce(Fix1(nmutations, arg), +, edges(arg))

function branchlength_tree(arg)
    lats = view(latitudes(arg), range(1, length = nleaves(arg) - 1))
    sum(lats, init = last(lats))
end

export branchlength
function branchlength(arg)
    ν = nrecombinations(arg)
    lats = view(latitudes(arg), range(nleaves(arg), length = 2ν))
    rec_branchlength = mapreduce(vs -> last(vs) - first(vs), +,
                                 Iterators.partition(lats, 2), init = 0.0)

    branchlength_tree(arg) + rec_branchlength
end

for (fun, par) ∈ Dict(:mut_rate => :(:μ_loc), :rec_rate => :(:ρ_loc))
    @eval begin
        export $fun

        function $fun(arg, scaled = true)
            rate_loc = getfield(arg.core, $par)
            mul = scaled ? 4eff_popsize(arg) : one(rate_loc)

            mul * rate_loc
        end
    end
end

#########################
# End of Arg definition #
#########################

##############################
# Plotting & pretty printing #
##############################

function argplot end
export argplot

function argplot(arg, x;
                 arrow_show = true,
                 wild_color = :blue,
                 derived_color = :red,
                 attributes...)
    n = nv(arg)
    nedges = ne(arg)

    edgecolor = fill(:gray, nedges)

    nodecolor = map(sequences(arg)) do vertex
        vertex[x] ? derived_color : wild_color
    end

    node_labels = string.(1:n)

    node_sizes = fill(25.0, n)

    p = graphplot(arg,
                  nlabels = node_labels,
                  nlabels_distance = 10,
                  node_size = node_sizes,
                  node_color = nodecolor,
                  edge_color = edgecolor,
                  edge_width = fill(3, nedges),
                  elabels = nothing,
                  arrow_show = arrow_show,
                  attributes...)


    deregister_interaction!(p.axis, :rectanglezoom)
    register_interaction!(p.axis, :nodedrag, NodeDrag(p.plot))
    register_interaction!(p.axis, :nodehover, NodeHoverHighlight(p.plot))

    p
end

function show(io::IO, ::MIME"text/plain", arg::Arg)
    println(io, "Ancestral recombination graph:")
    println(io, nleaves(arg), " leaves, ",
            nmarkers(arg), " markers")
    print(io, "tMRCA: ", tmrca(arg))
end

function show(io::IO, arg::Arg)
    print(io, "ARG(")
    print(io, nleaves(arg))
    print(io, " leaves, ")
    print(io, nmarkers(arg))
    print(io, " markers)")
end

#####################################
# End of plotting & pretty printing #
#####################################

####################
# ARG construction #
####################

function coalesce!(rng, arg, vertices, nlive)
    ## Select coalescing pair.
    arg.logprob += (log ∘ inv ∘ binomial)(length(vertices), 2)
    shuffle!(rng, vertices)

    ## Sample a latitude for the coalescence event.
    Δdist = Exponential(inv(nlive))
    Δ = rand(rng, Δdist)
    arg.logprob += logccdf(Δdist, Δ)

    ## Perform the coalescence.
    add_vertex!(arg) || @error "Could not add a vertex to ARG"
    parent = nv(arg)
    _children = Vector{eltype(arg)}(undef, 2)
    for child ∈ 1:2
        _children[child] = pop!(vertices)
        e = EdgeType(parent, _children[child])
        add_edge!(arg, e)

        arg.core.ancestral_interval[e] = Set([Ω(0, ∞)])
    end

    push!(vertices, parent)

    push!(sequences(arg),
          sequences(arg)[first(_children)] &
              sequences(arg)[last(_children)])

    push!(latitudes(arg), latitude(arg, parent - 1) + Δ)

    @info("$(first(_children)) and $(last(_children)) \
           coalesced into $parent at $(latitude(arg, parent))")
end

"""
    buildtree!(rng = GLOBAL_RNG, arg, idx = 1)

Build a marginal tree consistent with a specific index.
"""
function buildtree! end
export buildtree!

function buildtree!(rng::AbstractRNG, arg::Arg, idx = 1)
    nv(arg) > nleaves(arg) && @warn "ARG contains non-leaf vertices"

    ## Group wild/derived sequences together to make the algorithm simpler.
    sort!(sequences(arg), lt = (x, y) -> x[idx] <= y[idx])

    nlive = nleaves(arg)
    nlive_wild = map(s -> s[idx], sequences(arg)) .|> iszero |> sum
    nlive_derived = nlive - nlive_wild

    wild = (collect ∘ range)(1, length = nlive_wild)
    derived = (collect ∘ range)(nlive_wild + 1, length = nlive_derived)

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
        arg.logprob += logpdf(coalescence_dist, derived_coalescence)

        ## Perform the coalescence.
        coalesce!(rng, arg, cvertices, nlive)

        ## Update the number of live nodes.
        nlive -= 1
        nlive_derived -= derived_coalescence
        nlive_wild -= !derived_coalescence
    end

    cvertices = [derived; wild]
    while nlive > 1
        coalesce!(rng, arg, cvertices, nlive)
        nlive -= 1
    end

    arg
end
buildtree!(arg::Arg, idx = 1) = buildtree!(GLOBAL_RNG, arg, idx)

###########################
# End of ARG construction #
###########################

#######
# MMN #
#######

"""
    mmn_blocks_idx(arg, e, int = Ω(0, ∞))

Compute the indices of blocks to be used in a pass of the mmn algorithm.
"""
function mmn_blocks_idx end

function mmn_blocks_idx(arg::Arg{T}, e::EdgeType, int = Ω(0, ∞)) where T
    target_intervals = ancestral_intervals(arg, e) ∩ int

    n = nmarkers(arg)
    bs = blocksize(arg)

    (sort ∘ mapreduce)(union, target_intervals, init = Int[]) do i
        _endpoints = endpoints(i)

        postoidx_arg = Fix1(postoidx, arg)
        lbound = (postoidx_arg ∘ first)(_endpoints)
        rbound = (postoidx_arg ∘ prevfloat ∘ last)(_endpoints)

        blockidx(n, bs, lbound), blockidx(n, bs, rbound)
    end
end

mmn_blocks_idx(arg, σ::VertexType, δ::VertexType, int = Ω(0, ∞)) =
    mmn_blocks_idx(arg, EdgeType(σ, δ), int)

function first_inconsistent_position end
export first_inconsistent_position

function first_inconsistent_position(arg, int::Ω)
    n = nmarkers(arg)
    _mrca = mrca(arg)
    (iszero(_mrca) || iszero(n)) && return ∞, EdgeType[]

    ## Vertices to xor are stored in acc.
    acc = zeros(VertexType, nv(arg))
    acc_ptr = firstindex(acc)
    acc[acc_ptr] = _mrca
    acc_ptr += 1

    ## Contains mutation edges.
    mutation_edges = Dict{Int, Vector{EdgeType}}()

    bs = blocksize(arg)
    for k ∈ eachindex(acc)
        vertex = acc[k]
        iszero(vertex) && break

        for e ∈ map(Fix1(EdgeType, vertex), children(arg, vertex))
            blocks_idx = mmn_blocks_idx(arg, e, int)

            ## Non ancestral edge.
            isempty(blocks_idx) && continue

            ## Add child to acc.
            c = dst(e)
            if !isleaf(arg, c)
                (acc[acc_ptr] = dst(e))
                acc_ptr += 1
            end

            xored_sequences =
                sequences(arg, src(e)).data[blocks_idx] .⊻
                sequences(arg, dst(e)).data[blocks_idx]
            for (block_idx, block) ∈ zip(blocks_idx, xored_sequences)
                first_set_bit = leading_zeros(block) + 1
                if first_set_bit <= bs
                    idx = actualpos(n, bs, block_idx, first_set_bit)
                    haskey(mutation_edges, idx) || (mutation_edges[idx] = [])
                    push!(mutation_edges[idx], e)
                    break
                end
            end
        end
    end

    inconsistent_idx =
        minimum(p -> (length ∘ last)(p) > 1 ? first(p) : typemax(Int),
                mutation_edges)

    idxtopos(arg, inconsistent_idx), mutation_edges[inconsistent_idx]
end

first_inconsistent_position(arg, lbound::Real, ubound::Real = ∞) =
    first_inconsistent_position(arg, Ω(lbound, ubound))

##############
# End of MMN #
##############
