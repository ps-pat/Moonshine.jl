using Graphs: SimpleEdge,
              AbstractSimpleGraph,
              diameter,
              indegree, outdegree

import Graphs: edges, vertices, ne, nv,
               eltype, edgetype, is_directed,
               has_edge, has_vertex, inneighbors, outneighbors

import GraphMakie: graphplot

using LayeredLayouts

using IntervalSets

using GeometryBasics: Point

using Bumper

using DataStructures: Stack, DefaultDict, RBTree

import Base: IteratorSize, eltype

const VertexType = Int
const EdgeType = SimpleEdge{VertexType}
const ∞ = Inf

export Ω
const Ω = Interval{:closed, :open, Float64}

abstract type AbstractGenealogy <: AbstractSimpleGraph{VertexType} end

#############
# Interface #
#############

export graph
"""
    graph(genealogy)

Return the underlying graph of a genealogy.
"""
function graph end

"""
    describe(genealogy, long = true)

Return a string containing a long or short description of a genealogy.
Used internally by pretty printing functions.
"""
function describe end

export latitudes
"""
    latitudes(genealogy[, ivs])

Return the latitudes of (a subset of) the internal vertices of a genealogy.

See also [`latitude`](@ref) to get the latitude of a single vertex.

# Implementation
A default implementation for `latitudes(::AbstractGenealogy, ::Any)` is
available.
"""
function latitudes end

latitudes(genealogy, ivs) = getindex(latitudes(genealogy), ivs)

export latitude
"""
    latitude(genealogy, v)

Latitude of a vertex.

See also [`latitudes`](@ref) to get the latitudes of internal vertices stored
in a genealogy.
"""
function latitude end

export leaves, ivertices
"""
    leaves(genealogy)
    ivertices(genealogy_vertex)

An interable containing the leaves/internal vertices of a genealogy.

See also [`nleaves`](@ref) and [`nivertices`](@ref) for the number of leaves
and internal vertices.

# Implementation

Default implementations assume that the first `nleaves(genealogy)` vertices
are the leaves of the genealogy.
"""
function leaves end,
function ivertices end

leaves(genealogy) = Base.OneTo(nleaves(genealogy))
ivertices(genealogy) = UnitRange(nleaves(genealogy) + 1, nv(genealogy))

export sequences
"""
    sequences(genealogy[, e])
    sequences(genealogy[, vs])

Sequences of a genealogy. If an edge is specified, return the sequences
associated with the vertices incident to that edge. If an iterable of
vertices is specified, return the sequences associated with these vertices.

See also [`sequence`](@ref) to get the sequence associated with a vertex.

# Implementation

Custom types only need to implement `sequences(::T)`.
"""
function sequences end

sequences(genealogy, vs) = Iterators.map(v -> sequence(genealogy, v), vs)

sequences(genealogy, e::EdgeType) = (sequence(genealogy, src(e)),
                                     sequence(genealogy, dst(e)))

export mrca
"""
    mrca(genealogy[, vs])

Most recent common ancestor of a set of vertices. If omited, returns the mrca
of the whole genealogy.

See also [`tmrca`](@ref) for the time to the most recent common ancestor.
"""
function mrca end

function mrca(genealogy, pos::Real, vs = leaves(genealogy))
    length(vs) < 2 && return zero(VertexType)

    μ = argmax(latitudes(genealogy)) + nleaves(genealogy)
    vstack = Stack{VertexType}(ceil(Int, log(nv(genealogy))))
    push!(vstack, children(genealogy, μ, pos)...)

    ds = Set{VertexType}()
    while !isempty(vstack)
        μ2 = pop!(vstack)
        isleaf(genealogy, μ2) && continue
        vs ⊆ descendants!(ds, genealogy, μ2, pos) || continue

        μ = μ2
        push!(vstack, children(genealogy, μ, pos)...)
    end

    μ
end

function mrca(genealogy, vs::AbstractVector = leaves(genealogy))
    length(vs) < 2 && return zero(VertexType)

    μ = argmax(latitudes(genealogy)) + nleaves(genealogy)
    predicate = child ->
        isempty(ancestral_intervals(genealogy, Edge(μ, child)) ∩ Ω(0, 1))
    @inbounds while any(predicate, children(genealogy, μ))
        k = findfirst(!predicate, children(genealogy, μ))
        μ = children(genealogy, μ)[k]
    end

    μ
end

export prob
"""
    prob(genealogy; logscale = false)

Probability of the genealogy.

# Implementation,

Only mandatory if copulas are to be fitted on the genealogy.
"""
function prob end

export positions
"""
    positions(genealogy)

Positions of the markers.
"""
function positions end

export ancestral_intervals!, ancestral_intervals
"""
    ancestral_intervals!(ωs, genealogy, x)
    ancestral_intervals(genealogy, x)

Interval for which x is ancestral. Default implementation assumes that anything
is ancestral for [0, ∞).
"""
function ancestral_intervals end

function ancestral_intervals!(ωs, ::Any, ::Any)
    empty!(ωs)
    push!(ωs, Ω(0, ∞))
end

@generated ancestral_intervals(genealogy::Any, x::Any) =
    ancestral_intervals!(Set{Ω}(), genealogy, x)

## Recombinations ##

export nrecombinations
"""
    nrecombinations(genealogy)

Number of recombinations in a genealogy.

Default implementation returns 0.
"""
function nrecombinations end

@generated nrecombinations(::Any) = zero(Int)

export recombinations
"""
    recombinations(genealogy)

Iterator over the recombination vertices of a genealogy.

Default implementation returns an empty iterator.
"""
function recombinations end

@generated recombinations(::Any) = StepRange{Int, Int}(0, 1, 0)

export isrecombination
"""
    isrecombination(genealogy, v)

Returns true if `v` is a recombination for `genealogy`.

Default implementation always returns `false`.
"""
function isrecombination end

@generated isrecombination(::Any, ::Any) = false

#############
# Utilities #
#############

"""
    validate_positions(positions, nmarkers)

Validate a marker of positions. Tries to fix it if invalid. Return a vector
of uniformly spaced positions if it cannot be fixed.
"""
function validate_positions(positions, nmarkers)
    if length(positions) != nmarkers || !issorted(positions)
        @info((isempty(positions) ? "A" : "Invalid positions: a") *
              "ssuming equally spaced markers")
        positions = isone(nmarkers) ?
                    [0.0] : (collect ∘ range)(0, 1, length = nmarkers)
    end
    if minimum(positions) < 0
        @info("First position is < 0: shifting positions right")
        positions .-= minimum(positions)
    end
    if maximum(positions) > 1
        @info("Last position is > 0: scaling positions")
        positions ./= maximum(positions)
    end

    positions
end

"""
    idxtopos(genealogy, idx)

Return the position of the marker given its index.
"""
idxtopos(genealogy, idx) = getindex(positions(genealogy), idx)

"""
    postoidx(genealogy, pos)

Return the largest marker's index that is at a position lesser than the one
given.
"""
function postoidx(genealogy, pos)
    @inbounds for (k, p) ∈ enumerate(positions(genealogy))
        p > pos && return k - 1
    end

    nmarkers(genealogy)
end

"""
    ancestral_mask!(η, genealogy, x; wipe = true)
    ancestral_mask(genealogy, x)

Mask non ancestral positions to 0. If `wipe = true`, all markers in `η` wil be
initialized at 0.
"""
function ancestral_mask! end,
function ancestral_mask end

function ancestral_mask!(η, genealogy, ω::Ω; wipe = true)
    lpos, rpos = endpoints(ω)

    lidx = 1
    while lpos > positions(genealogy)[lidx]
        lidx += 1
    end

    ridx = postoidx(genealogy, rpos)

    wipe && wipe!(η)
    η.data[range(lidx, ridx)] .= true

    η
end

function ancestral_mask!(η, genealogy, ωs::Set{Ω}; wipe = true)
    wipe && wipe!(η)

    for ω ∈ ωs
        ancestral_mask!(η, genealogy, ω, wipe = false)
    end

    η
end

ancestral_mask(genealogy, ω) =
    ancestral_mask!(Sequence(falses(nmarkers(genealogy))), genealogy, ω,
                    wipe = false)

ancestral_mask!(η, ωs, genealogy, x::Union{VertexType, EdgeType}; wipe = true) =
    ancestral_mask!(η, genealogy, ancestral_intervals!(ωs, genealogy, x), wipe = wipe)

ancestral_mask(genealogy, x::Union{VertexType, EdgeType}) =
    ancestral_mask!(Sequence(undef, nmarkers(genealogy)), Set{Ω}(), genealogy, x)

function ancestral_mask!(η, genealogy, x::AbstractFloat; wipe = true)
    wipe && wipe!(η)
    η[postoidx(genealogy, x)] = true
    η
end

wipe!(η) = η.data.chunks .⊻= η.data.chunks

###################
# Pretty printing #
###################

function show(io::IO, ::MIME"text/plain", genealogy::AbstractGenealogy)
    if isempty(genealogy)
        println(io, "Empty " * describe(genealogy))
        return
    end

    println(io, describe(genealogy) * ":")
    println(io, nleaves(genealogy), " leaves, ",
            nmarkers(genealogy), " markers")
    print(io, "tMRCA: ", tmrca(genealogy))
end

function show(io::IO, genealogy::AbstractGenealogy)
    if isempty(genealogy)
        println(io, "Empty " * describe(genealogy, false))
        return
    end

    print(io, describe(genealogy, false) * "(")
    print(io, nleaves(genealogy))
    print(io, " leaves, ")
    print(io, nmarkers(genealogy))
    print(io, " markers)")
end

###########################
# AbstractGraph Interface #
###########################

## See https://juliagraphs.org/Graphs.jl/stable/ecosystem/interface/

for fun ∈ [:edges, :vertices, :ne, :nv]
    @eval function $fun(genealogy::AbstractGenealogy)
        $fun(graph(genealogy))
    end
end

for (fun, ret) ∈ Dict(:eltype => VertexType,
                      :edgetype => EdgeType,
                      :is_directed => :true)
    @eval begin
        @generated $fun(::Type{<:AbstractGenealogy}) = $ret
        @generated $fun(::AbstractGenealogy) = $ret
    end
end

let genealogy_edge = Expr(:(::), :e, EdgeType),
    genealogy_vertex = Expr(:(::), :v, VertexType)

    for (fun, a) ∈ Dict(:has_edge => genealogy_edge,
                        :has_vertex => genealogy_vertex,
                        :inneighbors => genealogy_vertex,
                        :outneighbors => genealogy_vertex)
        varname = first(a.args)
        @eval function $fun(genealogy::AbstractGenealogy, $a)
            $fun(graph(genealogy), $varname)
        end
    end
end

############
# Plotting #
############

"""
    maxdepth(genealogy, v)

Compute the depth of a vertex, that is the number of egdes between it and the
genealogy's mrca.
"""
function maxdepth(genealogy, v, depth = 0)
    _dads = dads(genealogy, v)
    isempty(_dads) && return depth

    mapreduce(d -> maxdepth(genealogy, d, depth + 1), max, _dads)
end

"""
    GenLayout(leaveslayer)

Layout function designed to plot genealogies. Leaves are all placed on the
same layer.
"""
GenLayout(leaveslayer, leaves) = function (genealogy_graph)
    xs, ys, _ = solve_positions(Zarate(),
                                genealogy_graph,
                                force_layer = Pair.(leaves, leaveslayer))

    ## Rotate by -π/2.
    Point.(zip(ys, -xs))
end

function graphplot(genealogy::AbstractGenealogy, ω;
                   wild_color = :blue,
                   derived_color = :red,
                   arrow_show = false,
                   edge_color = :gray,
                   edge_width = 3,
                   layout = nothing,
                   attributes...)
    if isnothing(layout)
        lastlayer = maximum(v -> maxdepth(genealogy, v), leaves(genealogy)) + 1
        layout = GenLayout(lastlayer, leaves(genealogy))
    end

    ## Color of the vertices.
    mask = fill(ancestral_mask(genealogy, ω), nv(genealogy))
    node_color = ifelse.(any.(sequences(genealogy) .& mask),
                         derived_color, wild_color)

    ## Hide non ancestral edges.
    ewidth = DefaultDict{EdgeType, Int}(edge_width)
    for e ∈ edges(genealogy)
        isdisjoint(ancestral_intervals(genealogy, e), ω) || continue
        ewidth[e] = 0
    end

    graphplot(graph(genealogy),
              layout = layout,
              ilabels = string.(range(1, nv(genealogy))),
              ilabels_color = :white,
              edge_color = edge_color,
              edge_width = ewidth,
              arrow_show = arrow_show,
              node_color = node_color,
              attributes...)
end

graphplot(genealogy::AbstractGenealogy; attributes...) =
    graphplot(genealogy, Ω(0, ∞); attributes...)

##################
# Common Methods #
##################

export isleaf, isroot, isivertex
"""
    isleaf(genealogy, v)
    isroot(genealogy, v)
    isivertex(genealogy, v)

True if `v` is a leaf/root in the given genealogy.
"""
function isleaf end,
function isroot end,
function isivertex end

isleaf(genealogy, v) = (iszero ∘ outdegree)(genealogy, v)

isroot(genealogy, v) = (iszero ∘ indegree)(genealogy, v)

isivertex(genealogy, v) = !isleaf(genealogy, v)

export nleaves, nivertices
"""
    nleaves(genealogy)
    nivertices(genealogy)

Number of leaves/internal vertices in a genealogy.

See also [`leaves`](@ref) and [`ivertices`](@ref) for an iterable over the
leaves/vertices.
"""
function nleaves end,
function nivertices end

nleaves(genealogy) = (length ∘ leaves)(genealogy)

nivertices(genealogy) = nv(genealogy) - nleaves(genealogy)

export sequence
"""
    sequence(genealogy, v)

Sequence of a genealogy associatex with a vertex.

See also [`sequences`](@ref) to get all the sequences of a genealogy.
"""
sequence(genealogy, v::VertexType) = getindex(sequences(genealogy), v)

export nmarkers
"""
    nmarkers(genealogy)

Number of markers in the sequences of a genealogy.
"""
nmarkers(genealogy) = (length ∘ first)(sequences(genealogy))

export branchlength
"""
    branchlength(genealogy[, e])

Total branch length of a genealogy. If an edge is specified, returns the length
of that edge.
"""
function branchlength end

branchlength(genealogy, e::EdgeType) =
    latitude(genealogy, src(e)) - latitude(genealogy, dst(e))

branchlength(genealogy) = mapreduce(e -> branchlength(genealogy, e), +, edges(genealogy))

branchlength(genealogy, ω) = mapreduce(e -> branchlength(genealogy, e), +,
                                       edges_interval(genealogy, ω))

export tmrca
function tmrca(genealogy)
    _mrca = mrca(genealogy)
    iszero(_mrca) && return zero(Float64)

    latitude(genealogy, _mrca)
end

tmrca(genealogy, vs) = latitude(genealogy, mrca(genealogy, vs))

export dads, dads!, children, children!
export descendants!, descendants, ancestors!, ancestors
"""
    dads!(buf, genealogy, v, ω)
    dads(genealogy, v[, ω])
    ancestors!(buf, genealogy, v[, ω])
    ancestors(genealogy, v[, ω])
    children!(buf, genealogy, v, ω)
    children(genealogy, v[, ω])
    descendants!(buf, genealogy, v[, ω])
    descendants(genealogy, v[, ω])

Parents/children/descendants of a vertex. ω can be either
- a number in [0, 1] representing a position;
- an Ω representing an interval of positions;
- a set of Ωs representing multiple interval of positions.

The following rules are used to decide if an edge `e` is ancestral:
- If ω is a number, the ancestral interval of `e` must cover ω.
- If ω is an Ω or a set of Ωs, the intersection of the ancestral
  interval of `e` with ω must be non-empty.
"""
function dads end, function dads! end,
function ancestors end, function ancestors! end,
function children end, function children! end,
function descendants end, function descendants! end

for (fun, graph_method) ∈ Dict(:dads => :inneighbors, :children => :outneighbors)
    @eval $fun(genealogy, v) = $graph_method(genealogy, v)
end

let funtransorder = Dict(:dads => (:ancestors, (x, y) -> (x, y)),
                         :children => (:descendants, (x, y) -> (y, x))),

    typesandfun = ((:Real, in), (:Ω, !isdisjoint), (:(Set{Ω}), !isdisjoint))
    for ((fun, (transfun, order)), (Argtype, testfun)) ∈
        Iterators.product(funtransorder, typesandfun)

        ## Parents & children
        fun! = Symbol(string(fun) * '!')

        @eval function $fun!(buf, genealogy, v, ω::$Argtype)
            resize!(buf, 2)
            ptr = firstindex(buf)

            @inbounds for u ∈ $fun(genealogy, v)
                ωs = ancestral_intervals(genealogy, Edge($order(u, v)))
                $testfun(ω, ωs) || continue

                buf[ptr] = u
                ptr += 1
            end

            resize!(buf, ptr - 1)
        end

        @eval $fun(genealogy, v::T, ω::$Argtype) where T =
            $fun!(sizehint!(Vector{T}(undef, 2), 2), genealogy, v, ω)

        ## Ancestors & descendants
        transfun! = Symbol(string(transfun) * '!')

        @eval function $transfun!(buf, genealogy, v, ω::$Argtype)
            funbuf = sizehint!(Vector{VertexType}(undef, 2), 2)
            writeptr = readptr = firstindex(buf)
            @inbounds for u ∈ $fun!(funbuf, genealogy, v, ω)
                buf[writeptr] = u
                writeptr += 1
            end

            @inbounds while readptr < writeptr
                v = buf[readptr]
                readptr += 1
                (isleaf(genealogy, v) || isroot(genealogy, v)) && continue

                resize!(funbuf, 2)
                for u ∈ $fun!(funbuf, genealogy, v, ω)
                    u ∈ view(buf, 1:(writeptr-1)) && continue
                    buf[writeptr] = u
                    writeptr += 1
    end
end

            resize!(buf, writeptr - 1)
        end

        @eval function $transfun(genealogy, v, ω::$Argtype)
            $transfun!(Vector{VertexType}(undef, nv(genealogy) - 1),
                       genealogy, v, ω)
        end

        @eval function $transfun!(buf, genealogy, v)
            writeptr = readptr = firstindex(buf)
            @inbounds for u ∈ $fun(genealogy, v)
                buf[writeptr] = u
                writeptr += 1
            end

            @inbounds while readptr < writeptr
                v = buf[readptr]
                readptr += 1
                (isleaf(genealogy, v) || isroot(genealogy, v)) && continue

                for u ∈ $fun(genealogy, v)
                    u ∈ view(buf, 1:(writeptr-1)) && continue
                    buf[writeptr] = u
                    writeptr += 1
                end
            end

            resize!(buf, writeptr - 1)
end

        @eval function $transfun(genealogy, v)
            $transfun!(Vector{VertexType}(undef, nv(genealogy) - 1),
                       genealogy, v)
        end
    end
end

export nmutations
"""
    nmutations(genealogy[, e])

Number of mutation on a genealogy. If an edge is specified, return only the
number of mutations on that edge.
"""
function nmutations end

function nmutations!(mask, ωs, genealogy, e)
    ret = zero(Int)
    ancestral_mask!(mask, ωs, genealogy, e)

    nchunks = (length(positions(genealogy)) - 1) ÷ blocksize(Sequence) + 1
    η1, η2 = sequences(genealogy, e)
    @inbounds for k ∈ range(1, length = nchunks)
        ret += count_ones((η1.data.chunks[k] ⊻ η2.data.chunks[k]) & mask.data.chunks[k])
    end

    ret
end

function nmutations(genealogy)
    mask = Sequence(undef, nmarkers(genealogy))
    ωs = Set{Ω}()
    mapreduce(e -> nmutations!(mask, ωs, genealogy, e), +, edges(genealogy),
              init = zero(Int))
end

## Edge functions.
for fun ∈ [:branchlength, :nmutations]
    @eval function $fun(genealogy, σ, δ)
        $fun(genealogy, Edge(σ, δ))
    end
end

"""
    mutable struct EdgeIntervalIter{G}

Iterate over the edges of a genealogy that have ancestral material in a given
interval.
"""
struct EdgeIntervalIter{G, O}
    genealogy::G
    ω::O
    vstack::Stack{VertexType}
    otheredge::Base.RefValue{EdgeType}
    visited::RBTree{VertexType}
    dads_buf::Vector{VertexType}
end

@generated IteratorSize(::EdgeIntervalIter) = Base.SizeUnknown()

@generated eltype(::EdgeIntervalIter) = EdgeType

function edges_interval(genealogy, ωs = Ω(0, ∞))
    vstack = Stack{VertexType}(ceil(Int, log(nv(genealogy))))
    root = argmax(latitudes(genealogy)) + nleaves(genealogy)
    push!(vstack, children(genealogy, root)...)

    EdgeIntervalIter(genealogy, ωs, vstack, Ref(Edge(0 => 0)),
                     RBTree{VertexType}(),
                     sizehint!(Vector{VertexType}(undef, 2), 2))
end

function iterate(eit::EdgeIntervalIter, state = 0)
    otheredge = eit.otheredge
    if otheredge[] != Edge(0 => 0)
        ret = otheredge[]
        otheredge[] = Edge(0 => 0)
        return ret, state + 1
    end

    vstack = eit.vstack
    isempty(vstack) && return nothing

    genealogy = eit.genealogy
    ω = eit.ω
    visited = eit.visited

    v = pop!(vstack)
    for child ∈ children(genealogy, v)
        if child ∉ visited
            push!(vstack, child)
            indegree(genealogy, child) > 1 && push!(visited, child)
        end
    end

    _dads = eit.dads_buf
    dads!(_dads, genealogy, v, ω)
    e = Edge(first(_dads), v)
    if length(_dads) > 1
        otheredge[] = e
        return Edge(last(_dads), v), state + 1
    else
        return e, state + 1
    end
end
