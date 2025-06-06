using Graphs

using Bumper: UnsafeArrays
using Graphs: AbstractSimpleGraph

import Graphs: edges, vertices, ne, nv,
               eltype, edgetype, is_directed,
               has_edge, has_vertex, inneighbors, outneighbors

using LayeredLayouts

using GeometryBasics: Point

using DataStructures: DefaultDict

import Base: IteratorSize, eltype, length

using NetworkLayout

using UnicodePlots: histogram

export AbstractGenealogy
"""
    $(TYPEDEF)

Abstract type for genealogies.

Implements [the Graphs.jl interface](https://juliagraphs.org/Graphs.jl/stable/ecosystem/interface/)
so that subtypes implementing the AbstractGenealogy interface can additionally
be treated as regular graphs.

[Tree](@ref) & [Arg](@ref) are subtypes that implement the AbstractGenealogy
interface.
"""
abstract type AbstractGenealogy <: AbstractSimpleGraph{VertexType} end

#############
# Interface #
#############

export graph
"""
    $(FUNCTIONNAME)(genealogy)

Return the underlying graph of a genealogy.

# Methods
$(METHODLIST)
"""
function graph end

"""
    $(FUNCTIONNAME)(genealogy, long = true)

Return a string containing a long or short description of a genealogy.
Used internally by pretty printing functions.

# Methods
$(METHODLIST)
"""
function describe end

export latitudes
"""
    $(FUNCTIONNAME)(genealogy[, vs])

Return the latitudes of (a subset of) the internal vertices of a genealogy.

See also [`latitude`](@ref) to get the latitude of a single vertex.

# Implementation
A default implementation for `latitudes(::AbstractGenealogy, ::Any)` is
available; only `latitudes(::T)` is required.

# Methods
$(METHODLIST)
"""
function latitudes end

latitudes(genealogy, vs) = getindex(latitudes(genealogy), vs)

export latitude
"""
    $(FUNCTIONNAME)(genealogy, v)

Latitude of a vertex.

See also [`latitudes`](@ref) to get the latitudes of internal vertices stored
in a genealogy.

# Methods
$(METHODLIST)
"""
function latitude end

export leaves
"""
    $(SIGNATURES)

Return an interable containing the leaves of a genealogy.

See also [`ivertices`](@ref) for internal vertices and [`nleaves`](@ref) for
the number of leaves.

# Implementation
Default implementations assume that the first `nleaves(genealogy)` vertices
are the leaves of the genealogy. If this is the case for your type, you do not
need to implement this method.
"""
leaves(genealogy) = Base.OneTo(nleaves(genealogy))

export ivertices
"""
    $(SIGNATURES)

Return an interable containing the internal vertices of a genealogy.

See also [`leaves`](@ref) for leaves and [`nivertices`](@ref) for the number
of leaves and internal vertices.

# Implementation
Default implementations assume that the first `nleaves(genealogy)` vertices
are the leaves of the genealogy. If this is the case for your type, you do not
need to implement this method.
"""
ivertices(genealogy) = UnitRange(nleaves(genealogy) + 1, nv(genealogy))

export sequences
"""
    $(FUNCTIONNAME)(genealogy[, e])
    $(FUNCTIONNAME)(genealogy[, vs])

Sequences of a genealogy. If an edge is specified, return the sequences
associated with the vertices incident to that edge. If an iterable of
vertices is specified, return the sequences associated with these vertices.

See also [`sequence`](@ref) to get the sequence associated with a vertex.

# Implementation
Custom types only need to implement `sequences(::T)`.

# Methods
$(METHODLIST)
"""
function sequences end

sequences(genealogy, vs) = Iterators.map(v -> sequence(genealogy, v), vs)

sequences(genealogy, e::Edge) = (sequence(genealogy, src(e)),
                                     sequence(genealogy, dst(e)))

export mrca
"""
    $(FUNCTIONNAME)(genealogy[, vs = leaves(genealogy), ωs = Ω(0, ∞)])

Most recent common ancestor of a set of vertices.

See also [`tmrca`](@ref) for the time to the most recent common ancestor.

# Methods
$(METHODLIST)
"""
function mrca end

export dens
"""
    $(FUNCTIONNAME)(genealogy; logscale = false)

Likelihood of a genealogy.

# Implementation
Only mandatory if copulas are to be fitted on the genealogy.

# Methods
$(METHODLIST)
"""
function dens end

export positions
"""
    $(FUNCTIONNAME)(genealogy)

Positions of the markers.

# Methods
$(METHODLIST)
"""
function positions end

export ancestral_intervals
"""
    $(FUNCTIONNAME)(genealogy, x)

Interval for which x is ancestral. Default implementation assumes that anything
is ancestral for [0, ∞).

See also [`ancestral_intervals!`](@ref) for a non-allocating alternative.

# Methods
$(METHODLIST)
"""
function ancestral_intervals end

@generated ancestral_intervals(genealogy::Any, x::Any) =
    ancestral_intervals!(AIsType(), genealogy, x)

export ancestral_intervals!
"""
    $(FUNCTIONNAME)(ωs, genealogy, x)

Interval for which x is ancestral. Default implementation assumes that anything
is ancestral for [0, ∞).

See also [`ancestral_intervals`](@ref).

# Methods
$(METHODLIST)
"""
function ancestral_intervals! end

function ancestral_intervals!(ωs, ::Any, ::Any)
    empty!(ωs)
    push!(ωs, Ω(0, ∞))
end

## Recombinations ##

export nrecombinations
"""
    $(FUNCTIONNAME)(genealogy)

Number of recombinations.

# Implementation
Default implementation returns 0.

# Methods
$(METHODLIST)
"""
function nrecombinations end

@generated nrecombinations(::Any) = zero(Int)

export recombinations
"""
    $(FUNCTIONNAME)(genealogy)

Iterator over the recombination vertices of a genealogy.

# Implementation
Default implementation returns an empty iterator.

# Methods
$(METHODLIST)
"""
function recombinations end

@generated recombinations(::Any) = StepRange{Int, Int}(0, 1, 0)

export isrecombination
"""
    $(FUNCTIONNAME)(genealogy, v)

Returns true if vertex `v` is a recombination for `genealogy`.

Default implementation always returns `false`.

# Methods
$(METHODLIST)
"""
function isrecombination end

@generated isrecombination(::Any, ::Any) = false

"""
    $(FUNCTIONNAME)(genealogy)

Layout function for genealogy plotting.

# Implementation
Defaults to `Spring()`.

# Methods
$(METHODLIST)

--*Internal*--
"""
function plot_layout end

plot_layout(::AbstractGenealogy) = Spring()

export plot_genealogy
"""
    $(FUNCTIONNAME)(genealogy[, ω]; kwargs...)

Plot a genealogy.

Only edges and vertices ancestral for `ω` are plotted.

Implemented as an extension to `Moonshine.jl`. To use this method, you have to
import [`GraphMakie`](@ref) (and a
[Makie](https://github.com/MakieOrg/Makie.jl) backend).

See also [`plot_layout`](@ref).

# Arguments
* `wild_color` (`:blue`): color of wild vertices, that is those having only wild
  markers in `ω`.
* `derived_color` (`:red`): color of derived (non-wild) vertices
* `arrow_show` (`false`): whether or not to draw arrows
* `edge_color` (`:gray`): color of edges
* `edge_width` (`3`): width of edges
* `layout` (`plot_layout(genealogy)`): layout function
* `attributes...`: attributes passed directly to
  [`GraphMakie.graphplot`](@extref)
"""
function plot_genealogy end

"""
    $(FUNCTIONNAME)(genealogy)

Maximum possible number of parents for a vertex. See also [`maxchildren`](@ref).

# Implementation
Must be a generated function.

# Methods
$(METHODLIST)

--*Internal*--
"""
function maxdads end

"""
    $(FUNCTIONNAME)(genealogy)

Maximum possible number of children for a vertex. See also [`maxdads`](@ref).

# Implementation
Must be a generated function.

# Methods
$(METHODLIST)

--*Internal*--
"""
function maxchildren end

#############
# Utilities #
#############

export idxtopos
"""
    $(FUNCTIONNAME)(genealogy, idx)

Return the position of the marker given its index.

# Methods
$(METHODLIST)
"""
function idxtopos end

export postoidx
"""
    $(FUNCTIONNAME)(genealogy, pos)

Return the largest marker's index that is at a position lesser than the one
given.

# Methods
$(METHODLIST)
"""
function postoidx end

for (f, arg) ∈ (:idxtopos => :idx, :postoidx => :pos)
    @eval $f(genealogy::AbstractGenealogy, $arg) = $f(sam(genealogy), $arg)
end

export ancestral_mask
"""
    $(FUNCTIONNAME)(reference, x; ωs_buf = Set{Ω}())

Mask non ancestral positions to 0. If `wipe = true`, all markers in `η` will be
initialized to 0.

# Methods
$(METHODLIST)
"""
function ancestral_mask end

ancestral_mask(genealogy::AbstractGenealogy, x) =
    ancestral_mask!(Sequence(falses(nmarkers(genealogy))), genealogy, x,
                    wipe = false)

export ancestral_mask!
"""
    $(FUNCTIONNAME)(η, reference, x; ωs_buf = Set{Ω}(), wipe = true)

Mask non ancestral positions to 0. If `wipe = true`, all markers in `η` will be
initialized to 0.

# Methods
$(METHODLIST)
"""
function ancestral_mask! end

ancestral_mask!(η, genealogy::AbstractGenealogy, x; wipe = true) =
    ancestral_mask!(η, sam(genealogy), x, wipe = wipe)

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
                      :edgetype => Edge{VertexType},
                      :is_directed => :true)
    @eval begin
        @generated $fun(::Type{<:AbstractGenealogy}) = $ret
        @generated $fun(::AbstractGenealogy) = $ret
    end
end

let genealogy_edge = Expr(:(::), :e, Edge),
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

## TODO: `yflip` doesn't work
export plot_latitudes
"""
    $(SIGNATURES)

Unicode histogram of a genealogy's latitudes. Additional keywords arguments are
passed directly to [`UnicodePlots.histogram`](https://github.com/JuliaPlots/UnicodePlots.jl#histogram).

See also [`latitudes`](@ref), [`latitude`](@ref)
"""
plot_latitudes(genealogy::AbstractGenealogy; kwargs...) =
    histogram(latitudes(genealogy),
              yflip = true,
              title = "Vertices' Latitudes",
              xlabel = "";
              kwargs...)

##################
# Common Methods #
##################

export isleaf
"""
    $(FUNCTIONNAME)(genealogy, v)

True if `v` is a leaf.

See also [`isivertex`](@ref), [`isroot`](@ref).
"""
function isleaf end

export isroot
"""
    $(FUNCTIONNAME)(genealogy, v)

True if `v` is the root.

See also [`isleaf`](@ref), [`isivertex`](@ref).
"""
function isroot end

for (fun, degfun) ∈ Dict(:isleaf => :outdegree, :isroot => :indegree)
    @eval $fun(genealogy, v) = (iszero ∘ $degfun)(genealogy, v)
end

export isivertex
"""
    $(SIGNATURES)

True if `v` is an internal vertex.

See also [`isleaf`](@ref), [`isroot`](@ref).
"""
isivertex(genealogy, v) = !isleaf(genealogy, v)

export nleaves
"""
    $(FUNCTIONNAME)(genealogy)

Number of leaves in a genealogy.

See also [`leaves`](@ref) for an iterable over leaves and [`nivertices`](@ref)
for an internal vertices counterpart.

# Methods
$(METHODLIST)
"""
function nleaves end

nleaves(genealogy) = (length ∘ leaves)(genealogy)

export nivertices
"""
    $(SIGNATURES)

Number of internal vertices in a genealogy.

See also [`ivertices`](@ref) for an iterable over internal vertices and
[`nleaves`](@ref) for a leaves counterpart.
"""
nivertices(genealogy) = nv(genealogy) - nleaves(genealogy)

export sequence
"""
    $(SIGNATURES)

Sequence of a genealogy associated with a vertex.

See also [`sequences`](@ref) to get all the sequences of a genealogy.
"""
sequence(genealogy, v) = getindex(sequences(genealogy), v)

export nmarkers
"""
    $(FUNCTIONNAME)(genealogy[, ωs])

Number of markers in the sequences of a genealogy. If an interval `ωs` is
specified, returns the number of markers contained in that interval.

# Methods
$(METHODLIST)
"""
function nmarkers end

nmarkers(genealogy, ωs...) = nmarkers(sam(genealogy), ωs...)

export branchlength
"""
    $(FUNCTIONNAME)(genealogy)
    $(FUNCTIONNAME)(genealogy, ωs)
    $(FUNCTIONNAME)(genealogy, e)

Total branch length of a genealogy. If an interval is specified, returns the
branch length of the associated marginal genealogy. If an edge is specified,
returns the length of that edge.

# Methods
$(METHODLIST)
"""
function branchlength end

branchlength(genealogy, e::Edge) =
    latitude(genealogy, src(e)) - latitude(genealogy, dst(e))

branchlength(genealogy) = mapreduce(e -> branchlength(genealogy, e), +, edges(genealogy))

function branchlength(genealogy, ωs; buffer = buffer)
    @no_escape buffer begin
        store = @alloc(Edge{VertexType}, nleaves(genealogy))
        visited = @alloc(Bool, nrecombinations(arg))
        ret = sum(e -> branchlength(genealogy, e),
                  edges_interval(genealogy, ωs, store, visited))
    end

    ret
end

export tmrca
"""
    $(FUNCTIONNAME)(genealogy[, vs])

Time to the most recent common ancestor (MRCA) of (a subset of) vertices.

# Methods
$(METHODLIST)
"""
function tmrca end

function tmrca(genealogy)
    _mrca = mrca(genealogy)
    iszero(_mrca) && return zero(Float64)

    latitude(genealogy, _mrca)
end

tmrca(genealogy, vs) = latitude(genealogy, mrca(genealogy, vs))

export dads
"""
    $(FUNCTIONNAME)(genealogy, v[, ωs])

Parents of a vertex, optionally restricted to a marginal genealogy. If you know
in advance that `v` has a single dad, use [̀̀`dad`](@ref) instead.

The following rules are used to decide if an edge `e` is ancestral:
- If ωs is a number, the ancestral interval of `e` must cover ωs.
- If ωs is an Ω or a set of Ωs, the intersection of the ancestral
  interval of `e` with ωs must be non-empty.

!!! danger
    Return a **reference** to the underlying adjacency lists. No touchy!

See also [`child`](@ref), [`children`](@ref), [`sibling`](@ref),
[`siblings`](@ref), [`descendants`](@ref) and [`ancestors`](@ref).

# Methods
$(METHODLIST)
"""
function dads end

export children
"""
    $(FUNCTIONNAME)(genealogy, v[, ωs])

Children of a vertex, optionally restricted to a marginal genealogy. If you know
in advance that `v` has a single child, use [̀̀`child`](@ref) instead.


The following rules are used to decide if an edge `e` is ancestral:
- If ωs is a number, the ancestral interval of `e` must cover ωs.
- If ωs is an Ω or a set of Ωs, the intersection of the ancestral
  interval of `e` with ωs must be non-empty.

!!! danger
    Return a **reference** to the underlying adjacency lists. No touchy!

See also [`dad`](@ref), [`dads`](@ref), [`sibling`](@ref),
[`siblings`](@ref), [`descendants`](@ref) and [`ancestors`](@ref).

# Methods
$(METHODLIST)
"""
function children end

export children3
"""
    $(FUNCTIONNAME)(genealogy, v, args...)

Childrens of a vertex, ignoring marginally degree 2 vertices.

All arguments are passed directly to [`children`](@ref).

See also [`dads3`](@ref).

!!! danger
    Return a **reference** to the underlying adjacency lists. No touchy!

# Methods
$(METHODLIST)
"""
function children3 end

export dads3
"""
    $(FUNCTIONNAME)(genealogy, v, args...)

Parents of a vertex, ignoring marginally degree 2 vertices.

All arguments are passed directly to [`dads`](@ref).

See also [`children3`](@ref).

!!! danger
    Return a **reference** to the underlying adjacency lists. No touchy!

# Methods
$(METHODLIST)
"""
function dads3 end

for (fun, list) ∈ Dict( :dads => Meta.quot(:badjlist), :children => Meta.quot(:fadjlist))
    @eval $fun(genealogy, v) = getfield(graph(genealogy), $list)[v]

    fun3 = Symbol(string(fun) * '3')
    @eval function $fun3(genealogy, v, args...)
        neig = $fun(genealogy, v, args...)
        @inbounds while (isone ∘ length)(neig)
            neig = $fun(genealogy, first(neig), args...)
        end

        neig
    end
end

export ancestors
"""
    $(FUNCTIONNAME)(genealogy, v[, ωs])
    $(FUNCTIONNAME)(genealogy, v[, ωs], buf_ptr)

Ancestors of a vertex, optionally restricted to a marginal genealogy.

A pointer to some kind of buffer (an array for instance) can be provided to
avoid allocation. In that case, an `UnsafeArray` wrapped around it will be
returned.

The following rules are used to decide if an edge `e` is ancestral:
- If ωs is a number, the ancestral interval of `e` must cover ωs.
- If ωs is an Ω or a set of Ωs, the intersection of the ancestral
  interval of `e` with ωs must be non-empty.

See also [`child`](@ref), [`dad`](@ref), [`children`](@ref), [`dads`](@ref),
[`sibling`](@ref), [`siblings`](@ref) and [`descendants`](@ref).

# Methods
$(METHODLIST)
"""
function ancestors end

export descendants
"""
    $(FUNCTIONNAME)(genealogy, v[, ωs])
    $(FUNCTIONNAME)(genealogy, v[, ωs], buf_ptr)

Descendants of a vertex, optionally restricted to a marginal genealogy.

A pointer to some kind of buffer (an array for instance) can be provided to
avoid allocation. In that case, an `UnsafeArray` wrapped around it will be
returned.

The following rules are used to decide if an edge `e` is ancestral:
- If ωs is a number, the ancestral interval of `e` must cover ωs.
- If ωs is an Ω or a set of Ωs, the intersection of the ancestral
  interval of `e` with ωs must be non-empty.

See also [`child`](@ref), [`dad`](@ref), [`children`](@ref), [`dads`](@ref),
[`sibling`](@ref), [`siblings`](@ref) and [`ancestors`](@ref).

# Methods
$(METHODLIST)
"""
function descendants end

let funtrans = Dict(:dads => :ancestors, :children => :descendants)
    for (fun, transfun) ∈ funtrans
        transfun! = Symbol(string(transfun) * '!')

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

        ## Ancestors & descendants
        @eval function $transfun!(buf_ptr, genealogy, v, ω)
            head = tail = zero(Int)
            @inbounds for u ∈ $fun(genealogy, v, ω)
                head += 1
                unsafe_store!(buf_ptr, u, head)
            end

            @inbounds while tail < head
                tail += 1
                v = unsafe_load(buf_ptr, tail)
                (isleaf(genealogy, v) || isroot(genealogy, v)) && continue

                for u ∈ $fun(genealogy, v, ω)
                    for k ∈ 1:head
                        u == unsafe_load(buf_ptr, k) && @goto skip
                    end

                    head += 1
                    unsafe_store!(buf_ptr, u, head)

                    @label skip
                end
            end

            UnsafeArray{eltype(buf_ptr), 1}(buf_ptr, (head,))
        end

        @eval function $transfun(genealogy, v, ω)
            $transfun!(pointer(Vector{VertexType}(undef, nv(genealogy) - 1)),
                       genealogy, v, ω)
        end
    end
end

export siblings
"""
    $(FUNCTIONNAME)(genealogy, v[, ωs])
    $(FUNCTIONNAME)(genealogy, v[, ωs], buf_ptr)

Siblings of a vertex, that is the other vertices in the genealogy that share at
least one parent, optionally restricted to a marginal genealogy.

A pointer to some kind of buffer (an array for instance) can be provided to
avoid allocation. In that case, an `UnsafeArray` wrapped around it will be
returned.

If you know in advance that `v` has a single sibling, you can use
[`sibling`](@ref) instead.

See also [`child`](@ref), [`dad`](@ref), [`children`](@ref), [`dads`](@ref),
[`descendants`](@ref) and [`ancestors`](@ref).

# Methods
$(METHODLIST)
"""
function siblings end

function siblings!(x::Ptr, genealogy, v, args...)
    len = 0

    for _dad ∈ dads(genealogy, v, args...)
        for _child ∈ children(genealogy, _dad, args...)
            _child == v && continue
            len += 1
            unsafe_store!(x, _child, len)
        end
    end

    UnsafeArray{VertexType, 1}(x, (len,))
end

function siblings!(x::AbstractArray, genealogy, v, args...)
    len = 0

    for _dad ∈ dads(genealogy, v, args...)
        for _child ∈ children(genealogy, _dad, args...)
            _child == v && continue
            len += 1
            x[len] = _child
        end
    end

    resize!(x, len)
end

function siblings(genealogy, v, args...)
    x = Vector{VertexType}(undef, 2)
    siblings!(x, genealogy, v, args...)
end

export sibling
"""
    $(SIGNATURES)

Sibling of a vertex, that is the other vertex in the genealogy that have the
same parent. It only makes sense to use this method if you know `v` has a
single sibling. Otherwise use [`siblings`](@ref).

See also [`child`](@ref), [`dad`](@ref), [`children`](@ref), [`dads`](@ref),
[`descendants`](@ref) and [`ancestors`](@ref).
"""
function sibling(genealogy, v)
    for _child ∈ children(genealogy, dad(genealogy, v))
        _child == v && continue
        return _child
    end
    zero(VertexType)
end

export dad
"""
    $(FUNCTIONNAME)(genealogy, v)

Return the parent of a vertex or 0 if none. It only makes sense to
use this method if you know `v` has a single parent. Otherwise use
[`dads`](@ref).

See also [`child`](@ref), [`descendants`](@ref), [`ancestors`](@ref) and
[̀`siblings`](@ref).
"""
function dad end

export child
"""
    $(FUNCTIONNAME)(genealogy, v)

Return the child of a vertex or 0 if none. It only makes sense to
use this method if you know `v` has a single child. Otherwise use
[`children`](@ref).

See also [`dad`](@ref), [`ancestors`](@ref), [`descendants`](@ref) and
[̀`siblings`](@ref).
"""
function child end

for (fun, nei) ∈ Dict(:dad => :inneighbors, :child => :outneighbors)
    @eval function $fun(genealogy, v)
        neig = $nei(genealogy, v)
        isempty(neig) && return zero(VertexType)
        first(neig)
    end
end

export nmutations
"""
    $(FUNCTIONNAME)(genealogy[, e])

Number of mutation on a genealogy. If an edge is specified, return the
number of mutations on that edge.

See also [`nmutations!`](@ref) for a non-allocating alternative.

# Methods
$(METHODLIST)
"""
function nmutations end

export nmutations!
"""
    $(SIGNATURES)

Number of mutation on a genealogy. If an edge is specified, return the
number of mutations on that edge.

See also [`nmutations`](@ref).
"""
function nmutations!(mask, genealogy, e)
    ret = zero(Int)
    ancestral_mask!(mask, genealogy, e)

    nchunks = (length(positions(genealogy)) - 1) ÷ blocksize(Sequence) + 1
    η1, η2 = sequences(genealogy, e)
    @inbounds for k ∈ range(1, length = nchunks)
        ret += count_ones((η1.data.chunks[k] ⊻ η2.data.chunks[k]) & mask.data.chunks[k])
    end

    ret
end

function nmutations(genealogy)
    mask = Sequence(undef, nmarkers(genealogy))
    mapreduce(e -> nmutations!(mask, genealogy, e), +, edges(genealogy),
              init = zero(Int))
end

function nmutations(genealogy, e)
    mask = Sequence(undef, nmarkers(genealogy))
    nmutations!(mask, genealogy, e)
end

## Edge functions.
for fun ∈ [:branchlength, :nmutations]
    @eval function $fun(genealogy, s, d)
        $fun(genealogy, Edge(s, d))
    end
end

"""
    $(TYPEDEF)

Flexible edge iterators that supports various constraints.

Possible constraints are any combination of the following:
* An interval of genetic positions. Any non-ancestral edge is ignored.
* A minimum latitude. Any edge under that latitude is ignored. An edge `e` is
  considered under a latitude `l` if `latitude(genealogy, dst(e)) < l`.
* A set of predicates. If any of the predicates evaluate to true for a given
  edge, that edge is ignored.

Do not construct directly, use [`edges_interval`](@ref) instead.

# Fields
$(TYPEDFIELDS)

# Constructors
$(METHODLIST)

--*Internal*--
"""
struct EdgesInterval{T, I, E}
    "Genealogy to iterate over"
    genealogy::T
    "Interval to consider"
    ωs::I
    "Edges buffer"
    buffer::CheapStack{E}
    "True is associated recombination vertex has been visited previously"
    visited::UnsafeArray{Bool, 1}
    "Only consider edges located above this latitude"
    min_latitude::Float64
    "Edges for which one of the predicates is true are blocked"
    block_predicates::Vector{FunctionWrapper{Bool, Tuple{E}}}
end

function EdgesInterval(genealogy::T, ωs::I, stack::CheapStack{E}, visited,
                       root = mrca(genealogy), min_latitude = zero(Float64),
                       block_predicates = []) where {T, I, E}
    fill!(visited, false)

    for d ∈ children(genealogy, root, ωs)
        e = Edge(root => d)
        any(p -> p(e), block_predicates) || push!(stack, e)
    end

    EdgesInterval{T, I, E}(genealogy, ωs, stack, visited,
                           convert(Float64, min_latitude), block_predicates)
end

EdgesInterval(genealogy::T, ωs::I, store::AbstractArray{E}, visited,
              root = mrca(genealogy), min_latitude = zero(Float64),
              block_predicates = []) where {T, I, E} =
    EdgesInterval(genealogy, ωs, CheapStack(store), visited,
                  root, min_latitude, block_predicates)

IteratorSize(::T) where T<:EdgesInterval = Base.SizeUnknown()
IteratorSize(::Type{<:EdgesInterval}) = Base.SizeUnknown()

eltype(::EdgesInterval) = Edge{VertexType}

_block(iter::EdgesInterval, e) = any(p -> p(e), iter.block_predicates)

function iterate(iter::EdgesInterval, state = 1)
    buffer = iter.buffer
    isempty(buffer) && return nothing

    genealogy = iter.genealogy
    ωs = iter.ωs
    visited = iter.visited
    min_latitude = iter.min_latitude
    n = nleaves(genealogy)
    block = e -> _block(iter, e)

    e = pop!(buffer)
    s = dst(e)
    latok = latitude(genealogy, s) >= min_latitude

    if isrecombination(genealogy, s)
        ridx = recidx(genealogy, s)
        visited[ridx] && return e, state + 1
        visited[ridx] = true

        ## Recombination vertex, no need to check ancestrality of downstream
        ## edge
        latok && push!(buffer, Edge(s => child(genealogy, s)))
    elseif latok
        for d ∈ children(genealogy, s, ωs)
            newe = Edge(s => d)
            block(newe) || push!(buffer, newe)
        end
    end

    e, state + 1
end

export edges_interval
"""
    $(FUNCTIONNAME)(genealogy, ωs)
    $(FUNCTIONNAME)(genealogy, ωs, buffer, visited, root = mrca(genealogy), min_latitude = zero(Float64), block_predicates = [])

Iterage over a genealogy's edges via [`EdgesInterval`](@ref).

`buffer` can be either a `CheapStack` or an `UnsafeArray` that will be used as
buffer for a newly constructed `CheapStack`.

# Methods
$(METHODLIST)
"""
function edges_interval end

edges_interval(genealogy, ωs, buffer, visited,
               root = mrca(genealogy), min_latitude = zero(Float64);
               block_predicates = []) =
    EdgesInterval(genealogy, ωs, buffer, visited, root, min_latitude, block_predicates)

function edges_interval(genealogy, ωs)
    ωs_e = AIsType()
    flt = function(e)
        ancestral_intervals!(ωs_e, genealogy, e)
        !isdisjoint(ωs_e, ωs)
    end

    Iterators.filter(flt, edges(genealogy))
end

export edgesmap
"""
    $(SIGNATURES)

Return a `Dict` that maps every edge of a genealogy to an integer in
`1:ne(genealogy)`.
"""
edgesmap(genealogy) = Dict(reverse.(enumerate(edges(genealogy))))

export nlive
"""
    $(SIGNATURES)

Number of live edges in a (marginal) genealogy at a given latitude.

`block_predicate` is passed directly to [`edges_interval`](@ref).

See also [`nlive!`](@ref).
"""
function nlive(genealogy, lat::Real, ωs;
               block_predicates = [], buffer = default_buffer())
    ## The grand MRCA is live forever
    lat > tmrca(genealogy) && return one(Int)

    live = zero(Int)

    @no_escape buffer begin
        store = @alloc(Edge{VertexType}, nleaves(genealogy) + nrecombinations(genealogy))
        visited = @alloc(Bool, nrecombinations(genealogy))
        @inbounds @simd ivdep for e ∈ edges_interval(genealogy, ωs, store, visited, mrca(genealogy),
                                                       lat, block_predicates = block_predicates)
            live += latitude(genealogy, dst(e)) <= lat <= latitude(genealogy, src(e))
        end
    end

    live
end

export nlive!
"""
    $(FUNCTIONNAME)(counts, genealogy, lats, ωs[, stack]; block_predicates = [], buffer = default_buffer())

Number of live edges in a (marginal) genealogy at a given latitude.

Counts are stored in `counts` which is filled with 0s initially. ̀`lats` must be
the same size as `counts`.

`block_predicate` and `stack` are passed directly to [`edges_interval`](@ref).

See also [`nlive`](@ref).

# Methods
$(METHODLIST)
"""
function nlive! end

function nlive!(counts, genealogy, lats::AbstractVector{<:Real}, ωs, stack;
                block_predicates = [], buffer = default_buffer())
    fill!(counts, 0)

    ## The grand MRCA is live forever
    first_above_tmrca_idx = findfirst(>(tmrca(genealogy)), lats)
    if isnothing(first_above_tmrca_idx)
        lastlat = lastindex(lats)
    else
        counts[first_above_tmrca_idx:end] .+= 1
        lastlat = first_above_tmrca_idx - 1
    end

    iszero(lastlat) && return counts

    @no_escape buffer begin
        visited = @alloc(Bool, nrecombinations(genealogy))
        @inbounds for e ∈ edges_interval(genealogy, ωs, stack, visited, mrca(genealogy), first(lats),
                                           block_predicates = block_predicates)
            @simd ivdep for k ∈ 1:lastlat
                counts[k] +=
                latitude(genealogy, dst(e)) <= lats[k] <= latitude(genealogy, src(e))
            end
        end
    end

    counts
end

nlive!(counts, genealogy, lats::AbstractVector{<:Real}, ωs;
       block_predicates = [], buffer = default_buffer()) =
    @no_escape buffer begin
        store = @alloc(Edge{VertexType}, nleaves(genealogy))
        stack = CheapStack(store)
        nlive!(counts, genealogy, lats, ωs, stack,
               block_predicates = block_predicates, buffer = buffer)
    end

export ismutation_edge
"""
    $(SIGNATURES)

Determines if an edge is a mutation edge for a given marker.
"""
function ismutation_edge(genealogy, e, idx)
    _chunkidx = chunkidx(Sequence, idx)
    hs = sequences(genealogy, e)
    chunk1, chunk2 = hs[1].data.chunks[_chunkidx], hs[2].data.chunks[_chunkidx]

    mask = one(UInt64) << (idxinchunk(Sequence, idx) - 1)
    !iszero((chunk1 ⊻ chunk2) & mask)
end
