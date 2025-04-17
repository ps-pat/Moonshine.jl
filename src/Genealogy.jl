using Graphs

using Bumper: UnsafeArrays
using Graphs: AbstractSimpleGraph

import Graphs: edges, vertices, ne, nv,
               eltype, edgetype, is_directed,
               has_edge, has_vertex, inneighbors, outneighbors

import GraphMakie: graphplot

using LayeredLayouts

using GeometryBasics: Point

using DataStructures: DefaultDict

import Base: IteratorSize, eltype, length

using NetworkLayout

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

sequences(genealogy, e::Edge) = (sequence(genealogy, src(e)),
                                     sequence(genealogy, dst(e)))

export mrca
"""
    mrca(genealogy[, vs = leaves(genealogy), ωs = Ω(0, ∞)])

Most recent common ancestor of a set of vertices.

See also [`tmrca`](@ref) for the time to the most recent common ancestor.
"""
function mrca end

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
    ancestral_intervals!(AIsType(), genealogy, x)

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

"""
    plot_layout(genealogy)

Layout function for genealogy plotting.
"""
function plot_layout end

plot_layout(::AbstractGenealogy) = Spring()

"""
    maxdads(Genealogy)
    maxchildren(Genealogy)

Maximum possible number of parents/children in a genealogy.

# Implementation
Must be a generated function.
"""
function maxdads end,
function maxchildren end

#############
# Utilities #
#############

"""
    idxtopos(genealogy, idx)

Return the position of the marker given its index.
"""
function idxtopos end

"""
    postoidx(genealogy, pos)

Return the largest marker's index that is at a position lesser than the one
given.
"""
function postoidx end

for (f, arg) ∈ (:idxtopos => :idx, :postoidx => :pos)
    @eval $f(genealogy::AbstractGenealogy, $arg) = $f(sam(genealogy), $arg)
end

"""
    ancestral_mask!(η, reference, x; ωs_buf = Set{Ω}(), wipe = true)
    ancestral_mask(reference, x; ωs_buf = Set{Ω}())

Mask non ancestral positions to 0. If `wipe = true`, all markers in `η` will be
initialized to 0.
"""
function ancestral_mask! end,
function ancestral_mask end

ancestral_mask!(η, genealogy::AbstractGenealogy, x; wipe = true) =
    ancestral_mask!(η, sam(genealogy), x, wipe = wipe)

ancestral_mask(genealogy::AbstractGenealogy, x) =
    ancestral_mask!(Sequence(falses(nmarkers(genealogy))), genealogy, x,
                    wipe = false)

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
                   layout = plot_layout(genealogy),
                   attributes...)
    vlabels = string.(range(1, nv(genealogy)))

    ## Color of the vertices.
    mask = fill(ancestral_mask(genealogy, ω), nv(genealogy))
    node_color = ifelse.(any.(sequences(genealogy) .& mask),
                         derived_color, wild_color)

    ## Hide non ancestral edges and vertices ##
    ewidth = DefaultDict{Edge{VertexType}, Int}(edge_width)
    for e ∈ edges(genealogy)
        isdisjoint(ancestral_intervals(genealogy, e), ω) || continue
        ewidth[e] = 0
    end

    vsize = DefaultDict{VertexType, Any}(30)
    for v ∈ ivertices(genealogy)
        isdisjoint(ancestral_intervals(genealogy, v), ω) || continue
        vsize[v] = 0
        vlabels[v] = ""
    end

    graphplot(graph(genealogy),
              layout = layout,
              node_size = vsize,
              ilabels = vlabels,
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
    nmarkers(genealogy, x)

Number of markers in the sequences of a genealogy. If an interval `x` is
specified, returns the number of markers contained in that interval.
"""
nmarkers(genealogy) = nmarkers(sam(genealogy))

nmarkers(genealogy, x) = nmarkers(sam(genealogy), x)

export branchlength
"""
    branchlength(genealogy[, e])

Total branch length of a genealogy. If an edge is specified, returns the length
of that edge.
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
function tmrca(genealogy)
    _mrca = mrca(genealogy)
    iszero(_mrca) && return zero(Float64)

    latitude(genealogy, _mrca)
end

tmrca(genealogy, vs) = latitude(genealogy, mrca(genealogy, vs))

export dads, children
export descendants!, descendants, ancestors!, ancestors
"""
    dads(genealogy, v[, ω])
    children(genealogy, v[, ω])
    ancestors!(buf, genealogy, v[, ω])
    ancestors(genealogy, v[, ω])
    descendants!(buf, genealogy, v[, ω])
    descendants(genealogy, v[, ω])

Parents/children/descendants/ancestors of a vertex. ω can be either
- a number in [0, 1] representing a position;
- an Ω representing an interval of positions;
- a set of Ωs representing multiple interval of positions.

The following rules are used to decide if an edge `e` is ancestral:
- If ω is a number, the ancestral interval of `e` must cover ω.
- If ω is an Ω or a set of Ωs, the intersection of the ancestral
  interval of `e` with ω must be non-empty.

Methods for `children` and `dads` return **references** to the underlying
adjacency lists. No touchy!
"""
function dads end, function children end,
function ancestors end, function ancestors! end,
function descendants end, function descendants! end

for (fun, list) ∈ Dict(:dads => Meta.quot(:badjlist),
                         :children => Meta.quot(:fadjlist))
    @eval $fun(genealogy, v) = getfield(graph(genealogy), $list)[v]
end

let funtransorder = Dict(:dads => (:ancestors, (x, y) -> (x, y)),
                         :children => (:descendants, (x, y) -> (y, x))),
    typesandfun = ((:Real, in), (:AI, !isdisjoint), (:(AIs{<:AbstractVector{<:AI}}), !isdisjoint))
    for (fun, (transfun, order)) ∈ funtransorder
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

        for (Argtype, testfun) ∈ typesandfun
            ## Parents & children
            @eval function $fun(genealogy, v, ω::$Argtype)
                idx = 0x06
                neig = $fun(genealogy, v)
                @inbounds @simd for k ∈ eachindex(neig)
                    ωs = ancestral_intervals(genealogy, Edge($order(neig[k], v)))
                    idx ⊻= (0x03 << 2(k - 1)) * $testfun(ω, ωs)
                end

                view(neig, range(idx & 0x03, idx >> 0x02))
            end

            ## Ancestors & descendants
            @eval function $transfun!(buf, genealogy, v, ω::$Argtype)
                writeptr = readptr = firstindex(buf)
                @inbounds for u ∈ $fun(genealogy, v, ω)
                    buf[writeptr] = u
                    writeptr += 1
                end

                @inbounds while readptr < writeptr
                    v = buf[readptr]
                    readptr += 1
                    (isleaf(genealogy, v) || isroot(genealogy, v)) && continue

                    resize!(funbuf, 2)
                    for u ∈ $fun(genealogy, v, ω)
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
        end
    end
end

export siblings, sibling
"""
    siblings!(x, genealogy, v, args...)
    siblings(genealogy, v, args...)
    sibling(genealogy, v)

Return the siblings of a vertex, that is the other vertices in the genealogy
that have the same parents.

`args` is splatted into internal calls to `dads` and `children`. It can be used
to constraint search to an interval.

If you are certain that `v` only has one sibling, you can use the `sibling`
method to avoid allocation.

Argument `x` in the allocating method can either be an `AbstractArray` or a
pointer.
"""
function siblings end,
function sibling end

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

function sibling(genealogy, v)
    for _child ∈ children(genealogy, dad(genealogy, v))
        _child == v && continue
        return _child
    end
    zero(VertexType)
end

## TODO: clean and generalize
export cousins
"""
    cousins(buf, genealogy, v, args...)

Return the cousins of a vertex, that is vertices sharing the same grandparents.

`args` is splatted into internal calls to `dads` and `children`. It can be used
to constraint search to an interval.

The maximum number of cousins is ``p^2 \\times c \\times (c - 1)`` where ``p``
and ``c`` are the maximum number of parents and children respectively.

See also [`iscousin`](@ref)
"""
function cousins(buf, genealogy::T, v, args...) where T
    head_buf = firstindex(buf)
    dads_v = zeros(MVector{maxdads(T), VertexType})

    @inbounds for (i, d) ∈ (enumerate ∘ dads)(genealogy, v, args...)
        dads_v[i] = d

        for gd ∈ dads(genealogy, d, args...)
            for u ∈ children(genealogy, gd, args...)
                ## Ensure that were dealing with a true uncle rather than a
                ## previously visited dad
                flag = zero(Int8)
                @simd ivdep for dad_v ∈ dads_v
                    ## Race condition is not an issue here. The ramainder of the
                    ## loop is skiped as soon as flag != 0.
                    flag += u == dad_v
                end

                if iszero(flag)
                    for c ∈ children(genealogy, u, args...)
                        c == v && continue

                        buf[head_buf] = c
                        head_buf += 1
                    end
                end
            end
        end
    end

    resize!(buf, head_buf - 1)
end

"""
    iscousin(genealogy, u, v, args...)

Determines if two vertices are cousins.

`args` is splatted into internal calls to `dads` and `children`. It can be used
to constraint search to an interval.

See also: [`cousins`](@ref)
"""
function iscousin(genealogy, u, v, args...)
    @inbounds for d ∈ dads(genealogy, u, args...)
        for gd ∈ dads(genealogy, d, args...)
            for uu ∈ children(genealogy, gd, args...)
                flag = zero(Int8)
                @simd ivdep for c ∈ children(genealogy, uu, args...)
                    ## Race condition is not an issue here. We return `true`
                    ## if flag != 0.
                    flag += c == v
                end

                iszero(flag) || return true
            end
        end
    end

    false
end

export dad, child
"""
    dad(genealogy, v)
    child(genealogy, v)

Return the parent/child of a vertex or 0 if none. It only makes sense to
use this method if you know `v` has a single parent/child.
"""
function dad end,
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
    nmutations(genealogy[, e])

Number of mutation on a genealogy. If an edge is specified, return only the
number of mutations on that edge.
"""
function nmutations end

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

export EdgesInterval
struct EdgesInterval{T, I}
    genealogy::T
    ωs::I
    buffer::CheapStack{Edge{VertexType}}
    visited::UnsafeArray{Bool, 1}
    min_latitude::Float64
end

function EdgesInterval(genealogy, ωs, store::AbstractArray, visited,
                       root = mrca(genealogy), min_latitude = zero(Float64))
    eibuffer = CheapStack(store)
    fill!(visited, false)

    for d ∈ children(genealogy, root, ωs)
        push!(eibuffer, Edge(root => d))
    end

    EdgesInterval(genealogy, ωs, eibuffer, visited, convert(Float64, min_latitude))
end

IteratorSize(::T) where T<:EdgesInterval = Base.SizeUnknown()
IteratorSize(::Type{<:EdgesInterval}) = Base.SizeUnknown()

eltype(::EdgesInterval) = Edge{VertexType}

function iterate(iter::EdgesInterval, state = 1)
    buffer = iter.buffer
    isempty(buffer) && return nothing

    genealogy = iter.genealogy
    ωs = iter.ωs
    visited = iter.visited
    min_latitude = iter.min_latitude
    n = nleaves(genealogy)

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
            push!(buffer, Edge(s => d))
        end
    end

    e, state + 1
end

export edges_interval

edges_interval(genealogy, ωs, store, visited,
               root = mrca(genealogy), min_latitude = zero(Float64)) =
    EdgesInterval(genealogy, ωs, store, visited, root, min_latitude)

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
    edgesmap(genealogy)

Return a `Dict` that maps every edge of a genealogy to an integer in
1:ne(genealogy).
"""
edgesmap(genealogy) = Dict(reverse.(enumerate(edges(genealogy))))

export nlive, nlive!
"""
    nlive([predicate,] genealogy, lat)
    nlive([predicate,] genealogy, lat, ωs; buffer = default_buffer())
    nlive!([predicate,] counts, genealogy, lats, ωs; buffer = default_buffer())

Number of live edges at a given latitude. Search can be restricted to an
interval. If a pedicate is passed as first argument, only edged satisfying it
will be counted.
"""
function nlive end,
function nlive! end

for (signature, test) ∈ Dict(
    :(nlive(genealogy, lat)) => :true,
    :(nlive(predicate, genealogy, lat)) => :(predicate(e)))
    @eval $signature = begin
        ## The grand MRCA is live forever
        lat > tmrca(genealogy) && return one(Int)

        live = zero(Int)

        @inbounds for e ∈ edges(genealogy)
            latitude(genealogy, dst(e)) <= lat <= latitude(genealogy, src(e)) || continue
            $test || continue
            live += 1
        end

        live
    end
end

for (signature, test) ∈ Dict(
    :(nlive(genealogy, lat::Real, ωs; buffer = default_buffer())) => :true,
    :(nlive(predicate, genealogy, lat::Real, ωs; buffer = default_buffer())) => :(predicate(e)))
    @eval $(signature) = begin
        ## The grand MRCA is live forever
        lat > tmrca(genealogy) && return one(Int)

        live = zero(Int)

        @no_escape buffer begin
            store = @alloc(Edge{VertexType}, nleaves(genealogy) + nrecombinations(genealogy))
            visited = @alloc(Bool, nrecombinations(genealogy))
            for e ∈ edges_interval(genealogy, ωs, store, visited, mrca(genealogy), lat)
                latitude(genealogy, dst(e)) <= lat <= latitude(genealogy, src(e)) || continue
                $test || continue
                live += 1
            end
        end

        live
    end
end

for (signature, test) ∈ Dict(
    :(nlive!(counts, genealogy, lats::AbstractVector{<:Real}, ωs; buffer = default_buffer())) =>
    :true,
    :(nlive!(predicate, counts, genealogy, lats::AbstractVector{<:Real}, ωs; buffer = default_buffer())) =>
    :(predicate(e)))

    @eval $(signature) = begin
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
            store = @alloc(Edge{VertexType}, nleaves(genealogy) + nrecombinations(genealogy))
            visited = @alloc(Bool, nrecombinations(genealogy))
            @inbounds for e ∈ edges_interval(genealogy, ωs, store, visited, mrca(genealogy), first(lats))
                @simd ivdep for k ∈ 1:lastlat
                    counts[k] +=
                    latitude(genealogy, dst(e)) <= lats[k] <= latitude(genealogy, src(e)) && $test
                end
            end
        end

        counts
    end
end

"""
    ismutation_edge(arg, e, idx)

Determines if an edge is a mutation edge for a given marker.
"""
function ismutation_edge(genealogy, e, idx)
    _chunkidx = chunkidx(Sequence, idx)
    hs = sequences(genealogy, e)
    chunk1, chunk2 = hs[1].data.chunks[_chunkidx], hs[2].data.chunks[_chunkidx]

    mask = one(UInt64) << (idxinchunk(Sequence, idx) - 1)
    !iszero((chunk1 ⊻ chunk2) & mask)
end
