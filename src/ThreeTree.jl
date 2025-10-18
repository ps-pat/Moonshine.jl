import Graphs:
    AbstractSimpleGraph,
    is_directed,
    inneighbors, outneighbors,
    nv, has_vertex,
    SimpleEdge, edges, edgetype, has_edge, ne, rem_edge!

"""
    $(TYPEDEF)

Simple data structure that can be used to record the topology of an ancestral
recombination graph.

The first `nleaves` entries of `neig` contain the parent vertices of the leaves.
The following `3nleaves - 3` store information about coalescence vertices.
The remaining entries store data about recombination and coalescence vertices
in an alternating fashion.

The first two entries of a triple corresponding to a coalescence/recoalescence
vertex contain its downstream vertices, and the remaining entry contains its
parent. The order is reversed for recombination vertices.

Types [`Tree`](@ref) and [`Arg`](@ref) use this data structure.

# Fields
$(TYPEDFIELDS)

# Constructors

$(METHODLIST)

--*Internal*--
"""
struct ThreeTree{T<:Integer} <: AbstractSimpleGraph{T}
    "Number of leaves (vertices with outdegree 0)"
    nleaves::T
    "Neighborhood data"
    neig::Vector{T}
end

ThreeTree(nleaves) = ThreeTree(nleaves, zeros(typeof(nleaves), nleaves))

## Helper Functions

_is_tt_leaf(tt, v) = v <= tt.nleaves

function _is_tt_coalescence(tt, v)
    n = tt.nleaves
    (v > 2n && isodd(v)) || n < v < 2n
end

_is_tt_recombination(tt, v) = v >= 2tt.nleaves && iseven(v)

@inbounds function _add_coalescence_vertex_idx_helper!(tt, child, new_vertex)
    n, neig = tt.nleaves, tt.neig

    if _is_tt_leaf(tt, child)
        iszero(neig[child]) && return child
    else
        i = 3child - 2n

        if _is_tt_coalescence(tt, child)
            iszero(neig[i]) && return i
        else
            ## Recombination vertex
            i -= 2
            j = i + 1

            if iszero(neig[i])
                if neig[j] > new_vertex
                    return i
                else
                    neig[i] = neig[j]
                    return j
                end
            elseif iszero(neig[j])
                if neig[i] < new_vertex
                    return j
                else
                    neig[j] = neig[i]
                    return i
                end
            end
        end
    end

    zero(new_vertex)
end

@inbounds function add_coalescence_vertex!(tt::ThreeTree, child1, child2)
    n, neig = tt.nleaves, tt.neig

    child1, child2 = minmax(child1, child2)
    new_vertex_idx = length(neig) + 1
    new_vertex =  n + div(new_vertex_idx - n, 3, RoundUp)

    idx1 = _add_coalescence_vertex_idx_helper!(tt, child1, new_vertex)
    iszero(idx1) && return false
    idx2 = _add_coalescence_vertex_idx_helper!(tt, child2, new_vertex)
    iszero(idx2) && return false

    push!(neig, child1, child2, 0)
    neig[idx1] = new_vertex
    neig[idx2] = new_vertex

    true
end

function _compute_idx(tt, u, v, predicate)
    n, neig = tt.nleaves, tt.neig

    uidx = _is_tt_leaf(tt, u) ? u : 3u - 2n
    if predicate(tt, u)
        uidx -= 2
        if neig[uidx] != v
            uidx += 1
        end
    end

    uidx
end

@inbounds function add_recombination_vertex!(tt::ThreeTree, redge)
    n, neig = tt.nleaves, tt.neig

    rvertex_idx = length(neig) + 1
    rvertex =  n + div(rvertex_idx - n, 3, RoundUp)

    s, d = src(redge), dst(redge)
    sidx = _compute_idx(tt, s, d, _is_tt_coalescence)
    didx = _compute_idx(tt, d, s, _is_tt_recombination)

    push!(neig, s, 0, d)
    neig[sidx] = rvertex
    neig[didx] = rvertex

    true
end

@inbounds function add_recoalescence_vertex!(tt::ThreeTree, cedge, rvertex)
    n, neig = tt.nleaves, tt.neig

    cvertex_idx = length(neig) + 1
    cvertex =  n + div(cvertex_idx - n, 3, RoundUp)

    s, d = src(cedge), dst(cedge)
    sidx = s == d ?
        (zero ∘ eltype)(tt) : _compute_idx(tt, s, d, _is_tt_coalescence)
    didx = _compute_idx(tt, d, s, _is_tt_recombination)

    push!(neig, d, rvertex)
    neig[3rvertex - 2n - 1] = cvertex
    neig[didx] = cvertex
    if iszero(sidx)
        push!(neig, 0)
    else
        push!(neig, s)
        neig[sidx] = cvertex
    end

    true
end

function add_rr_event!(tt::ThreeTree, redge, cedge)
    add_recombination_vertex!(tt, redge)
    rvertex = nv(tt)
    add_recoalescence_vertex!(tt, cedge, rvertex)

    true
end

# AbstractGraphs Interface

is_directed(::Type{<:ThreeTree}) = true
is_directed(::ThreeTree) = true

function inneighbors(tt::ThreeTree, v::Integer)
    n, neig = tt.nleaves, tt.neig

    ptr = pointer(neig, _is_tt_leaf(tt, v) ? v : 3v - 2n)
    len = 1
    if _is_tt_recombination(tt, v)
        ptr -= 2(sizeof ∘ eltype)(neig)
        len += 1
    elseif (iszero ∘ unsafe_load)(ptr)
        return UnsafeArray(convert(Ptr{eltype(neig)}, C_NULL), (0,))
    end

    UnsafeArray(ptr, (len,))
end

function outneighbors(tt::ThreeTree, v::Integer)
    n, neig = tt.nleaves, tt.neig

    _is_tt_leaf(tt, v) &&
        return UnsafeArray(convert(Ptr{eltype(neig)}, C_NULL), (0,))

    ptr = pointer(neig, 3v - 2n)
    len = 1
    if _is_tt_coalescence(tt, v)
        ptr -= 2(sizeof ∘ eltype)(neig)
        len += 1
    end

    UnsafeArray(ptr, (len,))
end

function nv(tt::ThreeTree)
    n, neig = tt.nleaves, tt.neig
    ret = n + ((length(neig)) - n) ÷ 3
    convert(VertexType, ret)
end

has_vertex(tt::ThreeTree, v) = 1 <= v <= nv(tt)

ne(tt::ThreeTree) = (3nv(tt) - 1) ÷ 2 - tt.nleaves

edgetype(::ThreeTree) = SimpleEdge{VertexType}

function edges(tt::ThreeTree)
    n = tt.nleaves

    build_edges = function(s)
        Iterators.map(d -> SimpleEdge(s => d), outneighbors(tt, s))
    end

    ## Coalescence vertices
    coal_range = range(n + one(n), length = VertexType(n - 1))
    recoal_range = range(VertexType(2n + 1), nv(tt), step = VertexType(2))
    coal_it = Iterators.flatten((coal_range, recoal_range))
    coal = Iterators.flatmap(build_edges, coal_it)

    ## Recombination vertices
    rec_it = range(VertexType(2n), nv(tt), step = VertexType(2))
    rec = Iterators.flatmap(build_edges, rec_it)

    Iterators.filter(e -> !(iszero ∘ src)(e) && (!iszero ∘ dst)(e),
        Iterators.flatten((coal, rec)))
end

function has_edge(tt::ThreeTree, e)
    s, d = src(e), dst(e)
    d ∈ outneighbors(tt, s)
end
