import Graphs:
    AbstractSimpleGraph,
    is_directed,
    inneighbors, outneighbors,
    nv, has_vertex,
    SimpleEdge, edges, edgetype, has_edge, ne, rem_edge!

struct ThreeTree{T<:Integer} <: AbstractSimpleGraph{T}
    nleaves::T
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
    end

    ## Root
    (iszero ∘ unsafe_load)(ptr) &&
        return UnsafeArray(convert(Ptr{eltype(neig)}, C_NULL), (0,))

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
    recoal_range = range(VertexType(2n + 1), nv(tt), step = 2)
    coal_it = Iterators.flatten((coal_range, recoal_range))
    coal = Iterators.flatmap(build_edges, coal_it)

    ## Recombination vertices
    rec_it = range(VertexType(2n), nv(tt), step = 2)
    rec = Iterators.flatmap(build_edges, rec_it)

    Iterators.filter(e -> !(iszero ∘ src)(e) && (!iszero ∘ dst)(e),
        Iterators.flatten((coal, rec)))
end

function has_edge(tt::ThreeTree, e)
    s, d = src(e), dst(e)
    d ∈ outneighbors(tt, s)
end
