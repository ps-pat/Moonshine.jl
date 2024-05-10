using Graphs

import Graphs: add_vertices!, add_edge!, rem_edge!

######################
# ArgCore Definition #
######################

struct ArgCore
    graph::SimpleDiGraph{VertexType}
    latitudes::Vector{Float64}
    sequences::Vector{Sequence}
    ancestral_intervals::Dict{EdgeType, Set{Ω}}

    positions::Vector{Float64}
    seq_length::Float64
    Ne::Float64
    μloc::Float64
    ρloc::Float64
end

ArgCore() = ArgCore(SimpleDiGraph{VertexType}(),
                    Float64[], Sequence[], Dict{EdgeType, Set{Ω}}(),
                    Float64[], typemin(Float64),
                    typemin(Float64), typemin(Float64), typemin(Float64))

ArgCore(treecore::TreeCore, ρloc = 1e-5) =
    ArgCore(treecore.graph, treecore.latitudes,
            treecore.sequences, Dict{EdgeType, Set{Ω}}(),
            treecore.positions, treecore.seq_length,
            treecore.Ne, treecore.μloc, ρloc)

function ArgCore(leaves::AbstractVector{Sequence};
                 positions = nothing,
                 seq_length = one(Float64),
                 Ne = one(Float64),
                 μloc = 1e-7,
                 ρloc = 1e-5)
    treecore = TreeCore(leaves,
                        positions = positions,
                        seq_length = seq_length,
                        Ne = Ne,
                        μloc = μloc)
    ArgCore(treecore, ρloc)
end

##################
# Arg Definition #
##################

mutable struct Arg <: AbstractGenealogy
    core::ArgCore
    logprob:: BigFloat
end
export Arg

Arg() = Arg(ArgCore(), zero(BigFloat))

Arg(tree::Tree, ρloc = 1e-5) =
    Arg(ArgCore(tree.core, ρloc), prob(tree, logscale = true))

# Arg(leaves::AbstractVector{Sequence}; genpars...) =
#     Arg(ArgCore(leaves, genpars...), zero(BigFloat))

###############################
# AbstractGenealogy Interface #
###############################

function nleaves(arg::Arg)
    ret = zero(Int)
    for η ∈ range(1, nv(arg))
        iszero(outdegree(arg, η)) || return ret
        ret += 1
    end

    ret
end

describe(::Arg, long = true) = long ? "Ancestral Recombination Graph" : "ARG"

###########################
# AbstractGraph Interface #
###########################

function add_vertices!(arg::Arg, H, lats)
    append!(latitudes(arg), lats)
    append!(sequences(arg), H)
    add_vertices!(graph(arg), length(H))
end

function add_edge!(arg::Arg, e, ints::Set{Ω})
    arg.core.ancestral_intervals[e] = ints
    add_edge!(graph(arg), e)
end

rem_edge!(arg::Arg, e) = rem_edge!(graph(arg), e)

##################
# Simple methods #
##################

rec_rate(arg::Arg, scaled = true) = arg.core.ρloc * (scaled ? 4 * Ne(arg) : 1.)

########################
# Ancestrality Methods #
########################

"""
    ancestral_intervals(arg, e)
    ancestral_intervals(arg, σ, δ)
    ancestral_intervals(arg, v)

Compute the interval for which an element is ancestral. Behavior is undefined
for non-existent element.
"""
function ancestral_intervals end
export ancestral_intervals

function ancestral_intervals!(ωs, arg::Arg, e::EdgeType; wipe = true)
    wipe && empty!(ωs)

    haskey(arg.core.ancestral_intervals, e) || return push!(ωs, Ω(0, ∞))

    for ω ∈ arg.core.ancestral_intervals[e]
        push!(ωs, ω)
    end

    ωs
end

ancestral_intervals(arg::Arg, e::EdgeType) =
    ancestral_intervals!(Set{Ω}(), arg, e)

ancestral_intervals!(ωs, arg::Arg, s::VertexType, d; wipe = true) =
    ancestral_intervals!(ωs, arg, Edge(s, d); wipe = wipe)

ancestral_intervals(arg::Arg, s::VertexType, d) =
    ancestral_intervals!(Set{Ω}(), arg, s, d)

function ancestral_intervals!(ωs, arg::Arg, v::VertexType; wipe = true)
    empty!(ωs)
    isleaf(arg, v) && return push!(ωs, Ω(0, ∞))

    for child ∈ children(arg, v)
        ancestral_intervals!(ωs, arg, v, child, wipe = false)
    end

    simplify!(ωs)
end

ancestral_intervals(arg::Arg, v::VertexType) =
    ancestral_intervals!(Set{Ω}(), arg, v)

#######
# MMN #
#######

function mutationsidx!(res, mask, ωs, arg, e, firstchunk, firstidx, lastchunk, lastidx)
    η1, η2 = sequences(arg, e)
    ancestral_mask!(mask, ωs, arg, e)
    marker_mask = one(UInt64) << (firstidx - 1)
    idx = one(Int)

    @inbounds for k ∈ range(firstchunk, lastchunk - 1)
        xored_chunks = (η1.data.chunks[k] ⊻ η2.data.chunks[k]) &
            mask.data.chunks[k]

        while !iszero(marker_mask)
            iszero(xored_chunks & marker_mask) || push!(res[idx], e)
            idx += 1
            marker_mask <<= 1
        end

        marker_mask = one(UInt64)
    end

    ## Process last chunk
    xored_chunks = (η1.data.chunks[lastchunk] ⊻ η2.data.chunks[lastchunk]) &
        mask.data.chunks[lastchunk]
    for _ ∈ 1:lastidx
        iszero(xored_chunks & marker_mask) || push!(res[idx], e)
        idx += 1
        marker_mask <<= 1
    end

    res
end

export mutation_edges
function mutation_edges(arg, ω::Ω)
    ## Compute the chunks and indices.
    lidx = postoidx(arg, leftendpoint(ω))
    ridx = postoidx(arg, rightendpoint(ω))
    firstchunk = chunkidx(Sequence, lidx)
    firstidx = idxinchunk(Sequence, lidx)
    lastchunk = chunkidx(Sequence, ridx)
    lastidx = idxinchunk(Sequence, ridx)

    m = ridx - lidx + 1
    mutations = [Vector{EdgeType}() for _ ∈ 1:m]

    mask = Sequence(undef, nmarkers(arg))
    ωs = Set{Ω}()
    for edge ∈ edges_interval(arg, ω)
        mutationsidx!(mutations, mask, ωs, arg, edge,
                      firstchunk, firstidx, lastchunk, lastidx)
    end

    mutations
end

##################
# Recombinations #
##################

export recombine!
function recombine!(arg, redge, cedge, breakpoint, rlat, clat)
    ## Add recombination and recoalescence vertices to arg.
    rvertex, cvertex = nv(arg) .+ (1, 2)
    add_vertices!(arg,
                  (Sequence(trues(nmarkers(arg))),
                   Sequence(trues(nmarkers(arg)))),
                  (rlat, clat))

    ## Replace recombination edge.
    ωr = ancestral_intervals(arg, redge)
    add_edge!(arg, Edge(rvertex, dst(redge)), ωr)
    add_edge!(arg, Edge(src(redge), rvertex), ωr ∩ Ω(0, breakpoint))
    rem_edge!(arg, redge)

    ## Replace recoalescence edge.
    ωc = ancestral_intervals(arg, cedge)
    add_edge!(arg, Edge(cvertex, rvertex), ωr ∩ Ω(breakpoint, ∞))
    add_edge!(arg, Edge(cvertex, dst(cedge)), ωc)
    rem_edge!(arg, cedge)
    isroot(arg, dst(cedge)) ||
        add_edge!(arg, Edge(src(cedge), cvertex), ωc ∪ (ωr ∩ Ω(breakpoint, ∞)))

    ## Compute sequence of new vertices.
    let mask = Sequence(undef, nmarkers(arg)),
        ωs = Set{Ω}()

        _compute_sequence!(arg, rvertex, mask, ωs)
        _compute_sequence!(arg, cvertex, mask, ωs)
    end

    ## Update sequences and ancetral intervals.
    update_upstream!(arg, src(redge))
    isroot(arg, dst(cedge)) || update_upstream!(arg, cvertex)

    arg
end

function _compute_sequence!(arg, v, mask, ωs)
    η = sequence(arg, v)
    η.data.chunks .⊻= .~η.data.chunks

    @inbounds for child ∈ children(arg, v)
        ancestral_mask!(mask, arg, ancestral_intervals!(ωs, arg, Edge(v, child)))
        η.data.chunks .&= sequence(arg, child).data.chunks .| .~mask.data.chunks
    end

    ancestral_mask!(mask, ωs, arg, v)
    η.data.chunks .&= mask.data.chunks
    η
end

function update_upstream!(arg, v)
    vstack = Stack{VertexType}(ceil(Int, log(nv(arg))))
    push!(vstack, v)

    mask = Sequence(undef, nmarkers(arg))
    ωs = Set{Ω}()
    oldωs = Set{Ω}()
    oldη = Sequence(undef, nmarkers(arg))

    while !isempty(vstack)
        v = pop!(vstack)

        ## Update ancestry of the parents of v
        for dad ∈ dads(arg, v)
            empty!(oldωs)
            ωdad = ancestral_intervals(arg, Edge(dad, v))
            while !isempty(ωdad)
                push!(oldωs, pop!(ωdad))
            end

            ancestral_intervals!(ωdad, arg, v)
            ωdad == oldωs && continue
            push!(vstack, dad)
        end

        ## Update sequence of v.
        oldη.data.chunks .= sequence(arg, v).data.chunks
        _compute_sequence!(arg, v, mask, ωs)

        sequence(arg, v) == oldη && continue
        push!(vstack, dads(arg, v)...)
    end

    arg
end
