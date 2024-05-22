using Graphs

import Graphs: add_vertices!, add_edge!, rem_edge!

using Random

using StatsBase: samplepair

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

export rec_rate
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
    get(() -> Set{Ω}((Ω(0, ∞),)), arg.core.ancestral_intervals, e)

ancestral_intervals!(ωs, arg::Arg, s::VertexType, d; wipe = true) =
    ancestral_intervals!(ωs, arg, Edge(s, d); wipe = wipe)

ancestral_intervals(arg::Arg, s::VertexType, d) =
    ancestral_intervals(arg, Edge(s, d))

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

function mutationsidx!(res, ptrs, mask, ωs, arg, e,
                       firstchunk, firstidx, lastchunk, lastidx)
    η1, η2 = sequences(arg, e)
    ancestral_mask!(mask, ωs, arg, e)
    marker_mask = one(UInt64) << (firstidx - 1)
    idx = one(Int)

    @inbounds for k ∈ range(firstchunk, lastchunk - 1)
        xored_chunks = (η1.data.chunks[k] ⊻ η2.data.chunks[k]) &
            mask.data.chunks[k]

        while !iszero(marker_mask)
            if !iszero(xored_chunks & marker_mask)
                k = ptrs[idx]
                if k ≤ length(res[idx])
                    res[idx][k] = e
                else
                    push!(res[idx], e)
                end
                ptrs[idx] += 1
            end
            idx += 1
            marker_mask <<= 1
        end

        marker_mask = one(UInt64)
    end

    ## Process last chunk
    xored_chunks = (η1.data.chunks[lastchunk] ⊻ η2.data.chunks[lastchunk]) &
        mask.data.chunks[lastchunk]
    @inbounds for _ ∈ 1:lastidx
        if !iszero(xored_chunks & marker_mask)
            k = ptrs[idx]
            if k ≤ length(res[idx])
                res[idx][k] = e
            else
                push!(res[idx], e)
            end
            ptrs[idx] += 1
        end
        idx += 1
        marker_mask <<= 1
    end

    res
end

export mutation_edges!, mutation_edges
function mutation_edges!(mutations, arg, ω::Ω)
    ## Compute the chunks and indices.
    lidx = postoidx(arg, leftendpoint(ω))
    ridx = postoidx(arg, rightendpoint(ω))
    firstchunk = chunkidx(Sequence, lidx)
    firstidx = idxinchunk(Sequence, lidx)
    lastchunk = chunkidx(Sequence, ridx)
    lastidx = idxinchunk(Sequence, ridx)

    resize!(mutations, ridx - lidx + 1)
    ptrs = fill(one(Int), length(mutations))

    mask = Sequence(undef, nmarkers(arg))
    ωs = Set{Ω}()
    for edge ∈ edges_interval(arg, ω)
        mutationsidx!(mutations, ptrs, mask, ωs, arg, edge,
                      firstchunk, firstidx, lastchunk, lastidx)
    end

    @inbounds for (k, ptr) ∈ enumerate(ptrs)
        resize!(mutations[k], ptr - 1)
    end

    mutations
end

function mutation_edges(arg, ω::Ω)
    m = postoidx(arg, rightendpoint(ω)) - postoidx(arg, leftendpoint(ω)) + 1
    ret = [Vector{EdgeType}(undef, nleaves(arg) ÷ 2 + 1) for _ ∈ 1:m]
    mutation_edges!(ret, arg, ω)
end

##################
# Recombinations #
##################

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
    ωs_buf1 = Set{Ω}()
    ωs_buf2 = Set{Ω}()
    oldη = Sequence(undef, nmarkers(arg))

    while !isempty(vstack)
        v = pop!(vstack)
        _dads = dads(arg, v)

        ## Update ancestry of edges.
        @inbounds for dad ∈ _dads
            arg.core.ancestral_intervals[Edge(dad, v)] =
                ancestral_intervals!(ωs_buf2, arg, Edge(dad, v)) ∩
                ancestral_intervals!(ωs_buf1, arg, v)
        end

        ## Update sequence of v.
        oldη.data.chunks .= sequence(arg, v).data.chunks
        _compute_sequence!(arg, v, mask, ωs_buf1)

        (sequence(arg, v) == oldη || isempty(_dads)) && continue
        push!(vstack, _dads...)
    end

    arg
end

export recombine!
"""
    recombine!(arg, redge, cedge, breakpoint, rlat, clat)

Add a recombination event to an ARG.
"""
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
    ## TODO: Wrong. Should add more edges.

    ## Compute sequence of new vertices.
    let mask = Sequence(undef, nmarkers(arg)),
        ωs = Set{Ω}()

        _compute_sequence!(arg, rvertex, mask, ωs)
        _compute_sequence!(arg, cvertex, mask, ωs)
    end

    ## Update sequences and ancetral intervals.
    update_upstream!(arg, src(redge))
    isroot(arg, dst(cedge)) || update_upstream!(arg, src(cedge))

    arg
end

function build!(rng, arg::Arg; winwidth = ∞)
    # λc = rec_rate(arg, true) / 2 * seq_length(arg)

    _mutation_edges = mutation_edges(arg, Ω(0, ∞))
    nextidx = findfirst(>(1) ∘ length, _mutation_edges)
    isnothing(nextidx) && return arg ## ARG is already consistent.
    breakpoint_lbound = idxtopos(arg, max(nextidx - 1, one(nextidx)))

    while !isnothing(nextidx)
        live_edges = _mutation_edges[nextidx + length(_mutation_edges) - nmarkers(arg)]

        ## Sample next breakpoint.
        nextpos = idxtopos(arg, nextidx)
        @debug "Recombination bounds" breakpoint_lbound nextpos
        if nextpos > breakpoint_lbound + eps(Float64)
            breakpoint_dist = LogUniform(breakpoint_lbound + eps(Float64),
                                         nextpos)
            breakpoint = rand(rng, breakpoint_dist)
            arg.logprob += logpdf(breakpoint_dist, breakpoint)
        else
            breakpoint = breakpoint_lbound
        end

        ## Create ancestral window.
        ω = Ω(breakpoint - winwidth / 2, breakpoint + winwidth / 2)

        ## Sample recombination and recoalescence edges.
        e1, e2 = samplepair(rng, length(live_edges))
        redge = live_edges[e1]
        cedge = live_edges[e2]
        if latitude(arg, dst(redge)) > latitude(arg, dst(cedge))
            redge, cedge = cedge, redge
        end

        ## Sample recombination latitude.
        ## It would be more accurate to use a truncated exponential
        ## distribution with rate proportional to the number of live
        ## vertices at the latitude of dst(redge). However, this is a
        ## substantial amount of trouble for very little gain...
        rlat_bounds = (latitude(arg, dst(redge)) + eps(Float64),
                       min(latitude(arg, src(redge)), latitude(arg, src(cedge))))
        if first(rlat_bounds) < last(rlat_bounds)
            rlat_dist = Uniform(rlat_bounds...)
            rlat = rand(rng, rlat_dist)
            arg.logprob += logpdf(rlat_dist, rlat)
        else
            rlat = first(rlat_bounds)
        end

        ## Sample recoalescence latitude.
        clat_bounds = (max(rlat, latitude(arg, dst(cedge))),
                       latitude(arg, src(cedge)))
        if first(clat_bounds) < last(clat_bounds)
            clat_dist = Uniform(clat_bounds...)
            clat = rand(rng, clat_dist)
            arg.logprob += logpdf(clat_dist, clat)
        else
            clat = first(clat_bounds)
        end

        ## Add recombination event to the graph.
        @debug "Recombination event" redge, cedge, breakpoint, rlat, clat
        recombine!(arg, redge, cedge, breakpoint, rlat, clat)

        ## Update variables.
        mutation_edges!(_mutation_edges, arg, Ω(breakpoint, ∞))
        nextidx_new = findfirst(>(1) ∘ length, _mutation_edges)

        isnothing(nextidx_new) && break
        nextidx_new += nmarkers(arg) - length(_mutation_edges)
        breakpoint_lbound = nextidx_new == nextidx ? breakpoint :
            idxtopos(arg, max(nextidx_new - 1, one(nextidx_new)))
        nextidx = nextidx_new
    end
    arg
end
