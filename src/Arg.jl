using Graphs

import Graphs: add_vertices!, add_edge!, rem_edge!

using Random

using StatsBase: samplepair, ProbabilityWeights

using SparseArrays

using Combinatorics: combinations

using Distributions

using StaticArrays: @SVector

##################
# Arg Definition #
##################

export Arg
struct Arg <: AbstractGenealogy
    graph::SimpleDiGraph{VertexType}
    latitudes::Vector{Float64}
    recombination_mask::Vector{AIsType}
    mrca::Base.RefValue{VertexType}
    sequences::Vector{Sequence}
    ancestral_intervals::Dict{Edge{VertexType}, AIsType}
    sample::Sample
    logprob::Base.RefValue{Float64x2}
end

function Arg(tree::Tree)
    ancestral_intervals = Dict{Edge{VertexType}, AIsType}()
    for e ∈ edges(tree)
        ancestral_intervals[e] = AIsType([Ω(0, ∞)])
    end

    Arg(graph(tree),
        latitudes(tree),
        Vector{AIsType}(undef, 0),
        Ref(mrca(tree)),
        sequences(tree),
        ancestral_intervals,
        sam(tree),
        Ref(prob(tree, logscale = true)))
end

function Arg(rng::AbstractRNG, n, μ, ρ, Ne, sequence_length)
    tree = Tree(rng, n, μ, ρ, Ne, sequence_length)
    build!(rng, tree)
    Arg(tree)
end

###############################
# AbstractGenealogy Interface #
###############################

nleaves(arg::Arg) = length(sam(arg).H)

describe(::Arg, long = true) = long ? "Ancestral Recombination Graph" : "ARG"

isrecombination(::Arg, v, n) = iseven(v) && v >= 2n

isrecombination(arg::Arg, v) = isrecombination(arg, v, nleaves(arg))

function recombinations(arg::Arg)
    isempty(arg.recombination_mask) && return StepRange{Int, Int}(0, 1, 0)

    start = 2nleaves(arg)
    step = 2
    stop = nv(arg)

    StepRange{Int, Int}(start, step, stop)
end

nrecombinations(arg::Arg) = ne(arg) - nv(arg) + 1

mrca(arg::Arg) = arg.mrca[]

mrca(arg, ωs) = mrca(arg, leaves(arg), ωs)

@generated maxdads(::Type{Arg}) = 2

@generated maxchildren(::Type{Arg}) = 2

###########################
# AbstractGraph Interface #
###########################

function add_vertices!(arg::Arg, H, lats)
    append!(latitudes(arg), lats)
    append!(sequences(arg), H)
    add_vertices!(graph(arg), length(H))
end

function add_edge!(arg::Arg, e, ints::AIsType)
    arg.ancestral_intervals[e] = ints
    add_edge!(graph(arg), e)
end

rem_edge!(arg::Arg, e) = rem_edge!(graph(arg), e)

function plot_layout(arg::Arg)
    initxs = rand(1:nleaves(arg), nv(arg))
    ys = vcat(zeros(nleaves(arg)), latitudes(arg))
    Spring(
        initialpos = (collect ∘ zip)(initxs, ys),
        # initialpos = vcat(
            # [(v, 0) for v ∈ leaves(arg)],
            # [(1, 10) for _ ∈ ivertices(arg)]
        # ),
        pin = [(false, true) for _ ∈ vertices(arg)])
end

########################
# Ancestrality Methods #
########################

function ancestral_intervals!(ωs, arg::Arg, e::Edge; wipe = true, simplify = true)
    wipe && empty!(ωs)

    union!(ωs, ancestral_intervals(arg, e), simplify = simplify)

    ωs
end

ancestral_intervals(arg::Arg, e::Edge) =
    get!(() -> AIsType([Ω(0, ∞)]), arg.ancestral_intervals, e)

function ancestral_intervals!(ωs, arg::Arg, v::VertexType; wipe = true)
    wipe && empty!(ωs)
    isleaf(arg, v) && return push!(ωs, Ω(0, ∞))

    for child ∈ children(arg, v)
        ancestral_intervals!(ωs, arg, Edge(v => child), wipe = false, simplify = false)
    end

    simplify!(ωs)
end

ancestral_intervals(arg::Arg, v::VertexType) = ancestral_intervals!(AIsType(), arg, v)

ancestral_mask!(η, arg::Arg, e::Edge{VertexType}; wipe = true) =
    ancestral_mask!(η, sam(arg), ancestral_intervals(arg, e), wipe = wipe)

function ancestral_mask!(h, arg::Arg, v::VertexType;
                         buffer = default_buffer(), wipe = true)
    ## Compute number of Ωs
    len = sum(c -> (length ∘ ancestral_intervals)(arg, Edge(v => c)),
              children(arg, v))

    @no_escape buffer begin
        ωs = @alloc(Ω, len)
        k = 1
        for c ∈ children(arg, v)
            for ω ∈ ancestral_intervals(arg, Edge(v, c))
                ωs[k] = ω
                k += 1
            end
        end
        ancestral_mask!(h, sam(arg), ωs, wipe = wipe)
        nothing # Workaround, see Bumper.jl's issue #49
    end

    h
end

ancestral_mask(arg::Arg, x::Union{VertexType, Edge{VertexType}}) =
    ancestral_mask!(Sequence(falses(nmarkers(arg))), arg, x)

recidx(arg, v) = (v - 2(nleaves(arg) - 1)) ÷ 2

function ancestral_mask(e::Edge, arg)
    s, d = src(e), dst(e)

    inc = s > otherdad(arg, s, d)
    arg.recombination_mask[2recidx(arg, d) - 1 + inc]
end

#          +----------------------------------------------------------+
#          |                      Other methods                       |
#          +----------------------------------------------------------+

"""
    otherdad(arg, s, d)
    otherdad(arg, e)

Return the parent of `d` that is not `s` for a recombination vertex `d`. If `d`
is not a recombination vertex, returns `s`. Can also take an edge as argument.
"""
function otherdad(arg, s, d)
    ret = s
    for dad ∈ dads(arg, d)
        ret == dad && continue
        ret = dad
        break
    end

    ret
end

otherdad(arg, e) = otherdad(arg, src(e), dst(e))

#######
# MMN #
#######

function mutationsidx!(res, mask, arg, e, firstchunk, firstidx, lastchunk;
                       ωs_buf = AIsType())
    η1, η2 = sequences(arg, e)
    ancestral_mask!(mask, arg, e)
    marker_mask = one(UInt64) << (firstidx - 1)
    idx = one(Int)

    @inbounds for k ∈ range(firstchunk, lastchunk)
        xored_chunk = (η1.data.chunks[k] ⊻ η2.data.chunks[k]) & mask.data.chunks[k]

        while !iszero(marker_mask)
            iszero(xored_chunk & marker_mask) || push!(res[idx], e)
            idx += 1
            marker_mask <<= 1
        end

        marker_mask = one(UInt64)
    end

    res
end

export mutation_edges!, mutation_edges
function mutation_edges!(mutations, arg, ω::Ω; buffer = default_buffer())
    ## Compute the chunks and indices.
    idx = postoidx(arg, ω)
    lidx, ridx = first(idx), last(idx)

    firstchunk = chunkidx(Sequence, lidx)
    firstidx = idxinchunk(Sequence, lidx)
    lastchunk = chunkidx(Sequence, ridx)

    @inbounds for k ∈ eachindex(mutations)
        resize!(mutations[k], 0)
    end

    mask = Sequence(undef, nmarkers(arg))
    ωs_buf = AIsType()
    @no_escape buffer begin
        store = @alloc(Edge{VertexType}, nleaves(arg) + nrecombinations(arg))
        visited = @alloc(Bool, nrecombinations(arg))
        @inbounds for e ∈ edges_interval(arg, ω, store, visited)
            mutationsidx!(mutations, mask, arg, e, firstchunk, firstidx, lastchunk,
                          ωs_buf = ωs_buf)
        end
    end

    mutations
end

function mutation_edges!(mutations, arg, idx::Int; buffer = default_buffer())
    pos = idxtopos(arg, idx)
    empty!(mutations)

    @no_escape buffer begin
        store = @alloc(VertexType, nleaves(arg))
        vstack = CheapStack(store)

        push!(vstack, mrca(arg))
        while !isempty(vstack)
            s = pop!(vstack)
            for d ∈ children(arg, s, pos)
                if sequence(arg, d)[idx]
                    push!(mutations, Edge(s => d))
                    continue
                end

                push!(vstack, d)
            end
        end
    end

    mutations
end

mutation_edges(arg, idx::Int; buffer = default_buffer()) =
    mutation_edges!(Vector{Edge{VertexType}}(undef, 0), arg, idx, buffer = buffer)

function mutation_edges(arg)
    ret = [Vector{Edge{VertexType}}(undef, 0) for _ ∈ 1:nmarkers(arg)]
    mutation_edges!(ret, arg, Ω(0, ∞))
end

function _compute_mutations_sequences(mutations_sequences, c_buffer, mask)
    @turbo for k ∈ 1:length(mutations_sequences)
        mutations_sequences[k] ⊻= c_buffer[k]
        mutations_sequences[k] &= mask
    end

    mutations_sequences
end

function next_inconsistent_idx(arg, idx;
                               mutations_edges = SVector{64}([Edge{VertexType}[] for _ ∈ 1:64]),
                               buffer = default_buffer())
    ## The mask restricts search to markers in (original) `idx` and
    ## `nmarkers(arg)` inclusively.
    mask = typemax(UInt64)
    mask <<= idxinchunk(Sequence, idx) - 1

    @inbounds while idx <= nmarkers(arg)
        empty!.(mutations_edges)

        ωlbound = idxtopos(arg, idx)

        idx_chunk = chunkidx(Sequence, idx)
        idx = 64idx_chunk + 1 # idx is now the first marker of the next chunk
        if idx > nmarkers(arg)
            mask >>>= 64 - idxinchunk(Sequence, nmarkers(arg))
            ωubound = ∞
        else
            ωubound = idxtopos(arg, idx)
        end

        base_ω = Ω(ωlbound, ωubound)

        @no_escape buffer begin
            store = @alloc(Edge{VertexType}, nleaves(arg) + nrecombinations(arg))
            visited = @alloc(Bool, nrecombinations(arg))
            ei_ptr = convert(Ptr{Edge{VertexType}}, @alloc_ptr(ne(arg) * sizeof(Edge{VertexType})))
            mutations_sequences_ptr = convert(Ptr{UInt}, @alloc_ptr((ne(arg)) * sizeof(UInt)))
            c_buffer_ptr = convert(Ptr{UInt}, @alloc_ptr((ne(arg)) * sizeof(UInt)))

            ## Traverse marginal graph and fill containers
            ne_interval = 0
            for e ∈ edges_interval(arg, base_ω, store, visited)
                ne_interval += 1

                unsafe_store!(ei_ptr, e, ne_interval)

                c1 = sequence(arg, src(e)).data.chunks[idx_chunk]
                c2 = sequence(arg, dst(e)).data.chunks[idx_chunk]
                unsafe_store!(mutations_sequences_ptr, c1, ne_interval)
                unsafe_store!(c_buffer_ptr, c2, ne_interval)
            end

            ## Wrap pointers into arrays
            mutations_sequences =
                UnsafeArray{UInt, 1}(mutations_sequences_ptr, (ne_interval,))
            c_buffer =
                UnsafeArray{UInt, 1}(c_buffer_ptr, (ne_interval,))

            ## Compute mutations sequences
            _compute_mutations_sequences(mutations_sequences, c_buffer, mask)

            ## Find mutation edges
            ei = UnsafeArray{Edge{VertexType}, 1}(ei_ptr, (ne_interval,))
            for (i, e) ∈ enumerate(ei)
                mutations_sequence = mutations_sequences[i]

                acc = 0
                while true
                    j = trailing_zeros(mutations_sequence) + 1
                    j > 64 && break
                    pos = idxtopos(arg, 64(idx_chunk - 1) + acc + j)
                    if pos ∈ ancestral_intervals(arg, e) && pos ∈ base_ω
                        push!(mutations_edges[acc + j], e)
                    end
                    mutations_sequence >>>= j
                    acc += j
                end
            end
        end

        idx_mutation_chunk = findfirst(>(1) ∘ length, mutations_edges)
        if !isnothing(idx_mutation_chunk)
            mutation_idx = 64(idx_chunk - 1) + idx_mutation_chunk
            mutation_edges = mutations_edges[idx_mutation_chunk]
            return mutation_idx, mutation_edges
        end

        empty!.(mutations_edges)
        mask = typemax(UInt64)
    end

    0, Edge{VertexType}[]
end

##################
# Recombinations #
##################
function _compute_sequence!(arg, v, mask)
    ## This function assumes that every marker is set to 1!
    η = sequence(arg, v)

    @inbounds for child ∈ children(arg, v)
        ancestral_mask!(mask, arg, Edge(v, child))
        η.data.chunks .&= sequence(arg, child).data.chunks .| .~mask
    end

    ancestral_mask!(mask, arg, v)
    η.data.chunks .&= mask
    η
end

function _update_ai!(vstack, arg, e, ωs, oldhash)
    oldhash = hash(oldhash, (hash ∘ ancestral_intervals)(arg, e))

    ## Update ancestral interval of e
    copy!(ancestral_intervals(arg, e), ωs)
    empty!(ωs)

    ## Update vstack
    newhash = hash((hash ∘ sequence)(arg, dst(e)),
                   (hash ∘ ancestral_intervals)(arg ,e))
    if (oldhash) != newhash
        push!(vstack, src(e))
    end

    ancestral_intervals(arg, e)
end

function update_upstream!(arg, v; buffer = default_buffer())
    @no_escape buffer begin
        store = @alloc(VertexType, nrecombinations(arg) + 1)
        vstack = CheapStack(store)
        push!(vstack, v)

        mask = @alloc(UInt, div(nmarkers(arg), blocksize(Sequence), RoundUp))
        ωsv = (AIsType(), AIsType())

        while !isempty(vstack)
            v = pop!(vstack)

            ## Update sequence of `v` ##
            h = sequence(arg, v)
            oldhash = hash(h)
            h.data.chunks .⊻= .~h.data.chunks
            _compute_sequence!(arg, v, mask)
            iszero(dad(arg, v)) && continue

            ## Update ancestral intervals of parental edges ##
            ancestral_intervals!(first(ωsv), arg, v)

            if isrecombination(arg, v)
                copy!(last(ωsv), first(ωsv))

                ## If `v` is a recombination vertex, the ancestral interval of
                ## one of its parental edge is the intersection of the
                ## appropriate interval associated with the breakpoint with ωv.

                @inbounds for (k, dad) ∈ enumerate(dads(arg, v))
                    e = Edge(dad => v)
                    intersect!(ωsv[k], ancestral_mask(e, arg), buffer = buffer)
                    _update_ai!(vstack, arg, e, ωsv[k], oldhash)
                end
            else
                ## If `v` is a coalescence vertex the ancestral interval of its
                ## parental edges is simply ωv.
                e = Edge(dad(arg, v) => v)

                _update_ai!(vstack, arg, e, first(ωsv), oldhash)
            end
        end
    end

    arg
end

export recombine!
"""
    recombine!(arg, redge, cedge, breakpoint, rlat, clat)

Add a recombination event to an ARG.
"""
function recombine!(arg, redge, cedge, breakpoint, rlat, clat;
                    buffer = default_buffer())
    ## Adjust recombination masks
    if isrecombination(arg, dst(redge)) && src(redge) < otherdad(arg, redge)
        idx = 2recidx(arg, dst(redge)) - 1
        arg.recombination_mask[idx], arg.recombination_mask[idx + 1] =
            arg.recombination_mask[idx + 1], arg.recombination_mask[idx]
    end

    if isrecombination(arg, dst(cedge)) && src(cedge) < otherdad(arg, cedge)
        idx = 2recidx(arg, dst(cedge)) - 1
        arg.recombination_mask[idx], arg.recombination_mask[idx + 1] =
            arg.recombination_mask[idx + 1], arg.recombination_mask[idx]
    end

    ## Add recombination and recoalescence vertices to arg ##
    rvertex, cvertex = nv(arg) .+ (1, 2)
    add_vertices!(
        arg,
        (Sequence(trues(nmarkers(arg))), Sequence(trues(nmarkers(arg)))),
        (rlat, clat))

    ## Replace recombination edge ##
    ωr = ancestral_intervals(arg, redge)
    ωr_left, ωr_right = copy(ωr), copy(ωr)
    intersect!(ωr_left, Ω(0, breakpoint))
    intersect!(ωr_right, Ω(breakpoint, ∞))

    add_edge!(arg, Edge(rvertex, dst(redge)), ωr)
    add_edge!(arg, Edge(src(redge), rvertex), ωr_left)
    rem_edge!(arg, redge)

    ## Replace recoalescence edge ##
    ωc = ancestral_intervals(arg, cedge)
    add_edge!(arg, Edge(cvertex, rvertex), ωr_right)
    add_edge!(arg, Edge(cvertex, dst(cedge)), ωc)
    root_recombination = !rem_edge!(arg, cedge)
    if root_recombination
        arg.mrca[] = cvertex
    else
        ωc_new = union(ωc, ωr_right)
        add_edge!(arg, Edge(src(cedge), cvertex), ωc_new)
    end

    ## Compute sequence of new vertices ##
    @no_escape buffer begin
        mask = @alloc(UInt, div(nmarkers(arg), blocksize(Sequence), RoundUp))
        _compute_sequence!(arg, rvertex, mask)
        _compute_sequence!(arg, cvertex, mask)
    end

    ## Update sequences and ancetral intervals ##
    update_upstream!(arg, src(redge), buffer = buffer)
    update_upstream!(arg, src(cedge), buffer = buffer)

    push!(arg.recombination_mask, AIsType([Ω(0, breakpoint)]))
    push!(arg.recombination_mask, AIsType([Ω(breakpoint, ∞)]))
    arg
end

"""
    extend_recombination!(arg, edge, otherdad, breakpoint; buffer = default_buffer())

Extend ancestral interval of a recombination vertex upstream edges:
* `edge` new interval is [`breakpoint`, ∞);
* the other branch's interval is intersected with [0, `breakpoint`).

This is mainly intended to be used within
[`sample_recombination_constrained!`](@ref).
"""
function extend_recombination!(arg, edge, otherdad, breakpoint; buffer = default_buffer())
    d = dst(edge)

    union!(ancestral_mask(edge, arg), Ω(breakpoint, ∞))

    edge = Edge(otherdad => d)
    intersect!(ancestral_mask(edge, arg), Ω(0, breakpoint))

    update_upstream!(arg, d, buffer = buffer)
end

# -- Edges Iterator ----------------------------------------------------

struct EdgesIntervalRec{I}
    genealogy::Arg
    ωs::I
    buffer::CheapStack{Edge{VertexType}}
    visited::BitVector
    min_latitude::Float64
    breakpoint::Float64
    nextidx::Int
end

function EdgesIntervalRec(arg, ωs, store, breakpoint, nextidx,
                          root = mrca(arg), min_latitude = zero(Float64))
    ##TODO: manage `visited` and `funbuffer` manually.
    eibuffer = CheapStack(store)
    visited = falses(nrecombinations(arg))

    for d ∈ children(arg, root, ωs)
        e = Edge(root => d)
        (breakpoint ∈ ancestral_intervals(arg, e) && !sequence(arg, dst(e))[nextidx]) && continue
        push!(eibuffer, e)
    end

    EdgesIntervalRec(arg, ωs, eibuffer, visited, min_latitude, breakpoint, nextidx)
end

IteratorSize(::EdgesIntervalRec) = Base.SizeUnknown()

eltype(::EdgesIntervalRec) = Edge{VertexType}

function iterate(iter::EdgesIntervalRec, state = 1)
    buffer = iter.buffer
    isempty(buffer) && return nothing

    arg = iter.genealogy
    ωs = iter.ωs
    visited = iter.visited
    min_latitude = iter.min_latitude
    breakpoint = iter.breakpoint
    nextidx = iter.nextidx

    e = pop!(buffer)
    s = dst(e)
    if isrecombination(arg, s)
        ridx = recidx(arg, s)
        visited[ridx] && return e, state + 1
        visited[ridx] = true
    end

    if latitude(arg, s) >= min_latitude
        for d ∈ children(arg, s, ωs)
            newe = Edge(s => d)
            if breakpoint ∈ ancestral_intervals(arg, newe)
                sequence(arg, dst(newe))[nextidx] || continue
            end

            if isrecombination(arg, d)
                breakpoint ∈ ancestral_mask(Edge(s => d), arg) || continue
            end

            push!(buffer, newe)
        end
    end

    e, state + 1
end

function _weight_edge(arg, e_ref, e, mask, h_buf, f)
    chunks = h_buf.data.chunks
    h_ref = sequence(arg, dst(e_ref)).data.chunks
    h = sequence(arg, dst(e)).data.chunks

    @inbounds @simd ivdep for k ∈ eachindex(chunks)
        chunks[k] ⊻= chunks[k]
        chunks[k] ⊻= h_ref[k]
        chunks[k] &= h[k]
        chunks[k] &= mask[k]
    end

    f(h_buf)
end

function _sample_cedge(rng, arg, lat, nextidx, window, live_edge::T, redge, buffer, λ = 0.3) where T
    nextpos = idxtopos(arg, nextidx)

    @no_escape buffer begin
        cedges_ptr = convert(Ptr{T}, @alloc_ptr(ne(arg) * sizeof(T)))
        ws_ptr = convert(Ptr{Float64}, @alloc_ptr(ne(arg) * sizeof(Float64)))
        len = 0

        mask = @alloc(UInt, div(nmarkers(arg), blocksize(Sequence), RoundUp))
        mask!(mask, sam(arg), range(nextidx + 1, min(nextidx + 11, nmarkers(arg))))
        ## TODO: Do that without allocating.
        x = Sequence(undef, nmarkers(arg))

        store = @alloc(T, nv(arg))
        @inbounds for e ∈ EdgesIntervalRec(arg, window, store, nextpos, nextidx, src(live_edge), lat)
            ## Compute distance & weight
            w = _weight_edge(arg, redge, e, mask, x, h -> λ * (1 - λ)^(1 - sum(h)))

            ## Store edge and weight
            len += 1
            unsafe_store!(cedges_ptr, e, len)
            unsafe_store!(ws_ptr, w, len)
        end

        cedges = UnsafeArray{T, 1}(cedges_ptr, (len,))
        ws = UnsafeArray{Float64, 1}(ws_ptr, (len,))
        sample(rng, cedges, ProbabilityWeights(ws))
    end
end

"""
    sample_redge(rng, arg, e, nextidx, rlat_ubound)

Sample a recombination edge conditional on a given edge.

The graph is traversed downstream starting from `dst(e)` as long as there is
only one derived child.
"""
function sample_redge(rng, arg, e, nextidx, rlat_ubound)
    nextpos = idxtopos(arg, nextidx)

    ## Compute the total length of valid branches
    s, d = src(e), dst(e)
    refstate = sequence(arg, d)[nextidx]

    valid_child = d
    total_length = rlat_ubound - latitude(arg, s)
    @inbounds while !iszero(valid_child)
        total_length += branchlength(arg, Edge(s => d))

        valid_child = 0
        for c ∈ children(arg, d, nextpos)
            sequence(arg, c)[nextidx] == refstate || continue
            if iszero(valid_child)
                valid_child = c
            else
                valid_child = 0
            end
        end

        if !iszero(valid_child)
            s = d
            d = valid_child
        end
    end

    ## Sample recombination location
    ## The larger α - β is, the stronger the bias towards ancient branches.
    Δlat_dist = Beta(3, 2)
    Δlat = rand(rng, Δlat_dist)
    arg.logprob[] += logpdf(Δlat_dist, Δlat)
    Δlat *= total_length
    lat = latitude(arg, d) + Δlat

    ## Compute recombination edge ##
    while Δlat > 0
        Δlat -= branchlength(arg, Edge(s => d))
        Δlat > 0 || break
        d = s
        s = (first ∘ dads)(arg, s, nextpos)
    end

    Edge(s => d), lat
end

# -- Constrained recombinations ----------------------------------------

"""
    sample_derived_recombination!(rng, arg, e1, e2, breakpoint, window, nextidx, nextpos, live_edges; buffer = default_buffer())

Sample a recombination & recoalescence event between two derived edges
constrained to eliminate one mutation.

Under certain conditions, it is not necessary to actually sample such an event.
This is the case when the recombination and recoalescence event are incident.
In this case, the recombination edge's ancestral interval is extended (see
[`extend_recombination!`](@ref).
"""
function sample_derived_recombination!(rng, arg, e1, e2,
                                       breakpoint, window, nextidx, nextpos,
                                       live_edges; buffer = default_buffer())
    ## This ensures that there is a least one edge available for recoalescence.
    if latitude(arg, dst(live_edges[e1])) > latitude(arg, dst(live_edges[e2]))
        e1, e2 = e2, e1
    end

    ## Sample recombination edge and latitude
    rlat_ubound = min(latitude(arg, src(live_edges[e1])),
                      latitude(arg, src(live_edges[e2])))
    redge, rlat = sample_redge(rng, arg, popat!(live_edges, e1), nextidx, rlat_ubound)
    if e2 > e1
        ## Accounts for length reduction of ̀`live_edges` following `popat!`
        e2 -= 1
    end

    ## Sample recoalescence edge ##
    cedge = _sample_cedge(rng, arg, rlat, nextidx, window, live_edges[e2],
                          redge, buffer)

    if  src(cedge) ∈ dads(arg, dst(redge)) && isrecombination(arg, dst(redge))
        ## In this situation, there is acutaly no need for a new recombination
        ## event. A mutation can be eliminated simply by extending `redge`'s
        ## ancestral interval.
        extend_recombination!(arg, Edge(src(cedge) => dst(redge)), src(redge),
                              breakpoint, buffer = buffer)

        newd = dst(cedge)
    else
        ## Sample recoalescence latitude ##
        clat_lbound = max(rlat, latitude(arg, dst(cedge)))
        clat_ubound = latitude(arg, src(cedge))
        clat_span = clat_ubound - clat_lbound

        ## Same strategy as for `rlat_dist` (see above).
        clat_dist = Beta(2)
        clat = rand(rng, clat_dist)
        arg.logprob[] += logpdf(clat_dist, clat)
        clat *= clat_span
        clat += clat_lbound

        ## Add recombination event to the graph ##
        @debug "Constrained recombination event" redge cedge breakpoint rlat clat
        recombine!(arg, redge, cedge, breakpoint, rlat, clat, buffer = buffer)

        newd = nv(arg)
    end

    ## Compute new live edge ##
    news = src(cedge)
    if news != src(live_edges[e2])
        while sequence(arg, news)[nextidx]
            newd = news
            news = (first ∘ dads)(arg, news, nextpos)
        end
    end

    for (k, e) ∈ enumerate(live_edges)
        src(e) == news || continue
        live_edges[k] = Edge(news => newd)
        break
    end

    live_edges
end

"""
    sample_recombination_constrained!(rng, arg, breakpoint, winwidth, live_edges; buffer = default_buffer())

Sample a recombination event constrained so as to reduce the number of
mutations for a given marker by one.
"""
function sample_recombination_constrained!(rng, arg, breakpoint, winwidth, live_edges;
                                           buffer = default_buffer())
    n = length(live_edges)
    nextidx = postoidx(arg, breakpoint)
    nextpos = idxtopos(arg, nextidx)

    ## It is necessary to extend the window as below when simulating type II
    ## recombination events.
    window = breakpoint ± winwidth / 2 ∪
        ClosedInterval(idxtopos(arg, nextidx - 1), nextpos)

    ## Sample recombination edge ##
    e1, e2 = sample(rng, eachindex(live_edges), 2; replace = false)
    arg.logprob[] += log(2) - log(n) - log(n - 1)

    sample_derived_recombination!(rng, arg, e1, e2,
                                  breakpoint, window, nextidx, nextpos,
                                  live_edges, buffer = buffer)

    breakpoint
end

function sample_recombination_unconstrained!(rng, arg, winwidth,
                                             buffer = default_buffer())
    ## Sample a recombination edge ##
    @no_escape buffer begin
        redges_ptr = convert(Ptr{Edge{VertexType}},
                             @alloc_ptr(ne(arg) * sizeof(Edge{VertexType})))
        redges_len = 0
        for e ∈ edges(arg)
            (isempty ∘ ancestral_intervals)(arg, e) && continue
            redges_len += 1
            unsafe_store!(redges_ptr, e, redges_len)
        end
        redges = UnsafeArray{Edge{VertexType}, 1}(redges_ptr, (redges_len,))

        ws_redges_data = @alloc(Float64, length(redges))
        map!(e -> branchlength(arg, e), ws_redges_data, redges)
        ws_redges = ProbabilityWeights(ws_redges_data)

        breakpoint = 0
        rlat = 0
        redge = Edge(0 => 0)
        valid_redge = false
        while !valid_redge
            redge_idx = sample(rng, eachindex(redges), ws_redges)
            ws_redges[redge_idx] = 0
            redge = redges[redge_idx]

            ## Sample a breakpoint ##
            bp_int = (closure ∘ ancestral_intervals)(arg, redge) ∩
                Ω(0, (last ∘ positions)(arg))

            breakpoint_dist = Uniform(endpoints(bp_int)...)
            breakpoint = rand(rng, breakpoint_dist)
            arg.logprob[] += logpdf(breakpoint_dist, breakpoint)

            ## Sample the recombination latitude ##
            rlat_dist = Beta(2)
            rlat = rand(rng, rlat_dist)
            arg.logprob[] += logpdf(rlat_dist, rlat)
            rlat *= branchlength(arg, redge)
            rlat += latitude(arg, dst(redge))

            ## Cancel the recombination event if it occurs above the mrca of the
            ## breakpoint.
            if nlive(arg, rlat, breakpoint, buffer = buffer) > 1
                valid_redge = true
            end
        end
    end

    window = breakpoint ± winwidth / 2

    ## Sample the recoalescence latitude ##
    Δclat_dist = Exponential(inv(nlive(arg, rlat, window, buffer = buffer)))
    Δclat = rand(rng, Δclat_dist)
    arg.logprob[] += logpdf(Δclat_dist, Δclat)
    clat = rlat + Δclat

    ## Sample the recoalescence edge ##
    @no_escape buffer begin
        store = @alloc(Edge{VertexType}, ne(arg))
        visited = @alloc(Bool, nrecombinations(arg))
        cedges_data = @alloc(Edge{VertexType}, ne(arg))
        cedges_ptr = firstindex(cedges_data)
        for e ∈ edges_interval(arg, window, store, visited, mrca(arg), clat)
            e == redge && continue
            cedges_data[cedges_ptr] = e
            cedges_ptr += 1
        end
        cedges = view(cedges_data, 1:(cedges_ptr-1))

        ws_cedges = @alloc(Float64, length(cedges))
        map!(e -> latitude(arg, dst(e)) < clat < latitude(arg, src(e)),
             ws_cedges, cedges)
        sum_ws = sum(ws_cedges)

        if clat > tmrca(arg)
            cedge = Edge(mrca(arg) => mrca(arg))
        else sum_ws > 0
            cedge = sample(rng, cedges, ProbabilityWeights(ws_cedges, sum_ws))
        end
    end

    ## Add recombination event to the graph ##
    @debug "Unconstrained recombination event" redge cedge breakpoint rlat clat
    recombine!(arg, redge, cedge, breakpoint, rlat, clat, buffer = buffer)

    breakpoint
end

function build!(rng, arg::Arg, ρ; winwidth = ∞, buffer = default_buffer())
    ## Unconstrained recombinations ##
    nrecs_dist = Poisson(ρ)
    nrecs = rand(rng, nrecs_dist)
    arg.logprob[] += logpdf(nrecs_dist, nrecs)

    @no_escape buffer begin
        for k ∈ 1:nrecs
            sample_recombination_unconstrained!(rng, arg, winwidth, buffer)
        end
    end

    ## Constrained recombinations ##
    @no_escape buffer begin
        mutation_edges_buffer = @SVector [Edge{VertexType}[] for _ ∈ 1:64]
        nextidx, live_edges = next_inconsistent_idx(arg, 1,
                                                    mutations_edges = mutation_edges_buffer,
                                                    buffer = buffer)

        while !iszero(nextidx)
            nbp = length(live_edges) - 1

            bp_lbound = isone(nextidx) ?
                zero(eltype(positions(arg))) : idxtopos(arg, nextidx - 1)
            bp_ubound = idxtopos(arg, nextidx)

            bp_dist = Uniform(bp_lbound, bp_ubound)
            breakpoints = (sort ∘ rand)(rng, bp_dist, nbp)
            arg.logprob[] -= nbp * log(bp_ubound - bp_lbound)

            for breakpoint ∈ breakpoints
                sample_recombination_constrained!(rng, arg, breakpoint,
                                                  winwidth, live_edges,
                                                  buffer = buffer)
            end

            previdx = nextidx
            nextidx, live_edges = next_inconsistent_idx(arg, nextidx + 1,
                                                        mutations_edges = mutation_edges_buffer,
                                                        buffer = buffer)
        end
    end

    arg
end
#
###########
# Algebra #
###########

"""
    _path_bfs_forward!(estack, argm s, d, vqueue, visited)

Forward pass of the path finding algorithm.
"""
function _path_bfs_forward!(estack, arg::Arg, s, d, vqueue, visited)
    n = nleaves(arg)
    empty!(estack)
    empty!(vqueue)
    empty!(visited)

    push!(visited, s)
    push!(vqueue.store, s) # Dafuq??

    while !isempty(vqueue)
        s = popfirst!(vqueue.store) # DAfuq?

        for v ∈ children(arg, s)
            v ∈ visited && continue

            ## Make sure that the edge is in the spanning tree. If `s` is not
            ## the largest parent of v, skip the edge s->v. This means the
            ## parental edge of a recombination vertex not adjacent to a
            ## recoalescence vertex is excluded from the spanning tree.
            if isrecombination(arg, v, n)
                any(>(s), dads(arg, v)) && continue
            end

            push!(visited, v)
            push!(estack, Edge(s, v))
            v == d && @goto done
            push!(vqueue.store, v) # Dafuq?
        end

        _dads = dads(arg, s)
        isempty(_dads) && continue

        ## Again, skip non-recoalescence edge.
        v = maximum(dads(arg, s))
        v ∈ visited && continue
        push!(visited, v)
        push!(estack, Edge(s, v))
        v == d && @goto done
        push!(vqueue.store, v) # Dafuq?
    end

    @label done
    estack
end

"""
    _update_vec!(vec, edgesid, e, lk)

Updates the variable `vec` in `_bfs_backtrack!.
"""
function _update_vec!(vec, edgesid, e, lk)
    val = 1
    idx = get(edgesid, e, zero(Int))

    if iszero(idx)
        idx = get(edgesid, reverse(e), zero(Int))
        val = -1
    end
    Threads.lock(lk) do
        vec[idx] = val
    end
end

function _bfs_backtrack!(vec, edgesid, estack, lk)
    e_prev = pop!(estack)
    _update_vec!(vec, edgesid, e_prev, lk)
    @inbounds while !isempty(estack)
        e = pop!(estack)
        dst(e) == src(e_prev) || continue

        _update_vec!(vec, edgesid, e, lk)
        e_prev = e
    end

    vec
end

export cbasis!, cbasis
"""
    cbasis!(vec, arg, v, lk; estack, edgesid, vqueue, visited)
    cbasis(arg, v, lk; estack, edgesid, vqueue, visited)
    cbasis!(mat, arg; edgesid)
    cbasis(arg; edgesid)

Compute the basis vector associated with recombination vertex `v`. If `v` is
not specified, return a matrix containing all basis vectors.
"""
function cbasis! end, function cbasis end

function cbasis!(vec, arg::Arg, v::VertexType, lk = Threads.ReentrantLock();
    estack = Stack{Edge{VertexType}}(ceil(Int, log(nv(arg)))),
    edgesid = edgesmap(arg),
    vqueue = Queue{VertexType}(ceil(Int, log(nv(arg)))),
    visited = Set{VertexType}())
    empty!(estack)
    empty!(vqueue)
    empty!(visited)

    ## Add edges upstream of `v` to vec.
    _dads = minmax(dads(arg, v)...)
    Threads.lock(lk) do
        vec[edgesid[Edge(first(_dads) => v)]] = -1
        vec[edgesid[Edge(last(_dads) => v)]] = 1
    end

    _path_bfs_forward!(estack, arg, _dads..., vqueue, visited)
    _bfs_backtrack!(vec, edgesid, estack, lk)
end

function cbasis(arg::Arg, v::VertexType, lk = Threads.ReentrantLock();
    estack = Stack{Edge{VertexType}}(ceil(Int, log(nv(arg)))),
    edgesid = edgesmap(arg),
    vqueue = Queue{VertexType}(ceil(Int, log(nv(arg)))),
    visited = Set{VertexType}())
    cbasis!(spzeros(Float64, ne(arg)), arg, v, lk,
        estack = estack, edgesid = edgesid,
        vqueue = vqueue, visited = visited)
end

function cbasis!(mat, arg::Arg; edgesid = edgesmap(arg))
    fill!(mat, 0)

    r = nrecombinations(arg)
    n = nleaves(arg)

    lk = Threads.ReentrantLock()
    rec_offset = ne(arg) - nv(arg) + 1 - nrecombinations(arg)

    tasks_chunks = chunks(range(1, length = r - rec_offset),
        n = Threads.nthreads(), split = :scatter)
    tasks = map(tasks_chunks) do ks
        Threads.@spawn begin
            local estack = Stack{Edge{VertexType}}(ceil(Int, log(nv(arg))))
            local vqueue = Queue{VertexType}(ceil(Int, log(nv(arg))))
            local visited = Set{VertexType}()

            for k ∈ ks
                v = 2(n + k + rec_offset - 1)
                cbasis!(view(mat, :, k), arg, v, lk,
                    edgesid = edgesid, vqueue = vqueue,
                    estack = estack, visited = visited)
            end
        end
    end

    fetch.(tasks)
    mat
end

cbasis(arg::Arg; edgesid = edgesmap(arg)) =
    cbasis!(spzeros(ne(arg), nrecombinations(arg)), arg, edgesid = edgesid)

"""
    _impedance_Z(arg, edgesmap, take_sqrt = true)

Diagonal matrix of edges impedance.
"""
function _impedance_Z(arg::Arg, edgesmap, take_sqrt = true)
    Z = Diagonal(Vector{Float64}(undef, ne(arg)))

    for (e, idx) ∈ edgesmap
        Δ = latitude(arg, src(e)) - latitude(arg, dst(e))
        ## Dirty trick to deal with numerical instability. Some impedances are
        ## <= 0...
        Z[idx, idx] = max(eps(Float64), Δ)
    end

    take_sqrt || return Z

    for k ∈ eachindex(Z.diag)
        Z.diag[k] = sqrt(Z.diag[k])
    end

    Z
end

"""
    _impedance_C(arg, p, edgesmap)

Return the generator of the cycle space of an ancestral recombination graph
expanded to include space for a pseudo edge. If `arg` has `r` recombinations and
`k` edges, a k x (r + p) matrix is returned. The generator is stored in
the upper k x r block. Remaining entries are initialized to 0.
"""
function _impedance_C(arg, p, edgesmap)
    r = nrecombinations(arg)
    k = ne(arg)

    C = spzeros(k, r + p)
    cbasis!(C, arg, edgesid = edgesmap)
    C
end

function _update_C_cycles!(C, arg, vs, k, edgesid, estack, vqueue, visited)
    r = nrecombinations(arg)
    v1, iter = Iterators.peel(vs)

    for v ∈ iter
        vec = @view C[:, r + k]
        fill!(vec, 0)

        _path_bfs_forward!(estack, arg, v1, v, vqueue, visited)
        _bfs_backtrack!(vec, edgesid, estack, Threads.ReentrantLock())

        k += 1
    end

    k
end

## TODO: parallelize?
"""
    _impedance_update_C!(C, arg, ss, ds, edgesmap, estack, vqueue, visited)

Store the pseudo-cycles generated by the pairs (d_1, d_2), ..., (d_1, s_pI) in
`C`. pI is the number of input nodes. Pseudo-cycle `k` is stored in column
`r + k`.
"""
function _impedance_update_C!(C, arg, ss, ds, edgesid, estack, vqueue, visited)
    r = nrecombinations(arg)
    k = 1

    ## Sources cycles ##
    k = _update_C_cycles!(C, arg, ss, k, edgesid, estack, vqueue, visited)

    ## Destination cycles ##
    k = _update_C_cycles!(C, arg, ds, k, edgesid, estack, vqueue, visited)

    ## Generator cycle ##
    vec = @view C[:, r + k]
    fill!(vec, 0)

    _path_bfs_forward!(estack, arg, first(ss), first(ds), vqueue, visited)
    _bfs_backtrack!(vec, edgesid, estack, Threads.ReentrantLock())

    C
end

export impedance!, impedance
"""
impedance(arg[, ss, ds]; estack, edgesmap, vqueue, visited)
impedance!(arg, ss, ds C, Z2; estack, edgesmap, vqueue, visited)
impedance_matrix(arg, estack, edgesmap)

Compute impedances. A set of source and destination vertices can be provided
through arguments `ss` and `ds` respectively. Otherwise, impedances are
computed pairwise between leaves.
"""
function impedance!(arg::Arg, ss, ds, C, Z2;
    edgesmap = edgesmap(arg),
    estack = Stack{Edge{VertexType}}(ceil(Int, log(nv(arg)))),
    vqueue = Queue{VertexType}(ceil(Int, log(nv(arg)))),
    visited = Set{VertexType}())
    empty!(estack)
    empty!(vqueue)
    empty!(visited)

    _impedance_update_C!(C, arg, ss, ds, edgesmap, estack, vqueue, visited)

    ## TODO: work out the ordering thing.
    U = (UpperTriangular ∘ Matrix)(qr(Z2 * C, ordering = 0).R)
    X = vcat(zeros((first ∘ size)(U) - 1), 1.0)
    ldiv!(U, ldiv!(U', X))
    w = last(X)
    inv(w)
end

function impedance(arg::Arg, ss, ds;
    estack = Stack{Edge{VertexType}}(ceil(Int, log(nv(arg)))),
    edgesmap = edgesmap(arg),
    vqueue = Queue{VertexType}(ceil(Int, log(nv(arg)))),
    visited = Set{VertexType}())
    empty!(estack)
    empty!(vqueue)
    empty!(visited)

    p = length(ss) + length(ds)

    ## Square root of impedances.
    Z2 = _impedance_Z(arg, edgesmap, true)

    ## Compute fundamental basis
    C = _impedance_C(arg, p - 1, edgesmap)

    impedance!(arg, ss, ds, C, Z2,
        edgesmap = edgesmap, estack = estack,
        vqueue = vqueue, visited = visited)
end

function impedance(arg::Arg;
    estack = Stack{Edge{VertexType}}(ceil(Int, log(nv(arg)))),
    edgesmap = edgesmap(arg),
    vqueue = Queue{VertexType}(ceil(Int, log(nv(arg)))),
    visited = Set{VertexType}())
    empty!(estack)

    C = _impedance_C(arg, 1, edgesmap)
    Z2 = _impedance_Z(arg, edgesmap, true)

    Iterators.map(combinations(1:nleaves(arg), 2)) do (s, d)
        impedance!(arg, s, d, C, Z2,
            estack = estack, edgesmap = edgesmap,
            vqueue = vqueue, visited = visited)
    end
end

export impedance_matrix
function impedance_matrix(arg::Arg,
    estack = Stack{Edge{VertexType}}(ceil(Int, log(nv(arg)))),
    edgesmap = edgesmap(arg))
    empty!(estack)

    n = nleaves(arg)

    C = _impedance_C(arg, 1, edgesmap)
    Z2 = _impedance_Z(arg, edgesmap, true)

    mat = Matrix{Float64}(undef, n, n)
    mat[diagind(mat)] .= 0
    for (j, i) ∈ combinations(1:n, 2)
        mat[i, j] =
            impedance!(arg, (i,), (j,), C, Z2, estack = estack, edgesmap = edgesmap)
    end

    Symmetric(mat, :L)
end

#          +----------------------------------------------------------+
#          |                      Thévenin Tree                       |
#          +----------------------------------------------------------+

export thevenin!, thevenin
"""
    thevenin!(rng, tree, arg)
    thevenin(rng, arg)

Sample a Thévenin tree conditional on an ARG.
"""
function thevenin! end, function thevenin end

function _thevenin_impedance_helper!(tree, arg, v1_leaves, v2, C, Z2,
    edgesmap, estack, vqueue, visited)
    r = nrecombinations(arg)
    v2_leaves = descendants(tree, v2) ∩ leaves(tree)
    if isempty(v2_leaves)
        push!(v2_leaves, v2)
    end

    p = length(v1_leaves) + length(v2_leaves)
    C_view = @view C[:, range(1, r + p - 1)]
    impedance!(arg, v1_leaves, v2_leaves, C_view, Z2,
        edgesmap = edgesmap, estack = estack,
        vqueue = vqueue, visited = visited)
end

function thevenin!(tree, arg::Arg; ϵ = 1e-5)
    n = nleaves(arg)
    r = nrecombinations(arg)

    ## Initialize tree ##
    empty!(graph(tree).fadjlist)
    graph(tree).ne = 0
    add_vertices!(graph(tree), 2n - 1)

    resize!(latitudes(tree), n - 1)

    resize!(sequences(tree), 2n - 1)
    sequences(tree)[1:n] .= sequences(arg, leaves(arg))

    tree.logprob[] = 0

    # ## Build tree ##
    edgesid = edgesmap(arg)
    estack = Stack{Edge{VertexType}}(ceil(Int, log(nv(arg))))
    vqueue = Queue{VertexType}(ceil(Int, log(nv(arg))))
    visited = Set{VertexType}()

    ## To limit allocations, we make `C` large enough to handle the worst case
    ## scenario. We can pass an appropriate `view` to `impedance!`.
    C = _impedance_C(arg, n - 1, edgesid)
    Z2 = _impedance_Z(arg, edgesid, true)

    live = collect(leaves(tree))

    ## Distance matrix. Used for clustering
    dists = impedance_matrix(arg, estack, edgesid)
    for k ∈ 1:n
        dists[k, k] = ∞
    end

    ## Descendant leaves of live vertices.
    dleaves = Vector{Set{VertexType}}(undef, n)
    for k ∈ eachindex(dleaves)
        dleaves[k] = Set(k)
    end

    ## Impedance between live vertices and their descendant leaves
    live_impedances = zeros(Float64, n)

    for nlive ∈ range(length(live), 2, step = -1)
        ## Find the next coalescing pair ##
        ii, jj = Tuple(findfirst(!isinf, dists))
        @inbounds for ii_it ∈ 1:(nlive-1)
            for jj_it ∈ (ii_it+1):nlive
                dists[ii_it, jj_it] <= dists[ii, jj] || continue
                ii, jj = ii_it, jj_it
            end
        end

        ## Create coalescence event ##
        i, j = live[ii], live[jj]
        k = 2n + 1 - nlive

        add_edge!(tree, k, i)
        add_edge!(tree, k, j)

        newlat = 0.5 * (dists[ii, jj] +
            latitude(tree, i) - live_impedances[ii] +
            latitude(tree, j) - live_impedances[jj])
        Δlat = newlat - latitude(tree, k - 1)
        if Δlat < ϵ
            newlat = latitude(tree, k - 1) + ϵ
        end
        latitudes(tree)[k - n] = newlat

        ## Compute descendant leaves ##
        union!(dleaves[ii], dleaves[jj])

        ## Update `live_impedances` ##
        live_impedances[ii] = inv(
            inv(live_impedances[ii] + branchlength(tree, k, i)) +
            inv(live_impedances[jj] + branchlength(tree, k, j)))

        ## Update live vertices and impedance matrix ##
        ## Replace i by k.
        live[ii] = k

        ## Make the distance from any vertex to j infinite.
        @inbounds for ii_it ∈ 1:(jj-1)
            dists.data[jj, ii_it] = ∞
        end
        @inbounds for ii_it ∈ (jj+1):n
            dists.data[ii_it, jj] = ∞
        end

        ## Compute distances from k.
        leaves_ii = dleaves[ii]
        n_ii = length(leaves_ii)

        @inbounds for jj_it ∈ 1:(ii-1)
            isinf(dists.data[ii, jj_it]) && continue
            leaves_jj_it = dleaves[jj_it]
            n_jj_it = length(leaves_jj_it)

            C_view = view(C, :, 1:(r + n_ii + n_jj_it - 1))
            dists.data[ii, jj_it] = impedance!(arg, leaves_ii, leaves_jj_it, C_view, Z2,
                                               edgesmap = edgesid, estack = estack,
                                               vqueue = vqueue, visited = visited)
        end

        @inbounds for jj_it ∈ (ii+1):n
            isinf(dists.data[jj_it, ii]) && continue
            leaves_jj_it = dleaves[jj_it]
            n_jj_it = length(leaves_jj_it)

            C_view = view(C, :, 1:(r + n_ii + n_jj_it - 1))
            dists.data[jj_it, ii] = impedance!(arg, leaves_jj_it, leaves_ii, C_view, Z2,
                                               edgesmap = edgesid, estack = estack,
                                               vqueue = vqueue, visited = visited)
        end
    end

    tree
end

thevenin(arg::Arg) = thevenin!(Tree(sam(arg)), arg)

function validate(arg::Arg; check_mutations = true)
    flag = true
    n = nleaves(arg)

    ## General properties ##
    if check_mutations && nmutations(arg) != nmarkers(arg)
        @info "Number of mutations not equal to the number of markers"
        flag = false
    end

    ## Leaves ##
    for v ∈ leaves(arg)
        if (!iszero ∘ length)(children(arg, v))
            @info "Leaf with children" v
            flag = false
        end

        if (!isone ∘ length)(dads(arg, v))
            @info "Leaf with number of parents not equal to 1" v
            flag = false
        end

        if ancestral_intervals(arg, v) != AIsType([Ω(0, ∞)])
            @info "Ancestral interval of leaf's parental edge not equal to [0, ∞)" v
            flag = false
        end
    end

    ## Non-root coalescence vertices ##
    for v ∈ Iterators.flatten((range(n + 1, 2n - 1), range(2n + 1, nv(arg), step = 2)))
        v == mrca(arg) && continue

        if length(children(arg, v)) != 2
            @info "Coalescence vertex with number of children not equal to 2" v
            flag = false
        end

        if (!isone ∘ length)(dads(arg, v))
            @info "Coalescence vertex with number of parents not equal to 1" v
            flag = false
        end

        ai_children = mapreduce(x -> ancestral_intervals(arg, Edge(v => x)), ∪, children(arg, v))
        if ancestral_intervals(arg, Edge(dad(arg, v) => v)) != ai_children
            msg = "Coalescence vertex whose parental edge's ancestral interval" *
                " is not equal to the union of the ancestral intervals of its" *
                " children's edge."
            @info msg v
            flag = false
        end
    end

    ## Recombination vertices ##
    for v ∈ range(2n, nv(arg), step = 2)
        if (!isone ∘ length)(children(arg, v))
            @info "Recombination vertex with number of children not equal to 1" v
            flag = false
        end

        if length(dads(arg, v)) != 2
            @info "Recombination vertex with number of parents not equal to 2" v
            flag = false
        end

        # ai_child = ancestral_intervals(arg, Edge(v => child(arg, v)))
        # ai_right = ancestral_intervals(arg, Edge(rightdad(arg, v) => v))
        # ai_left = ancestral_intervals(arg, Edge(leftdad(arg ,v) => v))
        # if ai_child != ai_left ∪ ai_right
        #     msg = "Recombination vertex whose children edge's ancestral" *
        #         " interval is not equal to the union of it parental edges's" *
        #         " ancestral intervals."
        #     @info msg v ai_child ai_left ai_right
        #     flag = false
        # end
    end

    ## Internal vertices ##
    for v ∈ ivertices(arg)
        ref = mapreduce(&, children(arg, v)) do _child
            sequence(arg, _child) | ~ancestral_mask(arg, Edge(v => _child))
        end
        ref &= ancestral_mask(arg, v)

        if sequence(arg, v) != ref
            @info "Inconsistent sequence" v sequence(arg, v) ref
            flag = false
        end

        if latitude(arg, v) < maximum(u -> latitude(arg, u), children(arg, v))
            @info "Inconsistent latitudes" v
            flag = false
        end
    end

    ## Edges ##
    # for e ∈ edges(arg)
    #     if isempty(ancestral_intervals(arg, e))
    #         @info "Empty ancestral interval" e
    #         flag = false
    #     end
    # end

    flag
end
