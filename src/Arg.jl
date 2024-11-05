using Graphs

import Graphs: add_vertices!, add_edge!, rem_edge!

using Random

using StatsBase: samplepair, ProbabilityWeights

using SparseArrays

using Combinatorics: combinations

using Distributions

##################
# Arg Definition #
##################

export Arg
struct Arg <: AbstractGenealogy
    graph::SimpleDiGraph{VertexType}
    latitudes::Vector{Float64}
    breakpoints::Vector{Float64}
    rightdads::Vector{VertexType}
    mrca::Base.RefValue{VertexType}
    sequences::Vector{Sequence}
    ancestral_intervals::Dict{Edge{VertexType}, Set{Ω}}
    sample::Sample
    logprob::Base.RefValue{Float64x2}
end

Arg(tree::Tree) = Arg(
    graph(tree),
    latitudes(tree),
    Vector{Float64}(undef, 0),
    Vector{VertexType}(undef, 0),
    Ref(mrca(tree)),
    sequences(tree),
    Dict{Edge{VertexType}, Set{Ω}}(),
    sam(tree),
    Ref(prob(tree, logscale = true))
)

###############################
# AbstractGenealogy Interface #
###############################

nleaves(arg::Arg) = length(sam(arg).H)

describe(::Arg, long = true) = long ? "Ancestral Recombination Graph" : "ARG"

function isrecombination(::Arg, v, n)
    v < 2n && return false
    isodd(v) && return false
    true
end

isrecombination(arg::Arg, v) = isrecombination(arg, v, nleaves(arg))

function recombinations(arg::Arg; dummy = false)
    isempty(arg.breakpoints) && return StepRange{Int, Int}(0, 1, 0)

    start = 2nleaves(arg)
    step = 2
    stop = nv(arg)

    if !dummy
        while isinf(recbreakpoint(arg, start))
            start += step
        end
    end

    StepRange{Int, Int}(start, step, stop)
end

nrecombinations(arg::Arg) = ne(arg) - nv(arg) + 1

mrca(arg::Arg) = arg.mrca[]

mrca(arg, ωs) = mrca(arg, leaves(arg), ωs)

###########################
# AbstractGraph Interface #
###########################

function add_vertices!(arg::Arg, H, lats)
    append!(latitudes(arg), lats)
    append!(sequences(arg), H)
    add_vertices!(graph(arg), length(H))
end

function add_edge!(arg::Arg, e, ints::Set{Ω})
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

function ancestral_intervals!(ωs, arg::Arg, e::Edge; wipe = true)
    wipe && empty!(ωs)

    haskey(arg.ancestral_intervals, e) || return push!(ωs, Ω(0, ∞))

    @inline for ω ∈ arg.ancestral_intervals[e]
        push!(ωs, ω)
    end

    ωs
end

## Use this `get!` method so the interval doesn't get constructed every call.
ancestral_intervals(arg::Arg, e::Edge) =
    get!(() -> Set{Ω}((Ω(0, ∞),)), arg.ancestral_intervals, e)

function ancestral_intervals!(ωs, arg::Arg, v::VertexType;
                              wipe = true, buffer = default_buffer())
    wipe && empty!(ωs)
    isleaf(arg, v) && return push!(ωs, Ω(0, ∞))

    for child ∈ children(arg, v)
        ancestral_intervals!(ωs, arg, Edge(v => child), wipe = false)
    end

    simplify!(ωs)
end

ancestral_intervals(arg::Arg, v::VertexType) = ancestral_intervals!(Set{Ω}(), arg, v)

ancestral_mask!(η, arg::Arg, x::Union{VertexType, Edge{VertexType}};
                ωs_buf = Set{Ω}(), wipe = true) =
    ancestral_mask!(η, sam(arg), ancestral_intervals!(ωs_buf, arg, x), wipe = wipe)

ancestral_mask(arg::Arg, x::Union{VertexType, Edge{VertexType}};
               ωs_buf = Set{Ω}()) =
    ancestral_mask!(Sequence(falses(nmarkers(arg))), arg, x,
                    ωs_buf = ωs_buf, wipe = false)

recidx(arg, v) = (v - 2(nleaves(arg) - 1)) ÷ 2

export recbreakpoint
"""
    recbreakpoint(arg, v)

Returns the breakpoint associated with vertex `v`. If `v` is not a
recombination vertex, returns ∞.
"""
function recbreakpoint end

export rightdad
"""
    rightdad(arg, v)
    leftdad(arg, v)

Returns the parent of recombination vertex `v` ancestral for material
to the left/right of the breakpoint associated with `v`. If `v` is not a
recombination vertex, returns 0.
"""
function rightdad end

for (fun, (def, field)) ∈ Dict(:recbreakpoint => (:∞, Meta.quot(:breakpoints)),
                                 :rightdad => (:(zero(VertexType)), Meta.quot(:rightdads)))
    @eval function $fun(arg::Arg, v)
        ret = $def
        isrecombination(arg, v) || return ret
        getfield(arg, $field)[recidx(arg, v)]
    end
end

function leftdad(arg, v)
    r = rightdad(arg, v)
    iszero(r) && return r

    for d ∈ dads(arg, v)
        r != d && return d
    end

    r
end

#######
# MMN #
#######

function mutationsidx!(res, mask, arg, e, firstchunk, firstidx, lastchunk;
                       ωs_buf = Set{Ω}())
    η1, η2 = sequences(arg, e)
    ancestral_mask!(mask, arg, e, ωs_buf = ωs_buf)
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
    ωs_buf = Set{Ω}()
    @no_escape buffer begin
        store = @alloc(Edge{VertexType}, nleaves(arg) + nrecombinations(arg))
        @inbounds for e ∈ edges_interval(arg, ω, store)
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

##################
# Recombinations #
##################

function _compute_sequence!(arg, v, mask; ωs_buf = Set{Ω}())
    η = sequence(arg, v)
    ## Set every markers to 1 ##
    η.data.chunks .⊻= .~η.data.chunks

    @inbounds for child ∈ children(arg, v)
        ancestral_mask!(mask, arg, Edge(v, child), ωs_buf = ωs_buf)
        η.data.chunks .&= sequence(arg, child).data.chunks .| .~mask.data.chunks
    end

    ancestral_mask!(mask, arg, v, ωs_buf = ωs_buf)
    η.data.chunks .&= mask.data.chunks
    η
end

function _update_ai!(vstack, arg, e, ωs, sequence_oldhash)
    ωs_oldhash = (hash ∘ ancestral_intervals)(arg, e)

    ## Update ancestral interval of e
    if haskey(arg.ancestral_intervals, e)
        empty!(arg.ancestral_intervals[e])
        while !isempty(ωs)
            push!(arg.ancestral_intervals[e], pop!(ωs))
        end
    else
        copy!(arg.ancestral_intervals[e], ωs)
    end

    ## Update vstack
    if (hash ∘ ancestral_intervals)(arg, e) != ωs_oldhash ||
        (hash ∘ sequence)(arg, dst(e)) != sequence_oldhash
        push!(vstack, src(e))
    end

    ancestral_intervals(arg, e)
end

function update_upstream!(arg, v; buffer = default_buffer())
    @no_escape buffer begin
        store = @alloc(VertexType, nrecombinations(arg) + 1)
        vstack = CheapStack(store)
        push!(vstack, v)

        mask = Sequence(undef, nmarkers(arg))
        ωsv = Set{Ω}()

        while !isempty(vstack)
            v = pop!(vstack)

            ## Update sequence of `v` ##
            sequence_oldhash = (hash ∘ sequence)(arg, v)
            _compute_sequence!(arg, v, mask, ωs_buf = ωsv)
            iszero(dad(arg, v)) && continue

            ## Update ancestral intervals of parental edges ##
            ancestral_intervals!(ωsv, arg, v, buffer = buffer)

            if isrecombination(arg, v)
                ## If `v` is a recombination vertex, the ancestral interval of
                ## one of its parental edge is the intersection of the
                ## appropriate interval associated with the breakpoint with ωv.
                bp = recbreakpoint(arg, v)

                ## Left dad
                e = Edge(leftdad(arg, v) => v)
                _update_ai!(vstack, arg, e, ωsv ∩ Ω(0, bp), sequence_oldhash)

                ## Right dad
                e = Edge(rightdad(arg, v) => v)
                intersect!(ωsv, Ω(bp, ∞))
                _update_ai!(vstack, arg, e, ωsv, sequence_oldhash)
            else
                ## If `v` is a coalescence vertex the ancestral interval of its
                ## parental edges is simply ωv.
                e = Edge(dad(arg, v) => v)

                _update_ai!(vstack, arg, e, ωsv, sequence_oldhash)
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

    ## Add recombination and recoalescence vertices to arg ##
    rvertex, cvertex = nv(arg) .+ (1, 2)
    add_vertices!(
        arg,
        (Sequence(trues(nmarkers(arg))), Sequence(trues(nmarkers(arg)))),
        (rlat, clat))

    ## Replace recombination edge ##
    ωr = ancestral_intervals(arg, redge)
    ωr_left = intersect(ωr, Ω(0, breakpoint), buffer = buffer)
    ωr_right = intersect(ωr, Ω(breakpoint, ∞), buffer = buffer)

    add_edge!(arg, Edge(rvertex, dst(redge)), ωr)
    add_edge!(arg, Edge(src(redge), rvertex), ωr_left)
    rem_edge!(arg, redge)
    if isrecombination(arg, dst(redge)) && src(redge) == rightdad(arg, dst(redge))
        arg.rightdads[recidx(arg, dst(redge))] = rvertex
    end

    ## Replace recoalescence edge ##
    ωc = ancestral_intervals(arg, cedge)
    add_edge!(arg, Edge(cvertex, rvertex), ωr_right)
    add_edge!(arg, Edge(cvertex, dst(cedge)), ωc)
    root_recombination = !rem_edge!(arg, cedge)
    if root_recombination
        arg.mrca[] = cvertex
    else
        ωc_new = union(ωc, ωr_right, buffer = buffer)
        add_edge!(arg, Edge(src(cedge), cvertex), ωc_new)
    end
    if isrecombination(arg, dst(cedge)) && src(cedge) == rightdad(arg, dst(cedge))
        arg.rightdads[recidx(arg, dst(cedge))] = cvertex
    end

    ## Compute sequence of new vertices ##
    let mask = Sequence(undef, nmarkers(arg)), ωs_buf = Set{Ω}()
        _compute_sequence!(arg, rvertex, mask, ωs_buf = ωs_buf)
        _compute_sequence!(arg, cvertex, mask, ωs_buf = ωs_buf)
    end

    ## Update sequences and ancetral intervals ##
    update_upstream!(arg, src(redge), buffer = buffer)
    # root_recombination || update_upstream!(arg, src(cedge), buffer = buffer)
    update_upstream!(arg, src(cedge), buffer = buffer)

    push!(arg.rightdads, cvertex)
    push!(arg.breakpoints, breakpoint)
    arg
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
                breakpoint < recbreakpoint(arg, d) && s == rightdad(arg, d) &&
                continue
                breakpoint >= recbreakpoint(arg, d) && s == leftdad(arg, d) &&
                continue
            end

            push!(buffer, newe)
        end
    end

    e, state + 1
end

function _sample_cedge(rng, arg, lat, nextidx, window, live_edge::T, redge, buffer, λ = 0.3) where T
    nextpos = idxtopos(arg, nextidx)

    @no_escape buffer begin
        cedges_ptr = convert(Ptr{T}, @alloc_ptr(ne(arg) * sizeof(T)))
        ws_ptr = convert(Ptr{Float64}, @alloc_ptr(ne(arg) * sizeof(Float64)))
        len = 0

        mask = Sequence(falses(nmarkers(arg)))
        mask.data[range(nextidx + 1, min(nextidx + 11, nmarkers(arg)))] .= true
        x = Sequence(undef, nmarkers(arg))

        store = @alloc(T, nv(arg))
        @inbounds for e ∈ EdgesIntervalRec(arg, window, store, nextpos, nextidx, src(live_edge), lat)
            ## Compute distance & weight
            x.data .⊻= x.data
            x.data .⊻= sequence(arg, dst(redge)).data .& sequence(arg, dst(e)).data
            x.data .&= mask.data
            w = λ * (1 - λ)^(1 - sum(x))

            ## Store edge and weight
            len += 1
            unsafe_store!(cedges_ptr, e, len)
            unsafe_store!(ws_ptr, w, len)
        end

        cedges = unsafe_wrap(Array, cedges_ptr, len)
        ws = unsafe_wrap(Array, ws_ptr, len)
        sample(rng, cedges, ProbabilityWeights(ws))
    end
end

function sample_redge(rng, arg, e, nextidx)
    nextpos = idxtopos(arg, nextidx)
    total_length = branchlength(arg, e)
    s, d = src(e), dst(e)

    ## Total length of valid branches ##
    valid_child = 0
    _children = children(arg, d, nextpos)
    @inbounds for child ∈ _children
        sequence(arg, child)[nextidx] || continue
        if iszero(valid_child)
            valid_child = child
        else
            valid_child = 0
        end
    end

    @inbounds while !iszero(valid_child)
        s = d
        d = valid_child
        total_length += branchlength(arg, Edge(s => d))
        _children = children(arg, d, nextpos)
        if (isone ∘ length)(_children)
            valid_child = first(_children)
        else
            valid_child = 0
        end
    end

    ## Sample recombination location ##
    rlat_dist = Beta(2)
    rlat = rand(rng, rlat_dist)
    arg.logprob[] += logpdf(rlat_dist, rlat)
    rlat *= total_length

    ## Compute recombination edge ##
    rlat -= branchlength(arg, Edge(s => d))
    while rlat > 0
        rlat -= branchlength(arg, Edge(s => d))
        rlat > 0 || break
        d = s
        s = (first ∘ dads)(arg, s, nextpos)
    end

    Edge(s => d)
end

function sample_recombination_constrained!(rng, arg, breakpoint, winwidth, live_edges;
                                           buffer = default_buffer())
    n = length(live_edges)
    nextidx = postoidx(arg, breakpoint)
    nextpos = idxtopos(arg, nextidx)

    ## Sample recombination edge ##
    e1, e2 = sample(rng, eachindex(live_edges), 2; replace = false)
    arg.logprob[] += log(2) - log(n) - log(n - 1)

    ## This ensures that there is a least one edge available for recoalescence.
    if latitude(arg, dst(live_edges[e1])) > latitude(arg, dst(live_edges[e2]))
        e1, e2 = e2, e1
    end

    redge = sample_redge(rng, arg, popat!(live_edges, e1), nextidx)
    if e2 > e1
        e2 -= 1
    end

    ## It is necessary to extend the window as below when simulating type II
    ## recombination events.
    window = breakpoint ± winwidth / 2 ∪
        ClosedInterval(idxtopos(arg, nextidx - 1), nextpos)

    ## Sample recombination latitude ##
    rlat_lbound = latitude(arg, dst(redge))
    rlat_ubound = min(latitude(arg, src(redge)), latitude(arg, src(live_edges[e2])))
    rlat_span = rlat_ubound - rlat_lbound
    rlat_dist = Beta(2)
    rlat = rand(rng, rlat_dist)
    arg.logprob[] += logpdf(rlat_dist, rlat)
    rlat *= rlat_span
    rlat += rlat_lbound

    ## Sample recoalescence edge ##
    cedge = _sample_cedge(rng, arg, rlat, nextidx, window, live_edges[e2],
                          redge, buffer)

    ## Sample recoalescence latitude ##
    clat_lbound = max(rlat, latitude(arg, dst(cedge)))
    clat_ubound = latitude(arg, src(cedge))
    clat_span = clat_ubound - clat_lbound
    clat_dist = Beta(2)
    clat = rand(rng, clat_dist)
    arg.logprob[] += logpdf(clat_dist, clat)
    clat *= clat_span
    clat += clat_lbound

    ## Add recombination event to the graph ##
    @debug "Constrained recombination event" redge cedge breakpoint rlat clat
    recombine!(arg, redge, cedge, breakpoint, rlat, clat, buffer = buffer)

    ## Compute new live edge ##
    news = src(cedge)
    newd = nv(arg)
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
        redges = unsafe_wrap(Array, redges_ptr, redges_len)

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
        cedges_data = @alloc(Edge{VertexType}, ne(arg))
        cedges_ptr = firstindex(cedges_data)
        for e ∈ edges_interval(arg, window, store, mrca(arg), clat)
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

function build!(rng, arg::Arg; winwidth = ∞, buffer = default_buffer())
    pbar = ProgressBar(total = nmarkers(arg), printing_delay = 1, unit = " markers")

    ## Unconstrained recombinations ##
    ρ = rec_rate(arg, true)
    nrecs_dist = Poisson(ρ)
    nrecs = rand(rng, nrecs_dist)
    arg.logprob[] += logpdf(nrecs_dist, nrecs)

    @no_escape buffer begin
        for _ ∈ 1:nrecs
            sample_recombination_unconstrained!(rng, arg, winwidth, buffer)
        end
    end

    ## Constrained recombinations ##
    @no_escape buffer begin
        _mutation_edges = [
            sizehint!(Vector{Edge{VertexType}}(undef, 0), nv(arg) ÷ 2) for
            _ ∈ 1:nmarkers(arg)
        ]
        mutation_edges!(_mutation_edges, arg, Ω(first(positions(arg)), ∞),
                        buffer = buffer)

        meidx = findfirst(>(1) ∘ length, _mutation_edges)
        nextidx = 0

        while !isnothing(meidx)
            nextidx += meidx
            live_edges = _mutation_edges[meidx]
            nbp = length(live_edges) - 1

            bp_lbound = isone(nextidx) ?
                zero(eltype(positions(arg))) : idxtopos(arg, nextidx - 1)
            bp_ubound = idxtopos(arg, nextidx)

            ## Constrained recombinations ##
            nbp = length(live_edges) - 1
            if nbp > 0
                bp_dist = Uniform(bp_lbound, bp_ubound)
                ## TODO: call to `sort` is not optimized for some reason.
                breakpoints = (sort ∘ rand)(rng, bp_dist, nbp)
                arg.logprob[] -= nbp * log(bp_ubound - bp_lbound)

                for breakpoint ∈ breakpoints
                    sample_recombination_constrained!(rng, arg, breakpoint,
                                                      winwidth, live_edges,
                                                      buffer = buffer)
                end
            end

            update(pbar, meidx)

            nextidx >= nmarkers(arg) && break
            mutation_edges!(_mutation_edges, arg, Ω(idxtopos(arg, nextidx + 1), ∞),
                            buffer = buffer)
            meidx = findfirst(>(1) ∘ length, _mutation_edges)
            isnothing(meidx) && break
        end
        update(pbar, nmarkers(arg) - nextidx)
    end

    arg
end

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

function validate(arg::Arg)
    flag = true
    n = nleaves(arg)

    ## General properties ##
    if nmutations(arg) != nmarkers(arg)
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

        if ancestral_intervals(arg, v) != Set([Ω(0, ∞)])
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

        ai_child = ancestral_intervals(arg, Edge(v => child(arg, v)))
        ai_right = ancestral_intervals(arg, Edge(rightdad(arg, v) => v))
        leftparent = dads(arg, v)[findfirst(!=(rightdad(arg, v)), dads(arg, v))]
        ai_left = ancestral_intervals(arg, Edge(leftparent => v))
        if ai_child != ai_left ∪ ai_right
            msg = "Recombination vertex whose children edge's ancestral" *
                " interval is not equal to the union of it parental edges's" *
                " ancestral intervals."
            @info msg v ai_child ai_left ai_right
            flag = false
        end
    end

    ## Internal vertices ##
    for v ∈ ivertices(arg)
        ref = mapreduce(&, children(arg, v)) do _child
            sequence(arg, _child) | ~ancestral_mask(arg, Edge(v => _child))
        end
        ref &= ancestral_mask(arg, v)

        if sequence(arg, v) != ref
            @info "Inconsistent sequence" v
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
