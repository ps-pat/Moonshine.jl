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
    visited::UnsafeArray{Bool, 1}
    min_latitude::Float64
    breakpoint::Float64
    nextidx::Int
end

function EdgesIntervalRec(arg, ωs, store::AbstractArray, visited, breakpoint, nextidx,
                          root = mrca(arg), min_latitude = zero(Float64))
    ##TODO: manage `visited` and `funbuffer` manually.
    eibuffer = CheapStack(store)
    fill!(visited, false)

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
        visited = @alloc(Bool, nrecombinations(arg))
        @inbounds for e ∈ EdgesIntervalRec(arg, window, store, visited,
                                             nextpos, nextidx, src(live_edge), lat)
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
    Δlat_dist = Beta(2, 2)
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
    es_idx = MVector{2, Int}(undef)
    sample!(rng, eachindex(live_edges), es_idx; replace = false)
    arg.logprob[] += log(2) - log(n) - log(n - 1)

    sample_derived_recombination!(rng, arg, first(es_idx), last(es_idx),
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
