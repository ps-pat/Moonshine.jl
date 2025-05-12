#          +----------------------------------------------------------+
#          |                        ARG Update                        |
#          +----------------------------------------------------------+

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
    newhash = hash((cheap_hash ∘ sequence)(arg, dst(e)),
                   (hash ∘ ancestral_intervals)(arg ,e))
    if (oldhash) != newhash
        push!(vstack, src(e))
    end

    ancestral_intervals(arg, e)
end

function update_upstream!(arg, v, stack; buffer = default_buffer())
    @no_escape buffer begin
        push!(stack, v)

        mask = @alloc(UInt, div(nmarkers(arg), blocksize(Sequence), RoundUp))
        ωsv = (AIsType(), AIsType())

        while !isempty(stack)
            v = pop!(stack)

            ## Update sequence of `v` ##
            h = sequence(arg, v)
            oldhash = cheap_hash(h)
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
                    _update_ai!(stack, arg, e, ωsv[k], oldhash)
                end
            else
                ## If `v` is a coalescence vertex the ancestral interval of its
                ## parental edges is simply ωv.
                e = Edge(dad(arg, v) => v)

                _update_ai!(stack, arg, e, first(ωsv), oldhash)
            end
        end
    end

    arg
end

#          +----------------------------------------------------------+
#          |                      Recombination                       |
#          +----------------------------------------------------------+

export recombine!
"""
    recombine!(arg, redge, cedge, breakpoint, rlat, clat)

Add a recombination event to an ARG.
"""
function recombine!(arg, redge, cedge, breakpoint, rlat, clat, stack;
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
    root_recombination = !rem_edge!(arg, cedge)
    if root_recombination
        ωc = AIsType([Ω(0, ∞)])
        arg.mrca[] = cvertex
    else
        ωc = ancestral_intervals(arg, cedge)
        let ωc_new = union(ωc, ωr_right)
            add_edge!(arg, Edge(src(cedge), cvertex), ωc_new)
        end
    end

    add_edge!(arg, Edge(cvertex, rvertex), ωr_right)
    add_edge!(arg, Edge(cvertex, dst(cedge)), ωc)

    ## Compute sequence of new vertices ##
    @no_escape buffer begin
        mask = @alloc(UInt, div(nmarkers(arg), blocksize(Sequence), RoundUp))
        _compute_sequence!(arg, rvertex, mask)
        _compute_sequence!(arg, cvertex, mask)
    end

    ## Update sequences and ancetral intervals ##
    update_upstream!(arg, src(redge), stack, buffer = buffer)
    update_upstream!(arg, src(cedge), stack, buffer = buffer)

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

#          +----------------------------------------------------------+
#          |                        Utilities                         |
#          +----------------------------------------------------------+

struct EdgesIntervalRec{I}
    genealogy::Arg
    ωs::I
    buffer::CheapStack{Edge{VertexType}}
    visited::UnsafeArray{Bool, 1}
    min_latitude::Float64
    breakpoint::Float64
    nextidx::Int
end

function EdgesIntervalRec(arg, ωs, stack, visited, breakpoint, nextidx,
                          root = mrca(arg), min_latitude = zero(Float64))
    ##TODO: manage `visited` and `funbuffer` manually.
    fill!(visited, false)

    for d ∈ children(arg, root, ωs)
        e = Edge(root => d)
        (breakpoint ∈ ancestral_intervals(arg, e) && !sequence(arg, dst(e))[nextidx]) && continue
        push!(stack, e)
    end

    EdgesIntervalRec(arg, ωs, stack, visited, min_latitude, breakpoint, nextidx)
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

#          +----------------------------------------------------------+
#          |                Constrained Recombinations                |
#          +----------------------------------------------------------+

function _sample_clat(rng, arg, rlat, possible_cedges)
    domain_lbound = max(rlat, minimum(e -> latitude(arg, dst(e)), possible_cedges))
    domain_ubound = maximum(e -> latitude(arg, src(e)), possible_cedges)
    pp_dist = Uniform(domain_lbound, domain_ubound)

    λbound = length(possible_cedges)
    t = (zero ∘ eltype)(pp_dist)
    accept = false
    @inbounds while !accept
        t = rand(rng, pp_dist)
        add_logdensity!(arg, pp_dist, t)

        ## `count` causes weird type instability for some reason
        λt = zero(Int)
        @simd ivdep for k ∈ eachindex(possible_cedges)
            e = possible_cedges[k]
            λt += latitude(arg, dst(e)) <= t <= latitude(arg, src(e))
        end

        accept_dist = Bernoulli(λt / λbound)
        accept = rand(rng, accept_dist)
        add_logdensity!(arg, accept_dist, accept)
    end

    t
end

function _sample_clat(rng, arg, minlat, redge, fedge, nextidx, stack;
                      buffer = default_buffer())
    nextpos = idxtopos(arg, nextidx)

    ## Homogeneous PP
    λ = nv(arg) - nrecombinations(arg)
    pp = Exponential(inv(λ))

    block_predicates = [
        e -> sequence(arg, dst(e))[nextidx],
        e -> e == fedge
    ]

    ret = zero(minlat)
    @no_escape buffer begin
        n = nleaves(arg)
        λts = @alloc(Int, n)
        clats = @alloc(Float64, n)

        while iszero(ret)
            Δclats = rand(rng, pp, n)
            cumsum!(clats, Δclats)
            clats .+= minlat

            nlive!(λts, arg, clats, nextpos, stack, block_predicates = block_predicates, buffer = buffer)

            for (k, (clat, λt)) ∈ (enumerate ∘ zip)(clats, λts)
                add_logdensity!(arg, pp, Δclats[k])
                Δt = clats[k] - (isone(k) ? minlat : clats[k - 1])

                ## Acceptation
                accept_dist = Bernoulli(λt / λ)
                accept = rand(rng, accept_dist)
                add_logdensity!(arg, accept_dist, accept)

                if accept
                    ret = clat
                    break
                end
            end

            minlat = last(clats)
        end
    end

    ret
end

"""
    rlat_interval(arg, e, nextidx, ubound)

Compute the valid recombination interval for an edge at a given position.
"""
function rlat_min end,
function rlat_min! end

function rlat_min(arg, e, nextidx)
    nextpos = idxtopos(arg, nextidx)
    s, d = src(e), dst(e)
    refstate = sequence(arg, d)[nextidx]

    @inbounds while true
        _children = children(arg, d, nextpos)
        (isone ∘ length)(_children) || break
        c = first(_children)
        sequence(arg, c)[nextidx] == refstate || break

        s = d
        d = c
    end

    latitude(arg, d)
end

function _find_actual_edge(arg, e, nextidx, lat)
    nextpos = idxtopos(arg, nextidx)

    while true
        latitude(arg, dst(e)) <= lat <= latitude(arg, src(e)) && return e
        e = Edge(dst(e) => (first ∘ children)(arg, dst(e), nextpos))
    end
end

function _update_live_edges_wild!(live_edges, arg, newd, nextidx, breakpoint, es_idx;
                                  buffer = default_buffer())
    e = _update_live_edges_derived!(live_edges, arg, newd, nextidx, breakpoint, es_idx,
                                    buffer = buffer)
    e = Edge(src(e) => sibling(arg, dst(e)))

    eidx = findfirst(==(e), live_edges)
    isnothing(eidx) || popat!(live_edges, eidx)

    live_edges
end

function _update_live_edges_derived!(live_edges, arg, newd, nextidx, breakpoint, es_idx;
                                     buffer = default_buffer())
    popat!(live_edges, last(es_idx))
    popat!(live_edges, first(es_idx))

    lastd = newd
    e = Edge((first ∘ dads)(arg, newd, breakpoint) => newd)
    @no_escape buffer begin
        ptr = convert(Ptr{VertexType}, @alloc_ptr(2 * sizeof(VertexType)))

        @inbounds while !ismutation_edge(arg, e, nextidx)
            let eidx = findfirst(==(e), live_edges)
                if !isnothing(eidx)
                    popat!(live_edges, eidx)
                else
                    ss = siblings!(ptr, arg, dst(e), breakpoint)
                    if !isempty(ss)
                        eidx = findfirst(==(Edge(src(e) => first(ss))), live_edges)
                        isnothing(eidx) || popat!(live_edges, eidx)
                    end
                end
            end

            lastd = newd
            newd = src(e)
            e = Edge((first ∘ dads)(arg, newd, breakpoint) => newd)
        end
    end

    e ∈ live_edges || push!(live_edges, e)
    Edge(dst(e) => lastd)
end

function _sample_cedge_wild(rng, arg, lat, fedge, nextidx, redge::T, stack;
                            buffer = default_buffer()) where T
    nextpos = idxtopos(arg, nextidx)

    block_predicates = [
        e -> sequence(arg, dst(e))[nextidx],
        e -> e == fedge
    ]

    @no_escape buffer begin
        cedges_ptr = convert(Ptr{T}, @alloc_ptr(ne(arg) * sizeof(T)))
        len = 0

        visited = @alloc(Bool, nrecombinations(arg))
        @inbounds for e ∈ edges_interval(arg, nextpos, stack, visited, mrca(arg), lat,
                                           block_predicates = block_predicates)
            latitude(arg, dst(e)) <= lat <= latitude(arg, src(e)) || continue

            len += 1
            unsafe_store!(cedges_ptr, e, len)
        end

        cedges = UnsafeArray{T, 1}(cedges_ptr, (len,))
        arg.logdensity[] -= log(len)
        sample(rng, cedges)
    end
end

"""
    sample_recombination_constrained!(rng, arg, breakpoint, winwidth, live_edges; buffer = default_buffer())

Sample a recombination event constrained so as to reduce the number of
mutations for a given marker by one.
"""
function sample_recombination_constrained!(rng, arg, breakpoint, winwidth,
                                           live_edges, estack, vstack;
                                           buffer = default_buffer())
    n = length(live_edges)
    nextidx = postoidx(arg, breakpoint)
    nextpos = idxtopos(arg, nextidx)
    window = breakpoint ± winwidth / 2

    @no_escape buffer begin
        ## Sample live edges ##
        es_idx = @alloc(Int, 2)
        sample!(rng, eachindex(live_edges), es_idx; replace = false)
        if first(es_idx) > last(es_idx)
            es_idx[1], es_idx[2] = es_idx[2], es_idx[1]
        end
        add_logdensity!(arg, log(2) - log(n) - log(n - 1))
        e1, e2 = @views live_edges[es_idx]

        ## Consider the possibility of a wild recombination
        eu = Edge(0 => 0) # Uncle edge
        newd = zero(VertexType)
        let dads_e1 = dads(arg, src(e1), breakpoint),
            dads_e2 = dads(arg, src(e2), breakpoint)

            if !isempty(dads_e1) && first(dads_e1) == src(e2)
                siblings_e1 = siblings(arg, dst(e1), breakpoint)
                if !isempty(siblings_e1)
                    eu = Edge(src(e1) => first(siblings_e1))
                end
                newd = src(e2)
            elseif !isempty(dads_e2) && first(dads_e2) == src(e1)
                siblings_e2 = siblings(arg, dst(e2), breakpoint)
                if !isempty(siblings_e2)
                    eu = Edge(src(e2) => first(siblings_e2))
                end
                newd = src(e1)
            end
        end

        ## Sample recombination location ##
        es_min_eu = (iszero ∘ dst)(eu) ? Inf : rlat_min(arg, eu, nextidx)
        es_mins = (rlat_min(arg, e1, nextidx), rlat_min(arg, e2, nextidx), es_min_eu)
        es_max = min(latitude(arg, src(e1)), latitude(arg, src(e2)))
        rints_width = @SVector [
            max(0, es_max - es_mins[1]),
            max(0, es_max - es_mins[2]),
            isinf(es_mins[3]) ? zero(Float64) : latitude(arg, src(eu)) - es_mins[3]]

        ## We proceed in two steps to avoid sampling a location too close to
        ## the ends of the recombination branch, which may cause numerical
        ## instability.
        ## First, we sample a recombination edge.
        recroot_idx = one(Int)
        @inbounds let t = sum(rints_width) * rand(rng)
            if t > first(rints_width)
                t -= first(rints_width)
                recroot_idx = t <= rints_width[2] ? 2 : 3
            end
end

        add_logdensity!(arg, log(rints_width[recroot_idx]) - (log ∘ sum)(rints_width))
        wildrec = recroot_idx == 3

        ## The larger α - β is, the stronger the bias towards ancient branches
        ## (0 = no bias)
        rlat_dist = Beta(2, 2)
        rlat = rand(rng, rlat_dist)
        add_logdensity!(arg, rlat_dist, rlat)
        rlat *= rints_width[recroot_idx]
        rlat += es_mins[recroot_idx]

        if wildrec
            redge = _find_actual_edge(arg, eu, nextidx, rlat)
            fedge =  Edge((first ∘ dads)(arg, src(eu)) => src(eu))

            ## Sample recoalescence location ##
            clat = _sample_clat(rng, arg, rlat, redge, fedge, nextidx, estack, buffer = buffer)

            if clat > tmrca(arg)
                cedge = Edge(mrca(arg) => mrca(arg))
            else
                cedge = _sample_cedge_wild(rng, arg, clat, fedge, nextidx, redge,
                                           estack, buffer = buffer)
            end
        else
            recroot, coalroot = isone(recroot_idx) ? (e1, e2) : (e2, e1)
            redge = _find_actual_edge(arg, recroot, nextidx, rlat)

            ## Sample recoalescence location ##
            visited = @alloc(Bool, nrecombinations(arg))
            possible_cedges_ptr = convert(Ptr{Edge{VertexType}},
                                          @alloc_ptr(ne(arg) * sizeof(Edge{VertexType})))
            npossible_cedges = zero(Int)
            for e ∈ EdgesIntervalRec(arg, window, estack, visited, nextpos, nextidx, src(coalroot), rlat)
                npossible_cedges += 1
                unsafe_store!(possible_cedges_ptr, e, npossible_cedges)
            end
            possible_cedges = UnsafeArray{Edge{VertexType}, 1}(possible_cedges_ptr, (npossible_cedges,))
            clat = _sample_clat(rng, arg, rlat, possible_cedges)

            ## Sample recoalescence edge
            ws = @alloc(Int, npossible_cedges)
            cumweight = zero(Int)
            for (k, e) ∈ enumerate(possible_cedges)
                cumweight += latitude(arg, dst(e)) <= clat <= latitude(arg, src(e))
                ws[k] = cumweight
            end

            cedge_idx = zero(Int)
            let k = rand(rng, 1:cumweight)
                cedge_idx = findfirst(==(k), ws)
            end
            add_logdensity!(arg, -log(cumweight))
            cedge = possible_cedges[cedge_idx]
        end

        @debug "Constrained recombination event" redge cedge breakpoint rlat clat
        recombine!(arg, redge, cedge, breakpoint, rlat, clat, vstack, buffer = buffer)

        if wildrec
            _update_live_edges_wild!(live_edges, arg, newd, nextidx, breakpoint, es_idx)
        else
            newd = nv(arg)
            _update_live_edges_derived!(live_edges, arg, newd, nextidx, breakpoint, es_idx)
        end
    end


    breakpoint
end

#          +----------------------------------------------------------+
#          |               Unconstrained Recombination                |
#          +----------------------------------------------------------+

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
            add_logdensity!(arg, breakpoint_dist, breakpoint)

            ## Sample the recombination latitude ##
            rlat_dist = Beta(2)
            rlat = rand(rng, rlat_dist)
            add_logdensity!(arg, rlat_dist, rlat)
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
    add_logdensity!(arg, Δclat_dist, Δclat)
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

#          +----------------------------------------------------------+
#          |                       ARG Building                       |
#          +----------------------------------------------------------+

"""
    $(SIGNATURES)

Build an ancestral recombination graph.

# Arguments
* `winwidth` (`∞`): width of the window of positions to consider
* `noprogress` (`false`): hide progress bar
"""
function build!(rng, arg::Arg; winwidth = ∞, buffer = default_buffer(), noprogress = false)
    progenabled = nleaves(arg) * nmarkers(arg) >= 10000000
    prog = Progress(nmarkers(arg), enabled = progenabled && !noprogress)

    @no_escape buffer begin
        estore = @alloc(Edge{VertexType}, max(100, nleaves(arg)))
        vstore = @alloc(VertexType, max(100, nleaves(arg)))
        estack = CheapStack(estore)
        vstack = CheapStack(vstore)

        ## Constrained recombinations ##
        mutation_edges_buffer = ntuple(_ -> Edge{VertexType}[], 8mmn_chunksize)
        nextidx, live_edges = next_inconsistent_idx(arg, 1, estack,
                                                    mutations_edges = mutation_edges_buffer,
                                                    buffer = buffer)

        while !iszero(nextidx)
            update!(prog, nextidx)

            nbp = length(live_edges) - 1

            bp_lbound = isone(nextidx) ?
                zero(eltype(positions(arg))) : idxtopos(arg, nextidx - 1)
            bp_ubound = idxtopos(arg, nextidx)

            bp_dist = Uniform(bp_lbound, bp_ubound)
            breakpoints = (sort! ∘ rand)(rng, bp_dist, nbp)
            arg.logdensity[] -= nbp * log(bp_ubound - bp_lbound)

            for breakpoint ∈ breakpoints
                sample_recombination_constrained!(rng, arg, breakpoint,
                                                  winwidth, live_edges,
                                                  estack, vstack,
                                                  buffer = buffer)
                (isone ∘ length)(live_edges) && break
            end

            previdx = nextidx
            nextidx, live_edges = next_inconsistent_idx(arg, nextidx + 1, estack,
                                                        mutations_edges = mutation_edges_buffer,
                                                        buffer = buffer)
        end
    end

    finish!(prog)
    arg
end
