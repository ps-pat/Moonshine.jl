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

    η
end

function update_upstream!(arg, v, stack; buffer = default_buffer())
    @no_escape buffer begin
        push!(stack, v)

        mask = @alloc(UInt, div(nmarkers(arg), blocksize(Sequence), RoundUp))

        while !isempty(stack)
            v = pop!(stack)

            if isrecombination(arg, v)
                sequences(arg)[v] = sequence(arg, child(arg, v))

                ## If `v` is a recombination vertex, the ancestral interval of
                ## one of its parental edge is the intersection of the
                ## appropriate interval associated with the breakpoint with ωv.

                @inbounds for dad ∈ dads(arg, v)
                    e = Edge(dad => v)

                    ancestral_intervals!(ancestral_intervals(arg, e), arg, v)
                    intersect!(ancestral_intervals(arg, e), recombination_mask(arg, e), buffer = buffer)

                    push!(stack, src(e))
                end
            else
                ## Update sequence of `v` ##
                h = sequence(arg, v)
                oldhash = cheap_hash(h)
                h.data.chunks .⊻= .~h.data.chunks
                _compute_sequence!(arg, v, mask)
                iszero(dad(arg, v)) && continue

                ## If `v` is a coalescence vertex the ancestral interval of its
                ## parental edges is simply ωv.
                e = Edge(dad(arg, v) => v)

                oldhash = hash(oldhash, (hash ∘ ancestral_intervals)(arg, e))

                ## Update ancestral interval of e
                ancestral_intervals!(ancestral_intervals(arg, e), arg, v)

                ## Update stack
                newhash = hash((cheap_hash ∘ sequence)(arg, dst(e)),
                               (hash ∘ ancestral_intervals)(arg ,e))
                if (oldhash) != newhash
                    push!(stack, src(e))
                end
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
    $(FUNCTIONNAME)(arg, redge, cedge, breakpoint, rlat, clat[, stack]; buffer = default_buffer())

Add a recombination event to an ARG.

# Methods
$(METHODLIST)
"""
function recombine! end

function recombine!(arg, redge, cedge, breakpoint, rlat, clat, stack;
                    buffer = default_buffer())
    ## Adjust recombination masks
    if isrecombination(arg, dst(redge)) && src(redge) < otherdad(arg, redge)
        idx = 2_recidx(arg, dst(redge)) - 1
        arg.recombination_mask[idx], arg.recombination_mask[idx + 1] =
            arg.recombination_mask[idx + 1], arg.recombination_mask[idx]
    end

    if isrecombination(arg, dst(cedge)) && src(cedge) < otherdad(arg, cedge)
        idx = 2_recidx(arg, dst(cedge)) - 1
        arg.recombination_mask[idx], arg.recombination_mask[idx + 1] =
            arg.recombination_mask[idx + 1], arg.recombination_mask[idx]
    end

    ## Add recombination and recoalescence vertices to arg
    rvertex, cvertex = VertexType.(nv(arg) .+ (1, 2))

    ## Store latitudes and haplotypes associated with new vertices
    push!(latitudes(arg), rlat, clat)
    push!(sequences(arg), sequence(arg, dst(redge)), Sequence(trues(nmarkers(arg))))

    ## Ancestral intervals of `rvertex`'s parental edges
    ωr = ancestral_intervals(arg, redge)
    ωr_left, ωr_right = copy(ωr), copy(ωr)
    intersect!(ωr_left, Ω(0, breakpoint))
    intersect!(ωr_right, Ω(breakpoint, ∞))

    arg.ancestral_intervals[Edge(rvertex, dst(redge))] = ωr
    arg.ancestral_intervals[Edge(src(redge), rvertex)] = ωr_left
    delete!(arg.ancestral_intervals, redge)

    ## Ancestral intervals of `cvertex`'s adjacent edges
    root_recombination = !has_edge(arg, cedge)
    if root_recombination
        ωc = AIsType([Ω(0, ∞)])
        arg.mrca[] = cvertex
    else
        ωc = ancestral_intervals(arg, cedge)
        delete!(arg.ancestral_intervals, cedge)
        let ωc_new = union(ωc, ωr_right)
            arg.ancestral_intervals[Edge(src(cedge), cvertex)] = ωc_new
        end
    end

    arg.ancestral_intervals[Edge(cvertex, rvertex)] = ωr_right
    arg.ancestral_intervals[Edge(cvertex, dst(cedge))] = ωc

    ## Update topology
    add_rr_event!(graph(arg), redge, cedge)

    ## Compute sequence of recoalescence vertex ##
    @no_escape buffer begin
        mask = @alloc(UInt, div(nmarkers(arg), blocksize(Sequence), RoundUp))
        _compute_sequence!(arg, cvertex, mask)
    end

    ## Update sequences and ancetral intervals ##
    update_upstream!(arg, src(redge), stack, buffer = buffer)
    update_upstream!(arg, src(cedge), stack, buffer = buffer)

    push!(arg.recombination_mask, AIsType([Ω(0, breakpoint)]))
    push!(arg.recombination_mask, AIsType([Ω(breakpoint, ∞)]))
    arg
end

function recombine!(arg, redge, cedge, breakpoint, rlat, clat; buffer = default_buffer())
    @no_escape buffer begin
        store = @alloc(VertexType, nleaves(arg) + nrecombinations(arg))
        stack = CheapStack(store)

        recombine!(arg, redge, cedge, breakpoint, rlat, clat, stack, buffer = buffer)
    end
end

"""
    extend_recombination!(arg, edge, otherdad, breakpoint; buffer = default_buffer())

Extend ancestral interval of a recombination vertex upstream edges:
* `edge` new interval is [`breakpoint`, ∞);
* the other branch's interval is intersected with [0, `breakpoint`).

This is mainly intended to be used within
[`sample_recombination_constrained!`](@ref).
"""
function gene_conversion!(arg, edge, otherdad, breakpoint, vstack;
    buffer = default_buffer())
    d = dst(edge)

    union!(recombination_mask(arg, edge), Ω(breakpoint, ∞))

    edge = Edge(otherdad => d)
    intersect!(recombination_mask(arg, edge), Ω(0, breakpoint))

    update_upstream!(arg, d, vstack, buffer = buffer)
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
        ridx = _recidx(arg, s)
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
                breakpoint ∈ recombination_mask(arg, Edge(s => d)) || continue
            end

            push!(buffer, newe)
        end
    end

    e, state + 1
end

#          +----------------------------------------------------------+
#          |                Constrained Recombinations                |
#          +----------------------------------------------------------+

function _sample_cedge(rng, arg, rlat, idx, window, coalroot, estack, buffer)
    @no_escape buffer begin
        cedges = @alloc(Edge{VertexType}, ne(arg))
        ws_ptr = convert(Ptr{Float64}, @alloc_ptr(ne(arg) * sizeof(Float64)))

        ncedges = 0
        visited = @alloc(Bool, nrecombinations(arg))
        for e ∈ EdgesIntervalRec(arg, window, estack, visited,
                                 idxtopos(arg, idx), idx, src(coalroot), rlat)
            ncedges += 1
            cedges[ncedges] = e
            w = latitude(arg, src(e)) - max(rlat, latitude(arg, dst(e)))
            unsafe_store!(ws_ptr, w, ncedges)
        end

        ws = ProbabilityWeights(UnsafeArray{Float64, 1}(ws_ptr, (ncedges,)))
        cedge_idx = sample(rng, ws)
        w = ws[cedge_idx]
        add_logdensity!(arg, log(w) - (log ∘ sum)(ws))
        cedge = cedges[cedge_idx]

        clat_dist = Beta(2, 2)
        clat = rand(rng, clat_dist)
        clat *= w
        clat = latitude(arg, src(cedge)) - clat
    end

    cedge, clat
end

struct EdgesIntervalArgCoal{T, I, E, P} <: AbstractEIterTD
    "Genealogy to iterate over"
    genealogy::T
    "Interval to consider"
    ωs::I
    "Edges buffer"
    stack::CheapStack{E}
    "True is associated recombination vertex has been visited previously"
    visited::UnsafeArray{Bool, 1}
    "Parameters for the block predicate"
    bp_pars::P
end

EdgesIntervalArgCoal(arg, ωs, stack, visited, bp_pars, root) =
    EIterTD(EdgesIntervalArgCoal, arg, ωs, stack, visited, bp_pars, root)

function block_predicate(iter::EdgesIntervalArgCoal, e)
    arg, nextidx, fedge, min_latitude = iter.bp_pars
    sequence(arg, dst(e))[nextidx] && return false
    latitude(iter.genealogy, src(e)) >= min_latitude
    e == fedge && return false
    true
end

function _sample_clat(rng, arg, minlat, fedge, nextidx, stack;
                      buffer = default_buffer())
    local ret
    @inbounds @no_escape buffer begin
        ## Compute quadrature latitudes & weights
        ws = @alloc(Float64, clat_gridsize)
        clats = @alloc(Float64, clat_gridsize)

        rg = logrange(minlat, tmrca(arg), clat_gridsize + 1)
        lat_prev, lat_it = Iterators.peel(rg)
        for (k, lat) ∈ enumerate(lat_it)
            clats[k] = lat_prev # "left" endpoint accounts for forbidden edges
            ws[k] = lat - lat_prev
            lat_prev = lat
        end

        ## Evaluate intensity function
        visited = @alloc(Bool, nrecombinations(arg))
        edges_iterator = EdgesIntervalArgCoal(arg, idxtopos(arg, nextidx), stack, visited,
                                              (arg, nextidx, fedge, first(clats)),
                                              mrca(arg))
        λts = @alloc(Int, clat_gridsize)
        nlive!(λts, arg, clats, edges_iterator)

        ## Sample latitude and pull back to ARG scale
        pp = Exponential()
        lat1 = rand(rng, pp)
        add_logdensity!(arg, pp, lat1)
        Λ = 0
        for (w, λt, rightlat) ∈ zip(ws, λts, rg[2:end])
            Λ += w * λt
            lat1 > Λ && continue

            ret = rightlat + (lat1 - Λ) / λt
            break
        end

        if lat1 > Λ
            ret = tmrca(arg) + lat1 - Λ
        end
    end

    ret
end

"""
    $(FUNCTIONNAME)(arg, e, nextidx)

Compute the minimum latitude for a recombination event rooted on `e`.

--*Internal*--
"""
function rlat_min end

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

function _update_live_edges!(live_edges, arg, newd, nextidx)
    nextpos = idxtopos(arg, nextidx)

    lastd = newd
    e = Edge((first ∘ dads)(arg, newd, nextpos) => newd)
    @inbounds while !ismutation_edge(arg, e, nextidx)
        for (k, live_edge) ∈ enumerate(live_edges)
            src(live_edge) == src(e) || continue
            popat!(live_edges, k)
            break
        end

        lastd = newd
        newd = src(e)
        e = Edge((first ∘ dads)(arg, newd, nextpos) => newd)
    end

    e ∈ live_edges || push!(live_edges, e)
    Edge(dst(e) => lastd)
end

function _sample_cedge_wild(rng, arg, lat, fedge, nextidx, redge::T, stack;
                            buffer = default_buffer()) where T
    nextpos = idxtopos(arg, nextidx)
    len = 0

    @no_escape buffer begin
        cedges_ptr = convert(Ptr{T}, @alloc_ptr(ne(arg) * sizeof(T)))

        visited = @alloc(Bool, nrecombinations(arg))
        edges_iterator = EdgesIntervalArgCoal(arg, nextpos, stack, visited,
                                              (arg, nextidx, fedge, lat), mrca(arg))

        for e ∈ edges_iterator
            latitude(arg, dst(e)) <= lat <= latitude(arg, src(e)) || continue

            len += 1
            unsafe_store!(cedges_ptr, e, len)
        end

        cedges = UnsafeArray{T, 1}(cedges_ptr, (len,))
        e = sample(rng, cedges)
    end

    add_logdensity!(arg, -log(len))
    e
end

function _breakpoint_bounds(arg, rightidx, redge, cedge, buffer)
    hr, hc = sequence(arg, dst(redge)).data.chunks, sequence(arg, dst(cedge)).data.chunks
    nchunks = length(hr)
    bs = blocksize(Sequence)

    @inbounds @no_escape buffer begin
        m = @alloc(UInt, nchunks)
        if has_edge(arg, cedge)
            ancestral_mask!(m, arg, cedge)
        else
            fill!(m, 0xffffffffffffffff)
        end

        iidx = idxinchunk(Sequence, rightidx)
        cidx = chunkidx(Sequence, rightidx)

        c = hr[cidx] ⊻ hc[cidx]
        c &= m[cidx]
        c |= ~m[cidx]
        c &= 0xffffffffffffffff >> (bs - iidx + 1)
        Δidx_chunk = leading_zeros(c) - bs + iidx
        Δidx = Δidx_chunk

        if Δidx_chunk == iidx
            Δidx_chunk = bs
        end
        while Δidx_chunk == bs && cidx > 1
            cidx -= 1

            c = hr[cidx] ⊻ hc[cidx]
            c &= m[cidx]
            c |= ~m[cidx]

            Δidx_chunk = leading_zeros(c)
            Δidx += Δidx_chunk
        end
    end

    idxtopos(arg, rightidx - Δidx), idxtopos(arg, rightidx)
end

"""
    sample_recombination_constrained!(rng, arg, breakpoint, winwidth, live_edges; buffer = default_buffer())

Sample a recombination event constrained so as to reduce the number of
mutations for a given marker by one.
"""
function sample_recombination_constrained!(
    rng, arg, nextidx, winwidth,
    live_edges, estack, vstack;
    buffer = default_buffer(),
    conversion = true)
    n = length(live_edges)
    nextpos = idxtopos(arg, nextidx)
    window = nextpos ± winwidth / 2

    @no_escape buffer begin
        ## Sample live edges ##
        es_idx = @alloc(Int, 2)
        sample!(rng, eachindex(live_edges), es_idx; replace = false)
        if first(es_idx) > last(es_idx)
            es_idx[1], es_idx[2] = es_idx[2], es_idx[1]
        end
        add_logdensity!(arg, log(2) - log(n) - log(n - 1))
        e1, e2 = @views live_edges[es_idx]
        deleteat!(live_edges, es_idx)

        ## Consider the possibility of a wild recombination
        eu = Edge(0 => 0) # Uncle edge
        newd = zero(VertexType)
        let dads_e1 = dads(arg, src(e1), nextpos),
            dads_e2 = dads(arg, src(e2), nextpos)

            if !isempty(dads_e1) && first(dads_e1) == src(e2)
                sibling_e1 = sibling(arg, dst(e1), src(e1), (nextpos,))
                eu = Edge(src(e1) => sibling_e1)
                newd = src(e2)
            elseif !isempty(dads_e2) && first(dads_e2) == src(e1)
                sibling_e2 = sibling(arg, dst(e2), src(e2), (nextpos,))
                eu = Edge(src(e2) => sibling_e2)
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
            clat = _sample_clat(rng, arg, rlat, fedge, nextidx, estack, buffer = buffer)

            if clat > tmrca(arg)
                cedge = Edge(mrca(arg) => mrca(arg))
            else
                cedge = _sample_cedge_wild(rng, arg, clat, fedge, nextidx, redge,
                                           estack, buffer = buffer)
            end
        else
            recroot, coalroot = isone(recroot_idx) ? (e1, e2) : (e2, e1)
            redge = _find_actual_edge(arg, recroot, nextidx, rlat)
            cedge, clat = _sample_cedge(rng, arg, rlat, nextidx, window, coalroot, estack, buffer)
            newd = nv(arg) + VertexType(2)
        end

        ## Sample breakpoint ##
        breakpoint_bounds = _breakpoint_bounds(arg, nextidx, redge, cedge, buffer)
        breakpoint_dist = Uniform(breakpoint_bounds...)
        breakpoint = rand(rng, breakpoint_dist)
        add_logdensity!(arg, breakpoint_dist, breakpoint)

        if conversion && src(cedge) ∈ dads(arg, dst(redge)) && isrecombination(arg, dst(redge))
            @debug "Gene conversion event" redge cedge breakpoint wildrec
            gene_conversion!(arg, Edge(src(cedge) => dst(redge)), src(redge),
                             breakpoint, vstack, buffer = buffer)
            if !wildrec
                newd = src(cedge)
            end
        else
            @debug "Constrained recombination event" redge cedge breakpoint rlat clat
            recombine!(arg, redge, cedge, breakpoint, rlat, clat, vstack, buffer = buffer)
        end

        _update_live_edges!(live_edges, arg, newd, nextidx)
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

function build!(rng, arg::Arg;
    winwidth = ∞,
    buffer = default_buffer(),
    noprogress = false,
    conversion = true)
    progenabled = nleaves(arg) * nmarkers(arg) >= 10000000
    prog = Progress(nmarkers(arg), enabled = progenabled && !noprogress)

    @no_escape buffer begin
        estore = @alloc(Edge{VertexType}, max(100, ne(arg)))
        vstore = @alloc(VertexType, max(100, convert(Int, nv(arg))))
        estack = CheapStack(estore)
        vstack = CheapStack(vstore)

        ## Constrained recombinations ##
        mutation_edges_buffer = ntuple(_ -> Edge{VertexType}[], 8mmn_chunksize)
        nextidx, live_edges = next_inconsistent_idx(arg, 1, estack,
                                                    mutations_edges = mutation_edges_buffer,
                                                    buffer = buffer)

        while !iszero(nextidx)
            update!(prog, nextidx)

            while !(isone ∘ length)(live_edges)
                sample_recombination_constrained!(
                    rng, arg, nextidx,
                    winwidth, live_edges,
                    estack, vstack,
                    buffer = buffer,
                    conversion = conversion)
            end

            nextidx, live_edges = next_inconsistent_idx(arg, nextidx + 1, estack,
                                                        mutations_edges = mutation_edges_buffer,
                                                        buffer = buffer)
        end
    end

    finish!(prog)
    arg
end
