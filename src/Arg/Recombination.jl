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
    rvertex, cvertex = VertexType.(nv(arg) .+ (1, 2))
    add_vertices!(
        arg,
        (sequence(arg, dst(redge)), Sequence(trues(nmarkers(arg)))),
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
    root_recombination = !has_edge(arg, cedge)
    if root_recombination
        ωc = AIsType([Ω(0, ∞)])
        arg.mrca[] = cvertex
    else
        ωc = ancestral_intervals(arg, cedge)
        rem_edge!(arg, cedge)
        let ωc_new = union(ωc, ωr_right)
            add_edge!(arg, Edge(src(cedge), cvertex), ωc_new)
        end
    end

    add_edge!(arg, Edge(cvertex, rvertex), ωr_right)
    add_edge!(arg, Edge(cvertex, dst(cedge)), ωc)

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

function _sample_cedge(rng, arg, rlat::T, possible_cedges, nextidx;
                       buffer = default_buffer()) where T
    cedge = zero(VertexType)
    clat = zero(T)

    @no_escape buffer begin
        ## Compute minimum latitudes & domain bounds
        domain_lbound = typemax(T)
        domain_ubound = typemin(T)
        min_latitudes = @alloc(T, length(possible_cedges))
        @inbounds for (k, e) ∈ enumerate(possible_cedges)
            max_latitude = latitude(arg, src(e))
            if domain_ubound < max_latitude
                domain_ubound = max_latitude
            end

            min_latitudes[k] = max(rlat, rlat_min(arg, e, nextidx))
            if domain_lbound > min_latitudes[k]
                domain_lbound = min_latitudes[k]
            end
        end

        ## Sample latitude
        ws = @alloc(Int, length(possible_cedges))
        pp_dist = Uniform(domain_lbound, domain_ubound)
        λbound = length(possible_cedges)
        clat = (zero ∘ eltype)(pp_dist)
        λt = zero(Int)
        accept = false
        @inbounds while !accept
            clat = rand(rng, pp_dist)
            add_logdensity!(arg, pp_dist, clat)

            ## `count` causes weird type instability for some reason
            λt = zero(Int)
            @simd ivdep for k ∈ eachindex(possible_cedges)
                λt += min_latitudes[k] <= clat <= latitude(arg, src(possible_cedges[k]))
                ws[k] = λt
            end

            accept_dist = Bernoulli(λt / λbound)
            accept = rand(rng, accept_dist)
            add_logdensity!(arg, accept_dist, accept)
        end

        cedge_idx = zero(Int)
        let k = rand(rng, 1:λt)
            cedge_idx = findfirst(==(k), ws)
        end
        add_logdensity!(arg, -log(λt))
        cedge = _find_actual_edge(arg, possible_cedges[cedge_idx], nextidx, clat)
    end

    cedge, clat
end

function _sample_clat(rng, arg, minlat, fedge, nextidx, stack;
                      buffer = default_buffer())
    nextpos = idxtopos(arg, nextidx)

    block_predicates = [
        e -> sequence(arg, dst(e))[nextidx],
        e -> e == fedge
    ]

    ret = zero(minlat)
    @no_escape buffer begin
        n = 100
        λts = @alloc(Int, n)
        clat_shifts = @alloc(Float64, n)
        c = (one ∘ eltype)(λts)
        pp = Exponential()

        while iszero(ret)
            rand!(rng, pp, clat_shifts)

            nlive!(λts, arg, clat_shifts .+ minlat, nextpos, stack,
                   block_predicates = block_predicates, buffer = buffer)

            for (clat_shift, λt) ∈ zip(clat_shifts, λts)
                add_logdensity!(arg, pp, clat_shift)

                ## Acceptation
                iszero(λt) && continue

                r = exp((1 - λt) * clat_shift)
                accept_dist = Bernoulli(r / c)
                accept = rand(rng, accept_dist)
                add_logdensity!(arg, accept_dist, accept)

                if accept
                    ret = clat_shift + minlat
                    break
                end

                c = max(c, r)
            end
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
    e = Edge(dad(arg, newd) => newd) ## `newd` is a coalescence vertex
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
function sample_recombination_constrained!(rng, arg, nextidx, winwidth,
                                           live_edges, estack, vstack;
                                           buffer = default_buffer())
    n = length(live_edges)
    nextpos = idxtopos(arg, nextidx)
    window = nextpos ± winwidth / 2

    local breakpoint
    let breakpoint_bounds = (idxtopos(arg, nextidx - 1), nextpos),
        breakpoint_dist = Uniform(breakpoint_bounds...)
        breakpoint = rand(rng, breakpoint_dist)
        add_logdensity!(arg, breakpoint_dist, breakpoint)
    end

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
        let dads_e1 = dads(arg, src(e1), breakpoint),
            dads_e2 = dads(arg, src(e2), breakpoint)

            if !isempty(dads_e1) && first(dads_e1) == src(e2)
                sibling_e1 = sibling(arg, dst(e1), src(e1), (breakpoint,))
                if !iszero(sibling_e1)
                    eu = Edge(src(e1) => sibling_e1)
                    newd = src(e2)
                end
            elseif !isempty(dads_e2) && first(dads_e2) == src(e1)
                sibling_e2 = sibling(arg, dst(e2), src(e2), (breakpoint,))
                if !iszero(sibling_e2)
                    eu = Edge(src(e2) => sibling_e2)
                    newd = src(e1)
                end
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
            # clat = _sample_clat(rng, arg, rlat, possible_cedges, nextidx, buffer = buffer)
            cedge, clat = _sample_cedge(rng, arg, rlat, possible_cedges, nextidx, buffer = buffer)
            newd = nv(arg) + VertexType(2)
        end

        @debug "Constrained recombination event" redge cedge breakpoint rlat clat
        recombine!(arg, redge, cedge, breakpoint, rlat, clat, vstack, buffer = buffer)
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

            while !(isone ∘ length)(live_edges)
                sample_recombination_constrained!(rng, arg, nextidx,
                                                  winwidth, live_edges,
                                                  estack, vstack,
                                                  buffer = buffer)
            end

            nextidx, live_edges = next_inconsistent_idx(arg, nextidx + 1, estack,
                                                        mutations_edges = mutation_edges_buffer,
                                                        buffer = buffer)
        end
    end

    finish!(prog)
    arg
end
