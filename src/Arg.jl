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

function mrca(arg::Arg, vs, ωs)
    μ = mrca(arg, vs)
    iszero(μ) && return μ

    while true
        @label beg
        for c ∈ children(arg, μ, ωs)
            vs ⊆ descendants(arg, c, ωs) || continue
            μ = c
            @goto beg
        end
        break
    end

    μ
end

########################
# Ancestrality Methods #
########################

function ancestral_intervals!(ωs, arg::Arg, e::Edge; wipe = true)
    if wipe
        copy!(ωs, get(() -> Set{Ω}((Ω(0, ∞),)), arg.ancestral_intervals, e))
    else
        haskey(arg.ancestral_intervals, e) || return push!(ωs, Ω(0, ∞))

        for ω ∈ arg.ancestral_intervals[e]
            push!(ωs, ω)
        end
    end

    ωs
end

## Use this `get!` method so the interval doesn't get constructed every call.
ancestral_intervals(arg::Arg, e::Edge) =
    get!(() -> Set{Ω}((Ω(0, ∞),)), arg.ancestral_intervals, e)

ancestral_intervals!(ωs, arg::Arg, s::VertexType, d; wipe = true) =
    ancestral_intervals!(ωs, arg, Edge(s, d); wipe = wipe)

ancestral_intervals(arg::Arg, s::VertexType, d) =
    ancestral_intervals(arg, Edge(s, d))

function ancestral_intervals!(ωs, arg::Arg, v::VertexType; wipe = true)
    wipe && empty!(ωs)
    isleaf(arg, v) && return push!(ωs, Ω(0, ∞))

    for child ∈ children(arg, v)
        ancestral_intervals!(ωs, arg, v, child, wipe = false)
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

Returns the parent of recombination vertex `v` ancestral for material
to the right of the breakpoint associated with `v`. If `v` is not a
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
        @inbounds for edge ∈ edges_interval(arg, ω, store)
            mutationsidx!(mutations, mask, arg, edge, firstchunk, firstidx, lastchunk,
                          ωs_buf = ωs_buf)
        end
    end

    mutations
end

function mutation_edges(arg, ω::Ω)
    ret = [Vector{Edge{VertexType}}(undef, 0) for _ ∈ 1:nmarkers(arg)]
    mutation_edges!(ret, arg, ω)
end

##################
# Recombinations #
##################

function _compute_sequence!(arg, v, mask; ωs_buf = Set{Ω}())
    η = sequence(arg, v)
    η.data.chunks .⊻= .~η.data.chunks

    @inbounds for child ∈ children(arg, v)
        ancestral_mask!(mask, arg, Edge(v, child), ωs_buf = ωs_buf)
        η.data.chunks .&= sequence(arg, child).data.chunks .| .~mask.data.chunks
    end

    ancestral_mask!(mask, arg, v, ωs_buf = ωs_buf)
    η.data.chunks .&= mask.data.chunks
    η
end

function _update_ai!(arg, e, ωs)
    if haskey(arg.ancestral_intervals, e)
        empty!(arg.ancestral_intervals[e])
    else
        arg.ancestral_intervals[e] = Set{Ω}()
    end

    while !isempty(ωs)
        push!(arg.ancestral_intervals[e], pop!(ωs))
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
            ancestral_intervals!(ωsv, arg, v)

            if isrecombination(arg, v)
                ## If `v` is a recombination vertex, the ancestral interval of
                ## one of its parental edge is the intersection of the
                ## appropriate interval associated with the breakpoint with ωv.
                bp = recbreakpoint(arg, v)
                _dads = dads(arg, v)

                for _dad ∈ _dads
                    e = Edge(_dad => v)
                    ωs_oldhash = (hash ∘ ancestral_intervals)(arg, e)

                    ωbp = _dad == rightdad(arg, v) ? Ω(bp, ∞) : Ω(0, bp)
                    ωs = ωsv ∩ ωbp
                    _update_ai!(arg, e, ωs)
                    if (hash ∘ ancestral_intervals)(arg, e) != ωs_oldhash ||
                        (hash ∘ sequence)(arg, v) != sequence_oldhash
                        push!(vstack, _dad)
                    end
                end
            else
                ## If `v` is a coalescence vertex the ancestral interval of its
                ## parental edges is simply ωv.
                _dad = dad(arg, v)
                e = Edge(_dad => v)
                ωs_oldhash = (hash ∘ ancestral_intervals)(arg, e)

                _update_ai!(arg, e, ωsv)
                if (hash ∘ ancestral_intervals)(arg, e) != ωs_oldhash ||
                    (hash ∘ sequence)(arg, v) != sequence_oldhash
                    push!(vstack, _dad)
                end
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
    if !root_recombination
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
    root_recombination || update_upstream!(arg, src(cedge), buffer = buffer)

    push!(arg.rightdads, cvertex)
    push!(arg.breakpoints, breakpoint)
    arg
end

# -- Edges Iterator ----------------------------------------------------

struct EdgesIntervalRec{I}
    genealogy::Arg
    ωs::I
    buffer::CheapStack{Edge{VertexType}}
    funbuffer::Vector{VertexType}
    visited::BitVector
    min_latitude::Float64
    breakpoint::Float64
    breakpointidx::Int
end

function EdgesIntervalRec(arg, ωs, store, breakpoint,
                          root = mrca(arg), min_latitude = zero(Float64))
    ##TODO: manage `visited` and `funbuffer` manually.
    eibuffer = CheapStack(store)
    funbuffer = Vector{VertexType}(undef, 2)
    visited = falses(nrecombinations(arg))

    breakpointidx = postoidx(arg, breakpoint)
    for d ∈ children!(funbuffer, arg, root, ωs)
        e = Edge(root => d)
        (breakpoint ∈ ancestral_intervals(arg, e) && !sequence(arg, dst(e))[breakpointidx]) && continue
        push!(eibuffer, e)
    end

    EdgesIntervalRec(arg, ωs, eibuffer, funbuffer, visited, min_latitude,
                     breakpoint, breakpointidx)
end

IteratorSize(::EdgesIntervalRec) = Base.SizeUnknown()

eltype(::EdgesIntervalRec) = Edge{VertexType}

function iterate(iter::EdgesIntervalRec, state = 1)
    buffer = iter.buffer
    isempty(buffer) && return nothing

    arg = iter.genealogy
    ωs = iter.ωs
    funbuffer = iter.funbuffer
    visited = iter.visited
    min_latitude = iter.min_latitude
    breakpoint = iter.breakpoint
    breakpointidx = iter.breakpointidx

    e = pop!(buffer)
    s = dst(e)
    if isrecombination(arg, s)
        ridx = recidx(arg, s)
        visited[ridx] && return e, state + 1
        visited[ridx] = true
    end

    resize!(funbuffer, 2)
    if latitude(arg, s) >= min_latitude
        for d ∈ children!(funbuffer, arg, s, ωs)
            newe = Edge(s => d)
            if breakpoint ∈ ancestral_intervals(arg, newe)
                sequence(arg, dst(newe))[breakpointidx] || continue
            end

            push!(buffer, newe)
        end
    end

    e, state + 1
end

function _sample_cedge(rng, arg, lat, idx, window, live_edge::T, taboo, buffer) where T
    pos = idxtopos(arg, idx)

    @no_escape buffer begin
        cedges = @alloc(T, ne(arg))
        cedges_ptr = firstindex(cedges)

        store = @alloc(T, nv(arg))
        @inbounds for e ∈ EdgesIntervalRec(arg, window, store, pos, src(live_edge), lat)
            ## Might happen if `src(live_edge) == src(taboo)`.
            src(e) == src(taboo) && continue
            ## Might happen if `dst(taboo)` is a recombination vertex.
            dst(e) == dst(taboo) && continue

            cedges[cedges_ptr] = e
            cedges_ptr += 1
        end

        ret = rand(rng, view(cedges, 1:(cedges_ptr-1)))
        ret
    end
end

function sample_recombination_constrained!(rng, arg, breakpoint, winwidth, live_edges;
                                           buffer = default_buffer())
    n = length(live_edges)

    ## Sample recombination edge ##
    @no_escape buffer begin
        ws_data = @alloc(Float64, length(live_edges))
        map!(e -> max(eps(Float64), branchlength(arg, e)), ws_data, live_edges)
        ws = ProbabilityWeights(ws_data)
        e1, e2 = sample(rng, eachindex(live_edges), ws, 2; replace = false)
    end
    arg.logprob[] += log(2) - log(n) - log(n - 1)

    ## This ensures that there is a least one edge available for recoalescence.
    if latitude(arg, dst(live_edges[e1])) > latitude(arg, dst(live_edges[e2]))
        e1, e2 = e2, e1
    end

    redge = popat!(live_edges, e1)
    if e2 > e1
        e2 -= 1
    end

    window = breakpoint ± winwidth

    ## Sample recombination latitude ##
    rlat_lbound = latitude(arg, dst(redge))
    rlat_ubound = min(latitude(arg, src(redge)), latitude(arg, src(live_edges[e2])))
    if rlat_lbound < rlat_ubound
        rlat_dist = Uniform(rlat_lbound, rlat_ubound)
        rlat = rand(rng, rlat_dist)
        arg.logprob[] += logpdf(rlat_dist, rlat)
    else
        rlat = rlat_lbound
        @debug "Recombination edge has length 0" redge
    end

    ## Sample recoalescence edge ##
    breakpoint_idx = postoidx(arg, breakpoint)
    cedge = _sample_cedge(rng, arg, rlat, breakpoint_idx,
                          window, live_edges[e2], redge, buffer)

    ## Sample recoalescence latitude ##
    clat_lbound = max(rlat, latitude(arg, dst(cedge)))
    clat_ubound = latitude(arg, src(cedge))
    if clat_lbound < clat_ubound
        ## It would be more accurate to sample from a shifted-truncated
        ## exponential distribution. However, it is highly numerically
        ## unstable to do so :(
        Δclat_dist = Uniform(clat_lbound - rlat, clat_ubound - rlat)
        Δclat = rand(rng, Δclat_dist)
        arg.logprob[] += logpdf(Δclat_dist, Δclat)
        clat = rlat + Δclat
    else
        clat = rlat + clat_lbound
        @debug "Interval available for recoalescence has length 0" redge
    end

    ## Add recombination event to the graph ##
    @debug "Recombination event" redge cedge breakpoint rlat clat
    recombine!(arg, redge, cedge, breakpoint, rlat, clat, buffer = buffer)

    ## Compute new live edge ##
    news = src(cedge)
    newd = nv(arg)
    if news != src(live_edges[e2])
        _dads = Vector{VertexType}(undef, 2)
        while sequence(arg, news)[breakpoint_idx]
            newd = news
            news = (first ∘ dads!)(_dads, arg, news, breakpoint)
        end
    end

    for (k, e) ∈ enumerate(live_edges)
        src(e) == news || continue
        live_edges[k] = Edge(news => newd)
    end

    breakpoint
end

function build!(rng, arg::Arg; winwidth = ∞, buffer = default_buffer())
    @no_escape buffer begin
        _mutation_edges = [
            sizehint!(Vector{Edge{VertexType}}(undef, 0), nv(arg) ÷ 2) for
            _ ∈ 1:nmarkers(arg)
        ]
        mutation_edges!(_mutation_edges, arg, Ω(first(positions(arg)), ∞),
                        buffer = buffer)

        meidx = findfirst(>(1) ∘ length, _mutation_edges)

        while !isnothing(meidx)
            nextidx = meidx
            live_edges = _mutation_edges[meidx]
            nbp = length(live_edges) - 1

            ## Sample breakpoints ##
            bp_lbound = isone(nextidx) ?
                zero(eltype(positions(arg))) : idxtopos(arg, nextidx - 1)
            bp_ubound = idxtopos(arg, nextidx)
            bp_dist = Uniform(bp_lbound, bp_ubound)
            breakpoints = (sort ∘ rand)(rng, bp_dist, nbp)
            arg.logprob[] -= nbp * log(bp_ubound - bp_lbound)

            ## Constrained recombinations ##
            for breakpoint ∈ breakpoints
                sample_recombination_constrained!(rng, arg, breakpoint,
                                                  winwidth, live_edges,
                                                  buffer = buffer)
            end

            mutation_edges!(_mutation_edges, arg, Ω(first(positions(arg)), ∞), buffer = buffer)
            meidx = findfirst(>(1) ∘ length, _mutation_edges)

            isnothing(meidx) && break
        end
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

    r = nrecombinations(arg, dummy = true)
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
expanded to include space for a pseudo edge. If `arg` has `r`recombinations and
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

function thevenin!(rng, tree, arg::Arg)
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

    # ## To limit allocations, we make `C` large enough to handle the worst case
    # ## scenario. We can pass an appropriate `view` to `impedance!`.
    C = _impedance_C(arg, n - 1, edgesid)
    Z2 = _impedance_Z(arg, edgesid, true)

    live = collect(leaves(tree))
    for nlive ∈ range(length(live), 2, step = -1)
        ## Sample first vertex ##
        v1_idx = sample(rng, 1:nlive)
        tree.logprob[] -= log(nlive)
        v1 = live[v1_idx]
        live[v1_idx] = live[nlive]
        nlive -= 1

        ## Sample second vertex ##
        v1_leaves = descendants(tree, v1) ∩ leaves(tree)
        if isempty(v1_leaves)
            push!(v1_leaves, v1)
        end

        potential = function (v)
            log(
                _thevenin_impedance_helper!(tree, arg, v1_leaves, v, C, Z2,
                    edgesid, estack, vqueue, visited))
        end

        v2_idx, (gumbel_x, gumbel_μ) = _sample_toilet(rng,
            live[1:nlive], potential, 1)
        v2 = live[v2_idx]
        v = 2n - nlive
        live[v2_idx] = v
        tree.logprob[] += logpdf(Gumbel(gumbel_μ), gumbel_x)

        ## Compute coalescence latitude ##
        Δcoal_dist = Exponential(inv(nlive))
        Δcoal = rand(rng, Δcoal_dist)
        tree.logprob[] += logpdf(Δcoal_dist, Δcoal)

        ## Add event to tree ##
        add_edge!(tree, v, v1)
        add_edge!(tree, v, v2)

        sequences(tree)[v] = sequence(tree, v1) & sequence(tree, v2)
        latitudes(tree)[n - nlive] = latitude(tree, v - 1) + Δcoal
    end

    tree
end

thevenin(rng, arg::Arg) = thevenin!(rng, Tree(sam(arg)), arg)

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
    for v ∈ Iterators.flatten((range(n + 1, 2n - 2), range(2n + 1, nv(arg), step = 2)))
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
    end

    ## Edges ##
    for e ∈ edges(arg)
        if isempty(ancestral_intervals(arg, e))
            @info "Empty ancestral interval" e
            flag = false
        end
    end

    flag
end
