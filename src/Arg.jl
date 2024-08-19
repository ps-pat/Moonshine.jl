using Graphs

import Graphs: add_vertices!, add_edge!, rem_edge!

using Random

using StatsBase: samplepair

using SparseArrays

using Combinatorics: combinations

##################
# Arg Definition #
##################

export Arg
mutable struct Arg <: AbstractGenealogy
    graph::SimpleDiGraph{VertexType}
    latitudes::Vector{Float64}
    sequences::Vector{Sequence}
    ancestral_intervals::Dict{Edge, Set{Ω}}
    sample::Sample
    logprob::Base.RefValue{Float64x2}
end

Arg(tree::Tree) = Arg(
    graph(tree),
    latitudes(tree),
    sequences(tree),
    Dict{Edge, Set{Ω}}(),
    sam(tree),
    Ref(prob(tree, logscale = true))
)

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

function nrecombinations(arg::Arg; dummy = false)
    ret = ne(arg) - nv(arg) + 1

    if !dummy
        v = 2nleaves(arg)
        while isinf(recbreakpoint(arg, v))
            ret -= 1
            v += 2
        end
    end

    ret
end

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

########################
# Ancestrality Methods #
########################

function ancestral_intervals!(ωs, arg::Arg, e::Edge; wipe = true)
    wipe && empty!(ωs)

    haskey(arg.ancestral_intervals, e) || return push!(ωs, Ω(0, ∞))

    for ω ∈ arg.ancestral_intervals[e]
        push!(ωs, ω)
    end

    ωs
end

ancestral_intervals(arg::Arg, e::Edge) =
    get(() -> Set{Ω}((Ω(0, ∞),)), arg.ancestral_intervals, e)

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

ancestral_mask!(η, ωs, arg::Arg, x::Union{VertexType, Edge{VertexType}}; wipe = true) =
    ancestral_mask!(η, sam(arg), ancestral_intervals!(ωs, arg, x), wipe = wipe)

ancestral_mask(arg::Arg, x::Union{VertexType, Edge{VertexType}}) =
    ancestral_mask!(Sequence(undef, nmarkers(arg)), Set{Ω}(), arg, x)

export recbreakpoint
"""
    recbreakpoint(arg, v)

Returns the breakpoint associated with vertex `v`. If `v` is not a
recombination vertex, returns ∞. If `v` is not a "true" recombination vertex
i.e. one of its outgoing edge is non ancestral, returns ∞.
"""
function recbreakpoint(arg::Arg, v)
    ret = ∞
    isrecombination(arg, v) || return ret

    _dads = dads(arg, v)
    ωs1 = ancestral_intervals(arg, Edge(first(_dads) => v))
    ωs2 = ancestral_intervals(arg, Edge(last(_dads) => v))

    (isempty(ωs1) || isempty(ωs2)) && return ∞

    for (p1, p2) ∈ Iterators.product(endpoints.((ωs1, ωs2))...)
        p1 == p2 && return p1
    end
end

#######
# MMN #
#######

function mutationsidx!(res, mask, ωs, arg, e,
                       firstchunk, firstidx, lastchunk, lastidx)
    η1, η2 = sequences(arg, e)
    ancestral_mask!(mask, ωs, arg, e)
    marker_mask = one(UInt64) << (firstidx - 1)
    idx = one(Int)

    @inbounds for k ∈ range(firstchunk, lastchunk)
        xored_chunk = (η1.data.chunks[k] ⊻ η2.data.chunks[k]) &
            mask.data.chunks[k]

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
function mutation_edges!(mutations, arg, ω::Ω)
    ## Compute the chunks and indices.
    idx = postoidx(arg, ω)
    lidx, ridx = first(idx), last(idx)

    firstchunk = chunkidx(Sequence, lidx)
    firstidx = idxinchunk(Sequence, lidx)
    lastchunk = chunkidx(Sequence, ridx)
    lastidx = idxinchunk(Sequence, ridx)

    @inbounds for k ∈ eachindex(mutations)
        resize!(mutations[k], 0)
    end

    mask = Sequence(undef, nmarkers(arg))
    ωs = Set{Ω}()
    @inbounds for edge ∈ edges_interval(arg, ω)
        mutationsidx!(mutations, mask, ωs, arg, edge,
                      firstchunk, firstidx, lastchunk, lastidx)
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
            arg.ancestral_intervals[Edge(dad, v)] =
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
    root_recombination = !rem_edge!(arg, cedge)
    root_recombination ||
        add_edge!(arg, Edge(src(cedge), cvertex), ωc ∪ (ωr ∩ Ω(breakpoint, ∞)))

    ## Compute sequence of new vertices.
    let mask = Sequence(undef, nmarkers(arg)),
        ωs = Set{Ω}()

        _compute_sequence!(arg, rvertex, mask, ωs)
        _compute_sequence!(arg, cvertex, mask, ωs)
    end

    ## Update sequences and ancetral intervals.
    update_upstream!(arg, src(redge))
    root_recombination || update_upstream!(arg, src(cedge))

    arg
end

function build!(rng, arg::Arg; winwidth = ∞)
    λc = rec_rate(arg, false)

    _mutation_edges = [sizehint!(Vector{Edge{VertexType}}(undef, 0), nv(arg) ÷ 2)
                       for _ ∈ 1:nmarkers(arg)]
    mutation_edges!(_mutation_edges, arg, Ω(first(positions(arg)), ∞))

    mepos = findfirst(>(1) ∘ length, _mutation_edges)
    isnothing(mepos) && return arg # ARG is already consistent.
    nextidx = mepos
    breakpoint = isone(nextidx) ? zero(Float64) : idxtopos(arg, nextidx - 1)

    while !isnothing(mepos)
        # live_edges = _mutation_edges[nextidx + length(_mutation_edges) - nmarkers(arg)]
        live_edges = _mutation_edges[mepos]

        ## Sample next breakpoint.
        nextpos = idxtopos(arg, nextidx)
        rint = Ω(breakpoint, nextpos)
        @debug "Recombination bounds" endpoints(rint)
        l = width(rint)

        if !iszero(l)
            b = branchlength(arg, rint)
            λ = λc * b
            Δbp_dist = truncated(Exponential(inv(λ)), upper = l)
            Δbp = rand(rng, Δbp_dist)
            breakpoint += Δbp
            arg.logprob[] += logpdf(Δbp_dist, Δbp)
        end

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
            arg.logprob[] += logpdf(rlat_dist, rlat)
        else
            rlat = first(rlat_bounds)
        end

        ## Sample recoalescence latitude.
        clat_bounds = (max(rlat, latitude(arg, dst(cedge))),
                       latitude(arg, src(cedge)))
        if first(clat_bounds) < last(clat_bounds)
            clat_dist = Uniform(clat_bounds...)
            clat = rand(rng, clat_dist)
            arg.logprob[] += logpdf(clat_dist, clat)
        else
            clat = first(clat_bounds)
        end

        ## Add recombination event to the graph.
        @debug "Recombination event" redge, cedge, breakpoint, rlat, clat
        recombine!(arg, redge, cedge, breakpoint, rlat, clat)

        ## Update variables.
        mutation_edges!(_mutation_edges, arg, Ω(breakpoint, ∞))
        mepos = findfirst(>(1) ∘ length, _mutation_edges)
        isnothing(mepos) && break
        nextidx = nextidx + mepos - 1

        if !isone(mepos)
            breakpoint = idxtopos(arg, nextidx - 1)
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

    # any(v -> isleaf(arg, v), (s, d)) ||
    #     push!(visited, first(children(arg, s) ∩ children(arg, d)))
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
function cbasis! end,
function cbasis end

function cbasis!(vec, arg::Arg, v::VertexType, lk = Threads.ReentrantLock();
                 estack = Stack{Edge{VertexType}}(ceil(Int, log(nv(arg)))),
                 edgesid = Dict(reverse.(enumerate(edges(arg)))),
                 vqueue = Queue{VertexType}(ceil(Int, log(nv(arg)))),
                 visited = Set{VertexType}())
    empty!(estack)
    empty!(vqueue)
    empty!(visited)

    ## Add the recoalescence edge to vec.
    _dads = minmax(dads(arg, v)...)
    Threads.lock(lk) do
        vec[edgesid[Edge(last(_dads) => v)]] = -1
    end

    _path_bfs_forward!(estack, arg, _dads..., vqueue, visited)
    _bfs_backtrack!(vec, edgesid, estack, lk)
end

cbasis(arg::Arg, v::VertexType, lk = Threads.ReentrantLock();
       estack = Stack{Edge{VertexType}}(ceil(Int, log(nv(arg)))),
       edgesid = Dict(reverse.(enumerate(edges(arg)))),
       vqueue = Queue{VertexType}(ceil(Int, log(nv(arg)))),
       visited = Set{VertexType}()) =
    cbasis!(zeros(Float64, ne(arg)), arg, v, lk,
            estack = estack, edgesid = edgesid,
            vqueue = vqueue, visited = visited)

function cbasis!(mat, arg::Arg;
                 edgesid = Dict(reverse.(enumerate(edges(arg)))))
    fill!(mat, 0)

    r = nrecombinations(arg, dummy = true)
    n = nleaves(arg)

    lk = Threads.ReentrantLock()
    rec_offset = ne(arg) - nv(arg) + 1 - nrecombinations(arg)

    tasks = map(chunks(range(1, length = r - rec_offset),
                       n = Threads.nthreads(),
                       split = :scatter)) do ks
            Threads.@spawn begin
                local estack = Stack{Edge}(ceil(Int, log(nv(arg))))
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

cbasis(arg::Arg; edgesid = Dict(reverse.(enumerate(edges(arg))))) =
    cbasis!(spzeros(ne(arg), nrecombinations(arg)), arg, edgesid = edgesid)

function _thevenin_R(arg::Arg, edgesmap, take_sqrt = true)
    R = Diagonal(Vector{Float64}(undef, ne(arg)))

    for (e, idx) ∈ edgesmap
        R[idx, idx] = impedance(arg, e)
    end

    take_sqrt || return R

    for k ∈ eachindex(R.diag)
        R.diag[k] = sqrt(R.diag[k])
    end

    R
end

"""
    _thevenin_C(arg, edgesmap)

Return the generator of the cycle space of an ancestral recombination graph
expanded to include space for a pseudo edge. If `arg` has r recombinations and
k edges, an (k+2) x r+1 matrix is returned. The generator is stored in the
upper k x r block. Remaining entriez are initialized to 0.
"""
function _thevenin_C(arg, edgesmap)
    C = spzeros(ne(arg), nrecombinations(arg) + 1)
    cbasis!(C, arg, edgesid = edgesmap)
    C
end

"""
    _thevenin_update_C!(C, arg, s, d, edgesmap, estack, vqueue, visited)

Store the pseudo-cycle generated by the pair of vertices (s, d) to the last
column of `C`.
"""
function _thevenin_update_C!(C, arg, s, d, edgesmap, estack, vqueue, visited)
    empty!(estack)
    empty!(vqueue)
    empty!(visited)

    vec = @view C[:,end]
    fill!(vec, 0)

    _path_bfs_forward!(estack, arg, s, d, vqueue, visited)
    _bfs_backtrack!(vec, edgesmap, estack, Threads.ReentrantLock())
end

export thevenin!, thevenin
function thevenin!(arg::Arg, s, d, C, R2;
                   edgesmap = Dict(reverse.(enumerate(edges(arg)))),
                   estack = Stack{Edge}(ceil(Int, log(nv(arg)))),
                   vqueue = Queue{VertexType}(ceil(Int, log(nv(arg)))),
                   visited = Set{VertexType}())
    empty!(estack)
    empty!(vqueue)
    empty!(visited)

    r = nrecombinations(arg)

    _thevenin_update_C!(C, arg, s, d, edgesmap, estack, vqueue, visited)

    U = qr(R2 * C).R
    invU = (inv ∘ Matrix)(U)
    current = (invU * invU')
    inv(last(current))
end

function thevenin(arg::Arg, s, d;
                  estack = Stack{Edge}(ceil(Int, log(nv(arg)))),
                  edgesmap = Dict(reverse.(enumerate(edges(arg)))),
                  vqueue = Queue{VertexType}(ceil(Int, log(nv(arg)))),
                  visited = Set{VertexType}())
    empty!(estack)
    empty!(vqueue)
    empty!(visited)

    ## Square root of impedances.
    R2 = _thevenin_R(arg, edgesmap, true)

    ## Compute fundamental basis
    C = _thevenin_C(arg, edgesmap)

    thevenin!(arg, s, d, C, R2,
              edgesmap = edgesmap, estack = estack,
              vqueue = vqueue, visited = visited)
end

function thevenin(arg::Arg;
                  estack = Stack{Edge}(ceil(Int, log(nv(arg)))),
                  edgesmap = Dict(reverse.(enumerate(edges(arg)))),
                  vqueue = Queue{VertexType}(ceil(Int, log(nv(arg)))),
                  visited = Set{VertexType}())
    empty!(estack)
    r = nrecombinations(arg)

    C = _thevenin_C(arg, edgesmap)
    R2 = _thevenin_R(arg, edgesmap, true)

    Iterators.map(combinations(1:nleaves(arg), 2)) do (s, d)
        thevenin!(arg, s, d, C, R2,
                  estack = estack, edgesmap = edgesmap,
                  vqueue = vqueue, visited = visited)
    end
end

export thevenin_matrix
function thevenin_matrix(arg::Arg,
                         estack = Stack{Edge}(ceil(Int, log(nv(arg)))),
                         edgesmap = Dict(reverse.(enumerate(edges(arg)))))
    empty!(estack)

    n = nleaves(arg)
    r = nrecombinations(arg)

    C = _thevenin_C(arg, edgesmap)
    R2 = _thevenin_R(arg, edgesmap, true)
    g = zeros(Float64, ne(arg))

    mat = Matrix{Float64}(undef, n, n)
    mat[diagind(mat)] .= 0
    for (j, i) ∈ combinations(1:n, 2)
        mat[i, j] = thevenin!(arg, i, j, C, R2,
                              estack = estack, edgesmap = edgesmap)
    end

    Symmetric(mat, :L)
end
