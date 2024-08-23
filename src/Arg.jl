using Graphs

import Graphs: add_vertices!, add_edge!, rem_edge!

using Random

using StatsBase: samplepair

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
    sequences::Vector{Sequence}
    ancestral_intervals::Dict{Edge{VertexType}, Set{Ω}}
    sample::Sample
    logprob::Base.RefValue{Float64x2}
end

Arg(tree::Tree) = Arg(
    graph(tree),
    latitudes(tree),
    sequences(tree),
    Dict{Edge{VertexType}, Set{Ω}}(),
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
    wipe && empty!(ωs)
    isleaf(arg, v) && return push!(ωs, Ω(0, ∞))

    for child ∈ children(arg, v)
        ancestral_intervals!(ωs, arg, v, child, wipe = false)
    end

    simplify!(ωs)
end

ancestral_intervals(arg::Arg, v::VertexType) = ancestral_intervals!(Set{Ω}(), arg, v)

#! format: off
ancestral_mask!( η, ωs, arg::Arg, x::Union{VertexType, Edge{VertexType}}; wipe = true) =
    ancestral_mask!(η, sam(arg), ancestral_intervals!(ωs, arg, x), wipe = wipe)
#! format: on

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

function mutationsidx!(res, mask, ωs, arg, e, firstchunk, firstidx, lastchunk)
    η1, η2 = sequences(arg, e)
    ancestral_mask!(mask, ωs, arg, e)
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
            firstchunk, firstidx, lastchunk)
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
    add_vertices!(
        arg,
        (Sequence(trues(nmarkers(arg))), Sequence(trues(nmarkers(arg)))),
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
    let mask = Sequence(undef, nmarkers(arg)), ωs = Set{Ω}()
        _compute_sequence!(arg, rvertex, mask, ωs)
        _compute_sequence!(arg, cvertex, mask, ωs)
    end

    ## Update sequences and ancetral intervals.
    update_upstream!(arg, src(redge))
    root_recombination || update_upstream!(arg, src(cedge))

    arg
end

function sample_recombination_constrained!(rng, arg, breakpoint, live_edges)
    ## Sample recombination and recoalescence edges ##
    e1, e2 = samplepair(rng, length(live_edges))
    redge = live_edges[e1]
    cedge = live_edges[e2]
    if latitude(arg, dst(redge)) > latitude(arg, dst(cedge))
        redge, cedge = cedge, redge
    end

    ## Sample recombination latitude ##
    rlat_bounds = (latitude(arg, dst(redge)) + eps(Float64),
        min(latitude(arg, src(redge)), latitude(arg, src(cedge))))
    if first(rlat_bounds) < last(rlat_bounds)
        rlat_dist = Uniform(rlat_bounds...)
        rlat = rand(rng, rlat_dist)
        arg.logprob[] += logpdf(rlat_dist, rlat)
    else
        rlat = first(rlat_bounds)
    end

    ## Sample recoalescence latitude ##
    clat_bounds =
        (max(rlat, latitude(arg, dst(cedge))), latitude(arg, src(cedge)))
    if first(clat_bounds) < last(clat_bounds)
        clat_dist = Uniform(clat_bounds...)
        clat = rand(rng, clat_dist)
        arg.logprob[] += logpdf(clat_dist, clat)
    else
        clat = first(clat_bounds)
    end

    ## Add recombination event to the graph ##
    @debug "Recombination event" redge, cedge, breakpoint, rlat, clat
    recombine!(arg, redge, cedge, breakpoint, rlat, clat)
end

function build!(rng, arg::Arg)
    λc = rec_rate(arg, false)

    _mutation_edges = [
        sizehint!(Vector{Edge{VertexType}}(undef, 0), nv(arg) ÷ 2) for
        _ ∈ 1:nmarkers(arg)
    ]
    mutation_edges!(_mutation_edges, arg, Ω(first(positions(arg)), ∞))

    mepos = findfirst(>(1) ∘ length, _mutation_edges)
    isnothing(mepos) && return arg # ARG is already consistent.
    nextidx = mepos
    breakpoint = isone(nextidx) ? zero(Float64) : idxtopos(arg, nextidx - 1)

    while !isnothing(mepos)
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

        sample_recombination_constrained!(rng, arg, breakpoint, live_edges)

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
