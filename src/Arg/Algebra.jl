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

