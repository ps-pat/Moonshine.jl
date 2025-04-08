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

