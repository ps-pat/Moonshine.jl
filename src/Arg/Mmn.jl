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

export next_inconsistent_idx
"""
    $(SIGNATURES)

Index of the next inconsistent marker, that is the next which mutates more than
once in a given ancestral recombination graph.

`stack` must be of type `CheapStack{Edge{VertexType}}` (see
[`CheapStack`](@ref)).
"""
function next_inconsistent_idx(arg, idx, stack;
                               mutations_edges = ntuple(_ -> Edge{VertexType}[], 8mmn_chunksize),
                               buffer = default_buffer())
    chunksize = 8mmn_chunksize
    ## The mask restricts search to markers in (original) `idx` and
    ## `nmarkers(arg)` inclusively.
    mask = typemax(mmn_chunktype)
    mask <<= idxinchunk(chunksize, idx) - 1

    @inbounds while idx <= nmarkers(arg)
        empty!.(mutations_edges)

        ωlbound = idxtopos(arg, idx)

        idx_chunk = chunkidx(chunksize, idx)
        idx = chunksize * idx_chunk + 1 # idx is now the first marker of the next chunk
        if idx > nmarkers(arg)
            mask >>>= chunksize - idxinchunk(chunksize, nmarkers(arg))
            ωubound = ∞
        else
            ωubound = idxtopos(arg, idx)
        end

        base_ω = Ω(ωlbound, ωubound)

        @inbounds @no_escape buffer begin
            visited = @alloc(Bool, nrecombinations(arg))
            ei_ptr = unsafe_convert(Ptr{Edge{VertexType}},
                                    @alloc_ptr(ne(arg) * sizeof(Edge{VertexType})))

            padded_ne = ne(arg) + ne(arg) % simd_vecsize
            chunks_s_ptr = unsafe_convert(Ptr{mmn_chunktype},
                                          @alloc_ptr(padded_ne * sizeof(VertexType)))
            chunks_d_ptr = unsafe_convert(Ptr{mmn_chunktype},
                                          @alloc_ptr(padded_ne * sizeof(VertexType)))

            ## Collect edges
            ne_interval = zero(Int)
            for e ∈ edges_interval(arg, base_ω, stack, visited)
                s_ptr = unsafe_convert(Ptr{mmn_chunktype},
                                       pointer(sequence(arg, src(e)).data.chunks))
                d_ptr = unsafe_convert(Ptr{mmn_chunktype},
                                       pointer(sequence(arg, dst(e)).data.chunks))
                schunk = unsafe_load(s_ptr, idx_chunk)
                dchunk = unsafe_load(d_ptr, idx_chunk)

                ne_interval += 1
                unsafe_store!(ei_ptr, e, ne_interval)
                unsafe_store!(chunks_s_ptr, schunk, ne_interval)
                unsafe_store!(chunks_d_ptr, dchunk, ne_interval)
            end

            ## Traverse marginal graph and compute mutation sequences
            padded_ne = ne_interval + ne_interval % simd_vecsize
            chunks_s = UnsafeArray{eltype(chunks_s_ptr), 1}(chunks_s_ptr, (padded_ne,))
            chunks_d = UnsafeArray{eltype(chunks_d_ptr), 1}(chunks_d_ptr, (padded_ne,))
            let lane = VecRange{simd_vecsize}(0)
                @inbounds for k ∈ 1:simd_vecsize:padded_ne
                    chunks_s[lane + k] ⊻= chunks_d[lane + k]
                    chunks_s[lane + k] &= mask
                end
            end

            ## Find mutation edges
            for i ∈ 1:ne_interval
                e = unsafe_load(ei_ptr, i)
                mutations_sequence = unsafe_load(chunks_s_ptr, i)
                acc = 0

                while true
                    j = trailing_zeros(mutations_sequence) + 1
                    j > chunksize && break
                    pos = idxtopos(arg, chunksize * (idx_chunk - 1) + acc + j)
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
            mutation_idx = chunksize * (idx_chunk - 1) + idx_mutation_chunk
            mutation_edges = mutations_edges[idx_mutation_chunk]
            return mutation_idx, mutation_edges
        end

        empty!.(mutations_edges)
        mask = typemax(mmn_chunktype)
    end

    0, Edge{VertexType}[]
end
