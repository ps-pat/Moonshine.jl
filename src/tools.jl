Base.@assume_effects :inaccessiblememonly :terminates_globally :effect_free :notaskstate function ⊙(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    c = size(A, 2)
    c == size(B, 2) ||
        throw(DimensionMismatch(lazy"matrix A has dimension $(size(A)), " *
                                    "matrix B has dimension $(size(B))"))

    p, q = size(A, 1), size(B, 1)

    ret = similar(A, T, (p * q, c))

    k = 1
    for j ∈ 1:c
        for iA ∈ 1:p
            @simd for iB ∈ 1:q
                @inbounds ret[k] = A[iA, j] * B[iB, j]
                k += 1
            end
        end
    end

    ret
end
