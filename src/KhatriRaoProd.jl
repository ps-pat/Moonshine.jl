function ⊙(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    r = size(A, 2)
    r == size(B, 2) ||
        throw(DimensionMismatch(lazy"matrix A has dimension $(size(A)), " *
                                    "matrix B has dimension $(size(B))"))

    p, q = size(A, 1), size(B, 1)

    ret = similar(A, T, (p * q, r))

    @inbounds for (k, (α, β)) ∈ (enumerate ∘ zip)(eachcol(A), eachcol(B))
        ret[:,k] .= kron(α, β)
    end

    ret
end
