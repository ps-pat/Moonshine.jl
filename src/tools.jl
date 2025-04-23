import Base:
union, union!,
intersect, intersect!,
join,
in,
issubset,
isdisjoint,
==, !=

import IntervalSets: leftendpoint, rightendpoint, endpoints, width

using IntervalSets: TypedEndpointsInterval

import SIMD._pointer

#          +----------------------------------------------------------+
#          |                    Khatri-Rao Product                    |
#          +----------------------------------------------------------+

function ⊙(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    c = size(A, 2)
    c == size(B, 2) ||
        throw(DimensionMismatch(lazy"matrix A has dimension $(size(A)), " *
                                "matrix B has dimension $(size(B))"))

    p, q = size(A, 1), size(B, 1)

    ret = similar(A, (p * q, c))

    k = 1
    for j ∈ 1:c
        for iA ∈ 1:p
            @inbounds @simd for iB ∈ 1:q
                ret[k] = A[iA, j] * B[iB, j]
                k += 1
            end
        end
    end

    ret
end

#          +----------------------------------------------------------+
#          |                     Set of intervals                     |
#          +----------------------------------------------------------+

function isdisconnected(A::AI, B::AI)
    AB = A ∩ B
    isempty(AB) && !=(endpoints(AB)...)
end

for idx ∈ 0:15
    L1 = isone(idx & 1) ? :closed : :open
    idx >>= 1
    R1 = isone(idx & 1) ? :closed : :open
    idx >>= 1
    L2 = isone(idx & 1) ? :closed : :open
    idx >>= 1
    R2 = isone(idx & 1) ? :closed : :open

    TypeA = TypedEndpointsInterval{L1, R1}
    TypeB = TypedEndpointsInterval{L2, R2}

    cmp1 = all(==(:closed), (R1, L2)) ? :(<) : :(<=)
    cmp2 = all(==(:closed), (L1, R2)) ? :(>) : :(>=)

    @eval isdisjoint(A::$TypeA, B::$TypeB) =
        $cmp1(rightendpoint(A), leftendpoint(B)) ||
        $cmp2(leftendpoint(A), rightendpoint(B))
end

#          +----------------------------------------------------------+
#          |                   UnsafeArrays & simd                    |
#          +----------------------------------------------------------+

Base.@propagate_inbounds _pointer(arr::UnsafeArray, i, I) =
    pointer(arr, LinearIndices(arr)[i, I...])
