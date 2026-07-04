import Base: ==, !=, in, isdisjoint, intersect!, join, union, union!

import IntervalSets: endpoints, leftendpoint, rightendpoint

using IntervalSets: TypedEndpointsInterval

import SIMD._pointer

#          +----------------------------------------------------------+
#          |                     Set of intervals                     |
#          +----------------------------------------------------------+

"""
    $(SIGNATURES)

True if intervals are disconnected, that is their intersection is empty *and*
they do not share an endpoint.

--*Internal*--
"""
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
