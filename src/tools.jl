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

######################
# Khatri-Rao Product #
######################

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

####################
# Set of intervals #
####################

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

    @eval function isdisjoint(A::$TypeA, B::$TypeB)
        $cmp1(rightendpoint(A), leftendpoint(B)) && return true
        $cmp2(leftendpoint(A), rightendpoint(B)) && return true
        false
    end
end

function simplify!(xs::Set{T}; buffer = default_buffer()) where T
    @no_escape buffer begin
        tmp_ptr = convert(Ptr{T}, @alloc_ptr(length(xs) * sizeof(T)))
        tmp_len = 0

        while !isempty(xs)
            x = pop!(xs)
            isempty(x) && continue

            newint = false
            for y ∈ xs
                isdisconnected(x, y) && continue
                newint = true
                yy = x ∪ pop!(xs, y)
                push!(xs, yy)
            end

            if !newint
                tmp_len += 1
                unsafe_store!(tmp_ptr, x, tmp_len)
            end
        end

        for k ∈ 1:tmp_len
            push!(xs, unsafe_load(tmp_ptr, k))
        end
    end

    xs
end

# -- Union -------------------------------------------------------------

function union!(As::Set{<:AI}, B::T;
                buffer = default_buffer(), simplify = true) where T<:AI
    @no_escape buffer begin
        tmp = @alloc(T, length(As))

        @inbounds for k ∈ eachindex(tmp)
            A = pop!(As)
            if isdisconnected(A, B)
                tmp[k] = A
            else
                tmp[k] = A ∪ B
            end
        end

        push!(As, B)

        @inbounds for k ∈ eachindex(tmp)
            push!(As, tmp[k])
        end
    end

    simplify && simplify!(As, buffer = buffer)
    As
end

union(As::Set{<:AI}, B::T; buffer = default_buffer()) where T<:AI =
    union!(copy(As), B, buffer = buffer)

function union!(As::Set{<:AI}, Bs::Set{<:AI}; buffer = default_buffer())
    for B ∈ Bs
        union!(As, B, buffer = buffer, simplify = false)
    end

    simplify!(As, buffer = buffer)
end

union(As::Set{<:AI}, B::Set{<:AI}; buffer = default_buffer()) =
    union!(copy(As), B, buffer = buffer)

# -- Intersection ------------------------------------------------------

function intersect!(As::Set{Ta}, B::Tb; buffer = default_buffer()) where {Ta<:AI, Tb<:AI}
    @no_escape buffer begin
        tmp = @alloc(Ta, length(As))

        @inbounds for k ∈ 1:length(As)
            tmp[k] = pop!(As) ∩ B
        end

        @inbounds for k ∈ eachindex(tmp)
            push!(As, tmp[k])
        end
    end

    simplify!(As, buffer = buffer)
end

intersect(As::Set{<:AI}, B::T; buffer = default_buffer()) where T<:AI =
    intersect!(copy(As), B, buffer = buffer)

function intersect!(As::Set{T}, Bs::Set{<:AI}; buffer = default_buffer()) where T<:AI
    @no_escape buffer begin
        tmp = @alloc(T, length(As) * length(Bs))
        tmp_ptr = firstindex(tmp)

        @inbounds while !isempty(As)
            A = pop!(As)
            for B ∈ Bs
                AB = A ∩ B
                isempty(AB) && continue
                tmp[tmp_ptr] = AB
                tmp_ptr += 1
            end
        end

        @inbounds for k ∈ 1:(tmp_ptr-1)
            push!(As, tmp[k])
        end
    end

    simplify!(As, buffer = buffer)
end

intersect(As::Set{<:AI}, Bs::Set{<:AI}; buffer = default_buffer()) =
    intersect!(copy(As), Bs, buffer = buffer)

# -- Inclusion ---------------------------------------------------------

in(x, As::Set{<:AI}) = any(A -> x ∈ A, As)

function issubset(As::Set{<:AI}, Bs::Set{<:AI})
    for A ∈ As
        A ⊆ Bs || return false
    end
    true
end

issubset(A::AI, Bs::Set{<:AI}) = any(B -> A ⊆ B, Bs)

issubset(As::Set{<:AI}, B::AI) = all(A -> A ⊆ B, As)

function isdisjoint(As::Set{<:AI{T}}, B::AI{T}) where T
    for A ∈ As
        isdisjoint(A, B) || return false
    end

    true
end

function isdisjoint(As::Set{<:AI}, Bs::Set{<:AI})
    for A ∈ As
        isdisjoint(Bs, A) || return false
    end

    true
end

function ==(As::Set{<:AI}, Bs::Set{<:AI})
    issubset(As, Bs) || return false
    issubset(Bs, As) || return false
    true
end

!=(As::Set{<:AI}, Bs::Set{<:AI}) = !(As == Bs)

# -- Helpers -----------------------------------------------------------

for fun ∈ [:union, :intersect, :isdisjoint]
    @eval $fun(A::T, Bs::Set{<:AI}) where T<:AI = $fun(Bs, A)
end

# -- Endpoints ---------------------------------------------------------

for (fun, op) ∈ Dict(:leftendpoint => :<, :rightendpoint => :>)
    @eval function $fun(As::Set{<:AbstractInterval})
        ret, As_rest = Iterators.peel(As)
        ret = $fun(ret)
        for A ∈ As_rest
            endpoint = $fun(A)
            if $op(endpoint, ret)
                ret = endpoint
            end
        end
        ret
    end
end

function endpoints(As::Set{<:AI})
    ω, ωs = Iterators.peel(As)
    left, right = endpoints(ω)

    for ω ∈ ωs
        newleft, newright = endpoints(ω)

        if newleft < left
            left = newleft
        end

        if newright > right
            right = newright
        end

    end

    left, right
end

function width(As::Set{<:AI})
    ep = endpoints(As)
    last(ep) - first(ep)
end

export closure
"""
    closure(x)

Mathematical closure of `x`
"""
function closure end

closure(As::Set{<:AI{T}}) where T = ClosedInterval{T}(endpoints(As)...)
