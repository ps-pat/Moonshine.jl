import Base:
union,
intersect, intersect!,
join,
in,
issubset,
isdisjoint

import IntervalSets: leftendpoint, rightendpoint, endpoints

const AI = AbstractInterval

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

function simplify!(xs::Set{T}; buffer = default_buffer()) where T
    @no_escape buffer begin
        tmp = @alloc(T, length(xs))
        tmpidx = firstindex(tmp)

        while !isempty(xs)
            x = pop!(xs)
            isempty(x) && continue

            newint = false
            for y ∈ xs
                isempty(x ∩ y) && continue
                newint = true
                yy = x ∪ pop!(xs, y)
                push!(xs, yy)
            end

            if !newint
                tmp[tmpidx] = x
                tmpidx += 1
            end
        end

        for k ∈ 1:(tmpidx-1)
            push!(xs, tmp[k])
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
    union!(deepcopy(As), B, buffer = buffer)

function union!(As::Set{<:AI}, Bs::Set{<:AI}; buffer = default_buffer())
    for B ∈ Bs
        union!(As, B, buffer = buffer, simplify = false)
    end

    simplify!(As, buffer = buffer)
end

union(As::Set{<:AI}, B::Set{<:AI}; buffer = default_buffer()) =
    union!(deepcopy(As), B, buffer = buffer)

# -- Intersection ------------------------------------------------------

function intersect!(As::Set{<:AI}, B::T; buffer = default_buffer()) where T<:AI
    @no_escape buffer begin
        tmp = @alloc(T, length(As))

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
    intersect!(deepcopy(As), B, buffer = buffer)

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
    intersect!(deepcopy(As), Bs, buffer = buffer)

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

function isdisjoint(As::Set{<:AI}, B::T) where T<:AI
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

# -- Helpers -----------------------------------------------------------

for fun ∈ [:union, :intersect, :isdisjoint]
    @eval $fun(A::T, Bs::Set{<:AI}) where T<:AI = $fun(Bs, A)
end

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

endpoints(As::Set{<:AbstractInterval}) = leftendpoint(As), rightendpoint(As)
