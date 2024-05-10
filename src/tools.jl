using IntervalSets: AbstractInterval, Interval

import Base: union,
    intersect,
    join,
    in,
    issubset,
    isdisjoint

using Base: Fix1

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

function simplify!(xs::Set{T}) where T
    tmp = Set{T}()

    while !isempty(xs)
        @label main
        x = pop!(xs)
        isempty(x) && continue

        for y ∈ xs
            isempty(x ∩ y) && continue
            yy = x ∪ pop!(xs, y)
            push!(xs, yy)
            @goto main
        end

        push!(tmp, x)
    end

    while !isempty(tmp)
        push!(xs, pop!(tmp))
    end

    xs
end

union(x::T, xs::Set{T}) where T = simplify!(Set([x]) ∪ xs)

function union(x::Set{T}, y::Set{T}) where T<:Interval
    ret = empty(x)
    for el ∈ x
        push!(ret, el)
    end

    for el ∈ y
        push!(ret, el)
    end

    simplify!(ret)
end

intersect(x::T, xs::Set{T}) where T = (simplify! ∘ Set ∘ broadcast)(Fix1(intersect, x), xs)

function intersect(xs::Set{T}, ys::Set{T}) where T
    (simplify! ∘ mapreduce)(x -> intersect(x, ys), ∪, xs, init = Set{T}())
end

for fun ∈ [:union, :intersect]
    @eval $fun(xs::Set{T}, x::T) where T = $fun(x, xs)
end

in(x, As::Set{<:AbstractInterval}) = any(A -> x ∈ A, As)

function issubset(As::Set{<:AbstractInterval}, Bs::Set{<:AbstractInterval})
    for A ∈ As
        A ⊆ Bs || return false
    end
    true
end

issubset(A::AbstractInterval, Bs::Set{<:AbstractInterval}) = any(B -> A ⊆ B, Bs)

issubset(As::Set{<:AbstractInterval}, B::AbstractInterval) = all(A -> A ⊆ B, As)

function isdisjoint(As::Set{<:AbstractInterval}, Bs::Set{<:AbstractInterval})
    for A ∈ As
        for B ∈ Bs
            isdisjoint(A, B) || return false
        end
    end

    true
end
