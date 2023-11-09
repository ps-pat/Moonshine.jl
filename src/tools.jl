using IntervalSets: AbstractInterval, Interval

import Base:
    union,
    intersect,
    join

######################
# Khatri-Rao Product #
######################

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

####################
# Set of intervals #
####################

function simplify!(xs)
    tmp = Set{eltype(xs)}()

    while !isempty(xs)
        @label main
        x = pop!(xs)
        isempty(x) && continue

        for y ∈ xs
            isempty(x ∩ y) && continue
            push!(xs, x ∪ pop!(xs, y))
            @goto main
        end

        push!(tmp, x)
    end

    while !isempty(tmp)
        push!(xs, pop!(tmp))
    end

    xs
end

union(x::T, xs::Set{T}) where T =
    simplify!(Set([x]) ∪ xs)

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

intersect(x::T, xs::Set{T}) where T =
    (simplify! ∘ Set ∘ broadcast)(Fix1(intersect, x), xs)

intersect(xs::Set{T}, ys::Set{T}) where T =
    (simplify! ∘ mapreduce)(x -> intersect(x, ys), ∪, xs)

for fun ∈ [:union, :intersect]
    @eval $fun(xs::Set{T}, x::T) where T = $fun(x, xs)
end

in(x::T, s::Set{<:AbstractInterval{T}}) where T = any(int -> x ∈ int, s)
