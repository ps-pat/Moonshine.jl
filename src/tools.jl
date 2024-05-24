using IntervalSets: AbstractInterval, Interval

import Base: union,
    intersect,
    join,
    in,
    issubset,
    isdisjoint

import IntervalSets: leftendpoint, rightendpoint, endpoints

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

function union(A::T, Bs::Set{<:T}) where T<:AbstractInterval
    ret = Set{T}()
    push!(ret, A)
    push!(ret, Bs...)
    simplify!(ret)
end

function union(As::Set{<:T}, Bs::Set{<:T}) where T<:AbstractInterval
    ret = Set{T}()

    for A ∈ As
        push!(ret, A)
    end

    for B ∈ Bs
        push!(ret, B)
    end

    simplify!(ret)
end

function intersect(A::T, Bs::Set{<:T}) where T<:AbstractInterval
    ret = Set{T}()
    for B ∈ Bs
        push!(ret, A ∩ B)
    end

    simplify!(ret)
end

function intersect(As::Set{T}, Bs::Set{T}) where T <:AbstractInterval
    ret = Set{T}()
    for A ∈ As
        for B ∈ Bs
            push!(ret, A ∩ B)
        end
    end

    simplify!(ret)
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

function isdisjoint(A::AbstractInterval, Bs::Set{<:AbstractInterval})
    for B ∈ Bs
        isdisjoint(A, B) || return false
    end

    true
end

function isdisjoint(As::Set{<:AbstractInterval}, Bs::Set{<:AbstractInterval})
    for A ∈ As
        isdisjoint(A, Bs) || return false
    end

    true
end

for fun ∈ [:union, :intersect, :isdisjoint]
    @eval $fun(As::Set{T}, B::T) where T<:AbstractInterval = $fun(B, As)
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
