import Base:
    ## Iteration
    iterate, length, size, eltype, IteratorSize,
    ## Indexing
    getindex, setindex!, firstindex, lastindex,
    ## AbstractArrays
    IndexStyle, deleteat!, resize!,
    strides, unsafe_convert, elsize,
    ## Set operations
    union!, union, intersect!, intersect, join, in, issubset, isdisjoint,
    copy, empty!, ==, !=

import IntervalSets: leftendpoint, rightendpoint, endpoints, width

export Ω
"""
    const Ω = Interval{:closed, :open, Float64}

Right semi-open interval.

See also [`AI`](@ref), [`AIs`](@ref) and [`AncestralIntervals`](@ref).
"""
const Ω = Interval{:closed, :open, Float64}

#          +----------------------------------------------------------+
#          |                     Type definition                      |
#          +----------------------------------------------------------+

export AncestralIntervals
"""
    $(TYPEDEF)

Collection of intervals.

Meant to represent the set of intervals an edge/vertex is ancestral for. You
might want to use the convenient shorthand [`AIs`](@ref) instead.

Implements [the iteration interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-iteration) and [the array interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array).

See also [`AI`](@ref) and [`Ω`](@ref).

# Fields
$(TYPEDFIELDS)

# Constructors
$(METHODLIST)

# Arguments
If `simplify = true`, intervals contained in data are simplified: see
[`simplify!`](@ref) for details.

!!! warning
    Many methods assume `AIs` to be simplified. You might want to disable
    simplification to optimize a sequence of operation, but you should probably
    simplify the final result.
"""
struct AncestralIntervals{T<:AbstractVector{<:AI}} <: AbstractVector{AI}
    data::T

    function AncestralIntervals{T}(data::T; simplify = true) where T
        ωs = new{T}(data)
        simplify && simplify!(ωs)
        ωs
    end
end

"""
    const AIs = AncestralIntervals

Alias for [`AncestralIntervals`](@ref).

See also [`AI`](@ref) and [`Ω`](@ref).
"""
const AIs = AncestralIntervals

AIs{V}() where V<:AbstractVector{T} where T = AIs(T[], simplify = false)

AIs(data::V; simplify = true) where V<:AbstractVector{T} where T =
    AIs{V}(data, simplify = simplify)

#          +----------------------------------------------------------+
#          |                        Interfaces                        |
#          +----------------------------------------------------------+

# -- Iteration ---------------------------------------------------------

for f ∈ (:iterate, :length, :size)
    @eval $f(ai::AIs) = $f(ai.data)
end

iterate(ai::AIs, state) = iterate(ai.data, state)

eltype(::Type{<:AIs{<:AbstractVector{T}}}) where T = T

IteratorSize(::AIs) = Base.HasLength()
IteratorSize(::Type{AIs}) = Base.HasLength()

# -- Indexing ----------------------------------------------------------

for f ∈ (:firstindex, :lastindex)
    @eval $f(ai::AIs) = $f(ai.data)
end

getindex(ai::AIs, i) = getindex(ai.data, i)

setindex!(ai::AIs, v, i) = setindex!(ai.data, v, i)

# -- AbstractArrays ----------------------------------------------------

IndexStyle(::Type{AIs}) = Base.IndexLinear()

deleteat!(ai::AIs, i) = deleteat!(ai.data, i)

resize!(ai::AIs, n) = resize!(ai.data, n)

# -- Strided Arrays ----------------------------------------------------

strides(ai::AIs) = strides(ai.data)

unsafe_convert(::Type{Ptr{T}}, ai::AIs) where T =
    unsafe_convert(Ptr{T}, ai.data)

elsize(::Type{AIs{T}}) where T = elsize(T)

#          +----------------------------------------------------------+
#          |                    Essential Methods                     |
#          +----------------------------------------------------------+

copy(ωs::AIs) = AIs(copy(ωs.data))

empty!(ωs::AIs) = empty!(ωs.data)

#          +----------------------------------------------------------+
#          |                      Set Operations                      |
#          +----------------------------------------------------------+

"""
    $(SIGNATURES)

Simplify an [`AIs`](@ref).

Two operations are performed:
* connected intervals are merged together (see [`isdisconnected`](@ref));
* intervals are sorted by left endpoint.

--*Internal*--
"""
function simplify!(ωs::AIs)
    filter!(!isempty, ωs)

    @label loop
    @inbounds for i ∈ length(ωs):-1:2
        for j ∈ (i-1):-1:1
            isdisconnected(ωs[i], ωs[j]) && continue
            ωs[j] = ωs[j] ∪ ωs[i]
            deleteat!(ωs, i)
            @goto loop
        end
    end

    sort!(ωs, alg = InsertionSort, by = leftendpoint)
end

# -- Union -------------------------------------------------------------

function union!(ai::AIs, ω::T; simplify = true) where T<:AI
    push!(ai.data, ω)
    simplify && simplify!(ai)
    ai
end

function union!(ωs1::AIs, ωs2::AIs; simplify = true)
    append!(ωs1.data, ωs2.data)
    simplify && simplify!(ωs1)
    ωs1
end

union(ωs1::AIs, ωs2::AIs; simplify = true) =
    union!(copy(ωs1), ωs2, simplify = simplify)

# -- Intersection ------------------------------------------------------

function intersect!(ωs::AIs, ω::T; simplify = true) where T<:AI
    @inbounds @simd for k ∈ eachindex(ωs)
        ωs[k] = ωs[k] ∩ ω
    end

    simplify && simplify!(ωs)
    ωs
end

function intersect!(ωs1::AIs{<:AbstractVector{T}}, ωs2::AIs;
                    simplify = true, buffer = default_buffer()) where T
    @no_escape buffer begin
        n1, n2 = length(ωs1), length(ωs2)
        ωbuf_ptr = convert(Ptr{T}, @alloc_ptr(n1 * n2 * sizeof(T)))
        ωbuf_len = 0

        @inbounds for k ∈ (n1*n2):-1:1
            newω = ωs1[((k - 1) % n1) + 1] ∩ ωs2[((k - 1) ÷ n1) + 1]
            isempty(newω) && continue
            ωbuf_len += 1
            unsafe_store!(ωbuf_ptr, newω, ωbuf_len)
        end

        resize!(ωs1, ωbuf_len)
        for k ∈ eachindex(ωs1)
            ωs1[k] = unsafe_load(ωbuf_ptr, k)
        end

        simplify && simplify!(ωs1)
        ωs1
    end
end

intersect(ωs1::AIs, ωs2::AIs; simplify = true, buffer = default_buffer()) =
    intersect!(copy(ωs1), ωs2, simplify = simplify, buffer = buffer)

# -- Inclusion & cie ---------------------------------------------------

in(x, ωs::AIs) = any(ω -> x ∈ ω, ωs)

issubset(x, ωs::AIs) = any(ω -> x ⊆ ω, ωs)

function isdisjoint(x, ωs::AIs)
    nb_not_disjoint = zero(Int)
    ## Possible race condition is not an issue
    @inbounds @simd ivdep for ω ∈ ωs
        nb_not_disjoint += !isdisjoint(x, ω)
    end

    iszero(nb_not_disjoint)
end

for f ∈ (:issubset, :isdisjoint)
    @eval $f(ωs::AIs, x::Union{<:AIs, AI}) = $f(x, ωs)
    @eval $f(ωs1::AIs, ωs2::AIs) = all(ω -> $f(ωs2, ω), ωs1)
end

# -- Allocating Methods ------------------------------------------------

for f ∈ (:union, :intersect)
    f! = Symbol(string(f) * '!')

    @eval $f(ωs::AIs, x::Union{<:AIs, AI}; simplify = true) =
        $f!(copy(ωs), x, simplify = simplify)

    @eval $f(x::Union{<:AIs, AI}, ωs::AIs; simplify = true) = $f(ωs, x, simplify = simplify)
end

#          +----------------------------------------------------------+
#          |                     Interval Methods                     |
#          +----------------------------------------------------------+

leftendpoint(ωs::AIs) = (leftendpoint ∘ first)(ωs)
rightendpoint(ωs::AIs) = maximum(rightendpoint, ωs)

endpoints(ωs::AIs) = leftendpoint(ωs), rightendpoint(ωs)

width(ωs::AIs) = rightendpoint(ωs) - leftendpoint(ωs)

export closure
"""
    $(FUNCTIONNAME)(x)

Mathematical closure of `x`

# Methods
$(METHODLIST)
"""
function closure end

closure(ωs::AIs{V}) where V<:AbstractVector{T} where T<:AI{S} where S =
    ClosedInterval{S}(endpoints(ωs)...)
