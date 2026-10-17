import Base: eltype, isempty, empty!, first, length, pop!, push!, show, iterate

export CheapStack
"""
    $(TYPEDEF)

Simple stack container.

# Functionalities

The following operations are supported:

  - `isempty`
  - `empty!`
  - `first`
  - `length`
  - `pop!`
  - `push!`

For convenience, `CheapStack` implements [the iteration interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-iteration).

# Fields

$(TYPEDFIELDS)
"""
mutable struct CheapStack{T}
    ptr::Ptr{T}
    beginning::Ptr{T}
end

"""
    $(FUNCTIONNAME)(ptr::Ptr{T})
    $(FUNCTIONNAME){T}(ptr::Ptr{Nothing})

Construct a [`CheapStack`](@ref) for storing objects of type `T`.
"""
CheapStack(ptr) = CheapStack(typeof(ptr)(0), ptr)

CheapStack{T}(ptr::Ptr{Nothing}) where T = CheapStack(unsafe_convert(Ptr{T}, ptr))

function update_store!(s::CheapStack, ptr)
    s.beginning = ptr
    empty!(s)
end

isempty(s::CheapStack) = s.ptr == typeof(s.ptr)(0)

function empty!(s::CheapStack)
    s.ptr = typeof(s.ptr)(0)
    s
end

first(s::CheapStack) = isempty(s) ? nothing : unsafe_load(s.ptr)

length(s::CheapStack) =
    isempty(s) ? zero(Int) :
    convert(Int8, s.ptr - s.beginning) ÷ (sizeof ∘ eltype)(s.ptr) + 1

function pop!(s::CheapStack)
    isempty(s) && return nothing

    x = unsafe_load(s.ptr)
    if s.ptr == s.beginning
        empty!(s)
    else
        s.ptr -= (sizeof ∘ eltype)(s.ptr)
    end

    x
end

push!(s::CheapStack) = s

function push!(s::CheapStack, x)
    if isempty(s)
        s.ptr = s.beginning
    else
        s.ptr += (sizeof ∘ eltype)(s.ptr)
    end

    unsafe_store!(s.ptr, x)
    s
end

function push!(s::CheapStack, x, y, rest...)
    push!(s, x)
    push!(s, y)

    for z ∈ rest
        push!(s, z)
    end
    s
end

push!(s::CheapStack, x::AbstractArray) = push!(s, x...)

#############
# Iteration #
#############

function iterate(s::CheapStack, state = 1)
    length(s) >= state || return nothing
    unsafe_load(s.ptr, 2 - state), state + 1
end

# @generated eltype(::Type{CheapStack{T}}) where T = T
# @generated eltype(::CheapStack{T}) where T = eltype(CheapStack{T})

###################
# Pretty printing #
###################

show(io::IO, s::CheapStack) = join(io, s, ", ")

function show(io::IO, ::MIME"text/plain", s::CheapStack{T}) where T
    print(io, "CheapStack{$T}:\n", s)
end
