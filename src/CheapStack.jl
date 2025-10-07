import Base: eltype,
             isempty,
             empty!,
             first,
             length,
             pop!,
             push!,
             show,
             iterate,
             vec

export CheapStack
"""
    $(TYPEDEF)

Simple stack container.

# Functionalities
The following operations are supported:
* `isempty`
* `empty!`
* `first`
* `length`
* `pop!`
* `push!`

For convenience, `CheapStack` implements [the iteration interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-iteration).

# Fields
$(TYPEDFIELDS)
"""
struct CheapStack{T}
    store::UnsafeArray{T, 1}
    ptr::Base.RefValue{Int}
end

"""
    $(SIGNATURES)

Construct a [`CheapStack`](@ref) that uses an `UnsafeArray` as store.
"""
CheapStack(store::UnsafeArray{T, 1}) where T =
    CheapStack{T}(store, Ref{Int}(0))

isempty(s::CheapStack) = iszero(s.ptr[])

empty!(s::CheapStack) = s.ptr = Ref(0)

first(s::CheapStack) = s.store[s.ptr[]]

length(s::CheapStack) = s.ptr[]

function pop!(s::CheapStack)
    @inline ret = s.store[s.ptr[]]
    s.ptr[] -= 1

    ret
end

push!(s::CheapStack) = s

function push!(s::CheapStack, x)
    s.ptr[] += 1

    @inbounds s.store[s.ptr[]] = x

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

#############
# Iteration #
#############

iterate(s::CheapStack) = isempty(s) ? nothing : (first(s.store), 1)
function iterate(s::CheapStack, state)
    state += 1
    state ≤ s.ptr[] ? (s.store[state], state) : nothing
end

@generated eltype(::Type{CheapStack{T}}) where T = T
@generated eltype(::CheapStack{T}) where T = eltype(CheapStack{T})

###################
# Pretty printing #
###################

show(io::IO, s::CheapStack) = join(io, s, ", ")

function show(io::IO, ::MIME"text/plain", s::CheapStack{T}) where T
    print(io, "CheapStack{$T}:\n", s)
end
