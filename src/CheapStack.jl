import Base:
    eltype,
    isempty,
    empty!,
    first,
    length,
    pop!,
    push!

mutable struct CheapStack{T}
    const store::Vector{T}
    ptr::Int
end

CheapStack(T, n) =
    CheapStack{T}(Vector{T}(undef, n), 0)

eltype(::Type{CheapStack{T}}) where T = T

isempty(s::CheapStack) = iszero(s.ptr)

empty!(s::CheapStack) = s.ptr = 0

first(s::CheapStack) = isempty(s) ? nothing : s.store[s.ptr]

length(s::CheapStack) = s.ptr

function pop!(s::CheapStack)
    isempty(s) && return nothing

    ret = s.store[s.ptr]
    s.ptr -= 1

    ret
end

function push!(s::CheapStack, x)
    s.ptr += 1
    s.store[s.ptr] = x

    s
end
