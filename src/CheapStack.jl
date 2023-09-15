import Base:
    eltype,
    isempty,
    empty!,
    first,
    length,
    pop!,
    push!,
    show,
    iterate

mutable struct CheapStack{T}
    const store::Vector{T}
    ptr::Int
end

CheapStack(T, n) =
    CheapStack{T}(Vector{T}(undef, n), 0)

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

function push!(s::CheapStack, xs...)
    for x ∈ xs
        push!(s, x)
    end
    s
end

#############
# Iteration #
#############

iterate(s::CheapStack) = isempty(s) ? nothing : (first(s.store), 1)
function iterate(s::CheapStack, state)
    state += 1
    state ≤ s.ptr ? (s.store[state], state) : nothing
end

@generated eltype(::Type{CheapStack{T}}) where T = T
@generated eltype(::CheapStack{T}) where T = eltype(CheapStack{T})

###################
# Pretty printing #
###################

show(io::IO, s::CheapStack) =
    join(io, s, ", ")

show(io::IO, ::MIME"text/plain", s::CheapStack{T}) where T =
    print(io, "CheapStack{$T}:\n", s)
