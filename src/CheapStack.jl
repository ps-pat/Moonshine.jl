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

mutable struct CheapStack{T}
    const store::Vector{T}
    ptr::Int
end

## TODO: Get rid of this constructor.
CheapStack(T, n) = CheapStack{T}(Vector{T}(undef, n), 0)

function CheapStack{T}(n) where T
    CheapStack{T}(Vector{T}(undef, n), 0)
end

isempty(s::CheapStack) = iszero(s.ptr)

empty!(s::CheapStack) = s.ptr = 0

first(s::CheapStack) = s.store[s.ptr]

length(s::CheapStack) = s.ptr

function pop!(s::CheapStack)
    @inline ret = s.store[s.ptr]
    s.ptr -= 1

    ret
end

push!(s::CheapStack) = s

function push!(s::CheapStack, x)
    s.ptr += 1
    @inline s.store[s.ptr] = x

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
    state ≤ s.ptr ? (s.store[state], state) : nothing
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

##############
# Conversion #
##############

vec(s::CheapStack) = resize!(s.store, length(s))
