import Base: empty,
             similar,
             ~, &, |, xor, >>>, >>, <<,
             ==,
             show,
             isempty,
             string,
             convert,
             hash,
             zeros, ones,
             firstindex, lastindex

using Random

import Base: bitcount

using LinearAlgebra

using SpecialFunctions: loggamma

export Sequence
## "Efficient storage of marker data."
struct Sequence
    data::BitVector
end

isempty(seq::Sequence) = isempty(seq.data)

@generated empty(::Sequence) = Sequence()

string(sequence::Sequence) = replace(bitstring(sequence.data), r"[ :]" => "")

function show(io::IO, ::MIME"text/plain", seq::Sequence)
    len = length(seq)
    header = "$len-markers Sequence"

    print(io, header)
    isempty(seq) || print(io, ":\n", string(seq))
end

show(io::IO, seq::Sequence) = print(io, string(seq))

"""
    blocksize(::Type{Sequence})
    blocksize(seq)

Size of the blocks of a sequence.
"""
function blocksize end

## Chunks of BitArrays are hardcoded as UInt64
@generated blocksize(::Type{Sequence}) = 8sizeof(UInt64)

@generated blocksize(::Sequence) = blocksize(Sequence)

for (fun, op) ∈ Dict(:chunkidx => :div, :idxinchunk => :mod)
    @eval $fun(::Type{Sequence}, x) = $op(x - 1, blocksize(Sequence)) + 1
    @eval $fun(::Sequence, x) = $fun(Sequence, x)
end

~(sequence::Sequence) = Sequence(broadcast(~, sequence.data))

for fun ∈ [:&, :|, :xor]
    @eval function $fun(sequence1::Sequence, sequence2::Sequence)
        Sequence(broadcast($fun, sequence1.data, sequence2.data))
    end
end

for fun ∈ [:<<, :>>, :>>>]
    @eval $fun(sequence::Sequence, k) = Sequence($fun(sequence.data, k))
end

bitcount(η::Sequence; init::T = 0) where T =
    bitcount(η.data.chunks, init = init)

function hash(η::Sequence, h::UInt)
    h = hash(η.data, h)
    hash(Sequence, h)
end

==(η1::Sequence, η2::Sequence) = η1.data == η2.data
==(seq::Sequence, str::AbstractString) = string(seq) == str
==(str::AbstractString, seq::Sequence) = seq == str

"""
    Sequence(undef, n)
    Sequence(data)
    Sequence([rng = GLOBAL_RNG,], n)
    Sequence([rng = GLOBAL_RNG,] minlength, maxlength)

  - `Sequence(undef, n)`: Uninitialized sequence of length `n`
  - `Sequence(data)`: sequence containing `data`
  - `Sequence(rng, n)`: random sequence of length n
  - `Sequence(rng, minlength, maxlength)`: random sequence of
    length uniformly sampled in `minlength:maxlength`

# Arguments

  - `n`: Number of markers.
  - `data`: Vector containing data.
  - `rng`: random number generator.
  - `minLength, maxLength`: bounds for sequence length.

# Notes

  - To construct an empty sequence, use `empty(Sequence)`
  - To construct a Sequence from a string of bits, use `convert(Sequence, str)`
"""
function Sequence end

Sequence() = Sequence(BitVector())

Sequence(::UndefInitializer, n) = Sequence(BitVector(undef, n))

Sequence(rng::AbstractRNG, n) = Sequence(bitrand(rng, n))

function Sequence(rng::AbstractRNG, minlength, maxlength)
    Sequence(rng, rand(rng, range(minlength, maxlength)))
end

Sequence(n::Integer) = Sequence(GLOBAL_RNG, n)

function Sequence(minlength::Integer, maxlength::Integer)
    Sequence(GLOBAL_RNG, minlength, maxlength)
end

function convert(::Type{Sequence}, str::AbstractString)
    data = zeros(length(str))

    @inbounds for (k, c) ∈ enumerate(str)
        if parse(Bool, c)
            data[k] = true
        end
    end

    Sequence(data)
end

for (f, fill_f) ∈ Dict(:zeros => :falses, :ones => :trues)
    @eval $f(::Type{Sequence}, n::Integer) = (Sequence ∘ $fill_f)(n)
end

#############
# Distances #
#############

export Distance

"""
    abstract type Distance{T}

Distance between two `Sequence`s.

# Implementation
The only required method for a distance `D<:Distance` to be usable for tree
contruction is
```
    distance(::D{T}, ::Sequence, ::Sequence) where T
```
It is also useful to implement a default constructor. For example, if
`D` is a discrete distance,
```
    D() = D{Int}()
```
"""
abstract type Distance{T} end

"""
    distance(Dist, η1, η2)
    distance(Dist, H)

Compute distance between sequences.
"""
function distance end

function distance(Dist::Distance{T}, H::AbstractVector) where {T}
    n = length(H)
    mat = Matrix{T}(undef, n, n)

    @inbounds for j ∈ 1:n
        mat[j, j] = 0
        for i ∈ 1:n
            mat[i, j] = distance(Dist, H[j], H[i])
        end
    end

    mat
end

## Hamming distance
export Hamming
struct Hamming{T} <: Distance{T} end

Hamming() = Hamming{Int64}()

function distance(::Hamming{T}, η1::Sequence, η2::Sequence) where T
    ## Works since unused bits of a BitArray are always set to 0.
    d = zero(T)
    nblocks = (length(η1) - 1) ÷ blocksize(η1) + 1

    @tturbo for k ∈ 1:nblocks
        d += convert(T, count_ones(η1.data.chunks[k] ⊻ η2.data.chunks[k]))
    end

    d
end

## Leftmost marker distance
export Left
struct Left{T} <: Distance{T} end

Left() = Left{Int64}()

distance(::Left{T}, η1::Sequence, η2::Sequence) where {T} =
    convert(T, first(η1) ⊻ first(η2))

#############
# Iteration #
#############

function Base.iterate(seq::Sequence, state = 1)
    state > lastindex(seq) ? nothing : (seq[state], state + 1)
end

@generated Base.eltype(::Type{Sequence}) = Bool

@generated Base.eltype(::Sequence) = eltype(Sequence)

Base.length(seq::Sequence) = length(seq.data)

Base.size(seq::Sequence) = size(seq.data)

## Indexing.
firstindex(seq::Sequence) = firstindex(seq.data)

lastindex(seq::Sequence) = lastindex(seq.data)

Base.getindex(sequence, I...) = getindex(sequence.data, I...)

## Array.
@generated Base.IndexStyle(::Type{Sequence}) = Base.IndexLinear()

Base.similar(seq::Sequence, args...) = similar(seq.data, args...)

Base.setindex!(sequence::Sequence, x, i) = setindex!(sequence.data, x, i)
