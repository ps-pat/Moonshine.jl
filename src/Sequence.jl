import Base: empty,
             similar,
             ~, &, |, xor, >>>, >>, <<,
             ==,
             show,
             isempty,
             string,
             convert,
             hash

using Random

export Sequence
## "Efficient storage of marker data."
struct Sequence
    data::BitVector
end

isempty(seq::Sequence) = isempty(seq.data)

@generated empty(::Sequence) = Sequence(BitVector())

string(sequence::Sequence) = filter(!isspace, bitstring(sequence.data))

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

~(sequence::Sequence) = Sequence(broadcast(~, sequence.data))

for fun ∈ [:&, :|, :xor]
    @eval $fun(sequence1::Sequence, sequence2::Sequence) =
        Sequence(broadcast($fun, sequence1.data, sequence2.data))
end

for fun ∈ [:<<, :>>, :>>>]
    @eval $fun(sequence, k) = Sequence($fun(sequence, k))
end

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

Sequence(::UndefInitializer, n) = Sequence(BitVector(undef, n))

Sequence(rng::AbstractRNG, n) = Sequence(bitrand(n))

Sequence(rng::AbstractRNG, minlength, maxlength) =
    Sequence(rng, rand(rng, range(minlength, maxlength)))

Sequence(n::Integer) = Sequence(GLOBAL_RNG, n)

Sequence(minlength::Integer, maxlength::Integer) =
    Sequence(GLOBAL_RNG, minlength, maxlength)

function convert(::Type{Sequence}, str::AbstractString)
    data = zeros(length(str))

    @inbounds for (k, c) ∈ enumerate(str)
        if parse(Bool, c)
            data[k] = true
        end
    end

    Sequence(data)
end

## Iteration.
function Base.iterate(seq::Sequence, state = 1)
    state > lastindex(seq) ? nothing : (seq[state], state + 1)
end

@generated Base.eltype(::Type{Sequence}) = Bool

@generated Base.eltype(::Sequence) = eltype(Sequence)

Base.length(seq::Sequence) = length(seq.data)

Base.size(seq::Sequence) = size(seq.data)

## Indexing.
Base.firstindex(seq::Sequence) = firstindex(seq.data)

Base.lastindex(seq::Sequence) = lastindex(seq.data)

Base.getindex(sequence, I...) = getindex(sequence.data, I...)

## Array.
@generated Base.IndexStyle(::Type{Sequence}) = Base.IndexLinear()

Base.similar(seq::Sequence, args...) = similar(seq.data, args...)

Base.setindex!(sequence::Sequence, x, i) = setindex!(sequence.data, x, i)
