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

import Base: bitcount

using LinearAlgebra

export Sequence
## "Efficient storage of marker data."
struct Sequence
    data::BitVector
end

isempty(seq::Sequence) = isempty(seq.data)

@generated empty(::Sequence) = Sequence()

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

#############
# Distances #
#############

export Distance
abstract type Distance{T} end

"""
    distance(Dist, η1, η2)
    distance(Dist, H)

Compute distance between sequences.
"""
function distance end

function distance(Dist::Type{<:Distance{T}}, H::AbstractVector) where T
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

distance(Dist::Type{<:Distance}, H::AbstractVector) = distance(Dist{Float64}, H)

function compute_types(H::AbstractVector{Sequence})
    types = unique(H)
    sort!(types, by = bitcount)
end

export coalescence_matrix
"""
    coalescence_matrix(Dist, H; μ = 1e-5, maximum_distance = ∞, bias0 = ∞)

Compute the transition matrix for the coalescence of sequences types. H is
expected to have been processed by [`compute_types`][@ref] or equivalent.

# Arguments

  - `Dist`: discrete distance
  - `H`: vector of sequences
  - `μ`: per locus mutation rate
  - `maximum_distance`: distance beyond which a type is considered
    inaccessible
  - `bias0`: biasing factor for 0->1 mutations

# Notes
  - Assuming more than one type, μ must be strictly greater than 0 for the
    chain to be recurrent.
  - Each entry of the lower triangle of the distance matrix is multiplied
    by bias0. In other words, the matrix used to compute transition
    probabilities is the Hadamard product of the distance matrix with a lower
    triangular matrix in which each entry is equal to bias0. Using bias0 = b is
    thus equivalent to stretching every 0->1 distances by a factor of b.
  - Relationship to maximum_distance is assessed *after* biasing.
  - For a given source type, if every destination type is beyond
    maximum_distance, the maximum distance for that source-destination pair is
    used instead.

See also [`compute_types`][@ref], [`Distance`][@ref], [`Sequence`][@ref].
"""
function coalescence_matrix end

function coalescence_matrix(Dist::Type{<:Distance{Float64}}, H::AbstractVector;
                            μ = 1e-5, maximum_distance = 1, bias0 = ∞)
    ## Compute the distance matrix and apply bias.
    distance_matrix = distance(Dist, H)

    for σ ∈ axes(distance_matrix, 2)[2:end]
        for δ ∈ (σ + 1):size(distance_matrix, 1)
            distance_matrix[δ, σ] *= bias0
        end
    end

    transition_matrix = similar(distance_matrix, Float64)

    for (σ, distances) ∈ (enumerate ∘ eachcol)(distance_matrix)
        vec_norm = 0

        ## If every other types is at a distance greater than the
        ## maximum distance, the type of the mrca of the sample is not
        ## well defined.
        maximum_distance_σ = max(maximum_distance,
                                 minimum(d -> iszero(d) ? typemax(d) : d,
                                         distance_matrix[:,σ]))

        for (δ, distance) ∈ enumerate(distances)
            ## 1:      No self-loop;
            ## (2, 3): A state at a distance greater than the maximum
            ##         allowed distance is not accessible.
            if σ == δ || distance == typemax(Int) || distance > maximum_distance_σ
                transition_matrix[δ, σ] = 0
                continue
            end

            d = distances[δ]
            prob = exp(d * log(μ) - sum(log, 2:d, init = zero(Float64)))
            transition_matrix[δ, σ] = prob
            vec_norm += transition_matrix[δ, σ]
        end

        ## If the norm of the vector is zero, we create an absorbing
        ## state.
        if iszero(vec_norm)
            transition_matrix[σ, σ] = 1
        else
            ## Otherwise, normalize the vector.
            transition_matrix[:, σ] ./= vec_norm
        end
    end

    transition_matrix
end

function coalescence_matrix(Dist::Type{<:Distance}, H::AbstractVector;
                            μ = 1e-7, maximum_distance = 1, bias0 = ∞)
    coalescence_matrix(Dist{Float64}, H,
                       μ = μ,
                       maximum_distance = maximum_distance,
                       bias0 = bias0)
end

## Hamming distance
export Hamming
struct Hamming{T} <: Distance{T} end

function distance(::Type{Hamming{T}}, η1::Sequence, η2::Sequence) where T
    ## Works since unused bits of a BitArray are always set to 0.
    d = zero(T)
    nblocks = (length(η1) - 1) ÷ blocksize(η1) + 1

    @inbounds for k ∈ 1:nblocks
        d += count_ones(η1.data.chunks[k] ⊻ η2.data.chunks[k])
    end

    d
end

function distance(::Type{Hamming}, η1::Sequence, η2::Sequence)
    distance(Hamming{Int}, η1, η2)
end

distance(::Type{Hamming}, H::AbstractVector) = distance(Hamming{Int}, H)

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
Base.firstindex(seq::Sequence) = firstindex(seq.data)

Base.lastindex(seq::Sequence) = lastindex(seq.data)

Base.getindex(sequence, I...) = getindex(sequence.data, I...)

## Array.
@generated Base.IndexStyle(::Type{Sequence}) = Base.IndexLinear()

Base.similar(seq::Sequence, args...) = similar(seq.data, args...)

Base.setindex!(sequence::Sequence, x, i) = setindex!(sequence.data, x, i)
