import Base: empty,
             similar,
             ~, &, |, xor, >>>, >>, <<, *,
             ==, isequal,
             show, display,
             isempty,
             string,
             convert,
             hash,
             zeros, ones,
             firstindex, lastindex,
             copy!,
             sum

using Random

using LinearAlgebra

using SpecialFunctions: loggamma

using UnicodePlots: heatmap, label!, annotate!

using ChunkSplitters: index_chunks

export Sequence
"""
    $(TYPEDEF)

Sequence of biallelic genetic markers (haplotype).

Implement [the iteration interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-iteration) as well as standard bitwise operations.

!!! info
    Random constructors sample marker's states via [`Random.bitrand`](@extref).

# Fields
$(TYPEDFIELDS)

# Constructors
$(METHODLIST)

where:
* `n`: number of markers
* `rng`: random number generator
* `minLength, maxLength`: bounds for sequence length
"""
struct Sequence
    data::BitVector
end

#          +----------------------------------------------------------+
#          |                         Equality                         |
#          +----------------------------------------------------------+

function hash(s::Sequence, h::UInt)
    h = hash(s.data.chunks, h)
    hash(:Sequence, h)
end

(==)(h1::Sequence, h2::Sequence) = h1.data == h2.data

# isequal(h1::Sequence, h2::Sequence) = isequal(h1.data, h2.data)

function isequal(h1::Sequence, h2::Sequence,
                 m1::AbstractVector{UInt64}, m2::AbstractVector{UInt64})
    @inbounds for k ∈ eachindex(m1)
        mask = m1[k] & m2[k]
        iszero((h1.data.chunks[k] ⊻ h2.data.chunks[k]) & mask) || return false
    end

    true
end

isempty(seq::Sequence) = isempty(seq.data)

"""
    $(FUNCTIONNAME)

Simple (and fast) hash function.

--*Internal*--
"""
function cheap_hash(s::Sequence)
    h = 0x5aacc77a902f6f2e
    o = 0x5246e7da0d562ddf
    ptr = unsafe_convert(Ptr{UInt32}, pointer(s.data.chunks))

    @inbounds for k ∈ eachindex(s.data.chunks)
        o += 0x5851f42d4c957f20
        c = unsafe_load(ptr, 2k - 1)
        x = c + o

        o += 0x14057b7ef7678140
        c = unsafe_load(ptr, 2k)
        x *= c + o

        h ⊻= x
    end

    h
end

@generated empty(::Sequence) = Sequence()

string(sequence::Sequence) = replace(bitstring(sequence.data), r"[ :]" => "")

show(io::IO, ::MIME"text/plain", h::Sequence) =
    print(io, plot_sequence(h, height = 2, colorbar = false))

show(io::IO, h::Sequence) =
    print(io, plot_sequence(h, height = 2, colorbar = false, title = "", labels = false))

"""
    $(FUNCTIONNAME)(::Type{Sequence})
    $(FUNCTIONNAME)(::Sequence)

Size (in *bits*) of the blocks (chunks) of a sequence.

# Methods
$(METHODLIST)

--*Internal*--
"""
function blocksize end

## Chunks of BitArrays are hardcoded as UInt64
@generated blocksize(::Type{Sequence}) = 8sizeof(UInt64)

@generated blocksize(::Sequence) = blocksize(Sequence)

"""
    $(FUNCTIONNAME)(n, x)
    $(FUNCTIONNAME)(::Type{Sequence}, x)
    $(FUNCTIONNAME)(::Sequence, x)

Index of the chunk associated with marker `x`.

See also [`idxinchunk`](@ref).

--*Internal*--
"""
function chunkidx end

"""
    $(FUNCTIONNAME)(n, x)
    $(FUNCTIONNAME)(::Type{Sequence}, x)
    $(FUNCTIONNAME)(::Sequence, x)

Index of marker `x` in its associated chunk.

See also [`chunkidx`](@ref).

--*Internal*--
"""
function idxinchunk end

function chunkidx(n, x)
    d = div(x, n)
    d + (d * n != x)
end

idxinchunk(n, x) = mod1(x, n)

for (fun, op) ∈ Dict(:chunkidx => :fld1, :idxinchunk => :mod1)
    @eval begin
        $fun(::Type{Sequence}, x) = $fun(blocksize(Sequence), x)
        $fun(::Sequence, x) = $fun(Sequence, x)
    end
end

~(sequence::Sequence) = Sequence(broadcast(~, sequence.data))

for fun ∈ [:&, :|, :xor]
    @eval function ($fun)(h1::Sequence, h2)
        h = Sequence(undef, length(h1))
        @inbounds @simd for k ∈ eachindex(h.data.chunks)
            h.data.chunks[k] = ($fun)(h1.data.chunks[k], h2[k])
        end

        h
    end

    @eval ($fun)(h1, h2::Sequence) = ($fun)(h2, h1)

    @eval ($fun)(h1::Sequence, h2::Sequence) = ($fun)(h1, h2.data.chunks)

end

for fun ∈ [:<<, :>>, :>>>]
    @eval $fun(sequence::Sequence, k::Integer) = Sequence($fun(sequence.data, k))
end

function *(h::Sequence, v::AbstractVector{T}) where T
    acc = zero(T)
    k = 1
    @inbounds @simd for chunk ∈ h.data.chunks
        Δ = trailing_zeros(chunk)
        k += Δ
        while Δ < blocksize(h)
            acc += v[k]
            chunk >>>= Δ + 1

            Δ = trailing_zeros(chunk)
            k += Δ
        end
    end

    acc
end

*(v::AbstractVector, h::Sequence) = h * v

Sequence() = Sequence(BitVector())

Sequence(undef::UndefInitializer, n) = Sequence(BitVector(undef, n))

Sequence(rng::AbstractRNG, n) = Sequence(bitrand(rng, n))

function Sequence(rng::AbstractRNG, minlength, maxlength)
    Sequence(rng, rand(rng, range(minlength, maxlength)))
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

"""
    $(SIGNATURES)

Set all markers to wild state (0).

--*Internal*--
"""
function wipe!(h::Sequence)
    h.data.chunks .⊻= h.data.chunks
    h
end

function copy!(h1::Sequence, h2::Sequence)
    h1.data .⊻= h1.data .⊻ h2.data
    h1
end

function sum(h::Sequence)
    ret = zero(Int)
    @inbounds @simd for chunk ∈ h.data.chunks
        ret += count_ones(chunk)
    end

    ret
end

"""
    $(FUNCTIONNAME)(mask, h, idx; wipe = true)
    $(FUNCTIONNAME)(mask, sample, ω; wipe = true)
    $(FUNCTIONNAME)(mask, genealogy, x; wipe = true)

`&`-mask for non-ancestral markers.
Construct a mask for a range of markers in a given haplotype.

The first form returns a mask that `&` markers outside of the range `idx` to 0
when applied to `h.data.chunks`. `h` can either be a sequence or an integer
meant to be the number of markers in a sequence. The second forms acts
similarly for an interval (or set of intervals). If `x` is a vertex, the third
form returns a mask that `&` the markers of the associated haplotype not
ancestral for it to 0. The ancestral intervals of `x` is defined as the union
of its downstream edges' intervals. Behaviour is similar if `x` is an edge, but
the mask is meant for the intersection of chunks between its adjacent vertices.

If `wipe = true`, `mask` is wiped beforehand.

For convenience, each method has an allocating counterpart; see
[`ancestral_mask`](@ref).

See also [`wipe!`](@ref).

# Methods
$(METHODLIST)

--*Internal*--
"""
function ancestral_mask! end

function ancestral_mask!(mask, h::Int, idx; wipe = true)
    bs = blocksize(Sequence)
    nchunks = div(h, bs, RoundUp)

    wipe && wipe!(mask)

    chunks_idx = chunkidx(Sequence, idx.start), chunkidx(Sequence, idx.stop)

    first(chunks_idx) > last(chunks_idx) && return mask

    ## Evrything is in reverse in line with the way data is stored in
    ## BitVectors.

    ## First chunk
    ormask = typemax(UInt)
    let nleading0 = idxinchunk(Sequence, idx.start) - 1,
        ntrailing0 = bs - (first(chunks_idx) == last(chunks_idx) ?
            idxinchunk(Sequence, idx.stop) : bs)

        ormask >>>= ntrailing0
        ormask &= typemax(UInt) << nleading0

        mask[first(chunks_idx)] |= ormask
    end

    ## Early termination in case only one chunk needs to be modified
    first(chunks_idx) == last(chunks_idx) && return mask

    ## Middle chunks
    let nchunks = last(chunks_idx) - first(chunks_idx) - 1,
        extra_chunks = nchunks % simd_chunksize,
        lane = VecRange{simd_chunksize}(0)

        for k ∈ range(first(chunks_idx) + 1, last(chunks_idx) - 1 - extra_chunks, step = simd_chunksize)
            mask[lane + k] |= typemax(UInt)
        end

        for k ∈ range(last(chunks_idx) - extra_chunks, last(chunks_idx) - 1)
            mask[k] |= typemax(UInt)
        end
    end

    ## Last chunk
    ormask = typemax(UInt)
    let ntrailing0 = bs - idxinchunk(Sequence, idx.stop)
        ormask >>>= ntrailing0

        mask[last(chunks_idx)] |= ormask
    end

    mask
end

ancestral_mask!(mask, h::Sequence, idx; wipe = true) =
    ancestral_mask!(mask, length(h), idx, wipe = wipe)

"""
Allocating version of [`ancestral_mask!`](@ref)

# Methods
$(METHODLIST)

--*Internal*--
"""
function ancestral_mask end

ancestral_mask(h::Sequence, idx) =
    ancestral_mask!(similar(h.data.chunks), h, idx)

ancestral_mask(h::Int, idx) =
    ancestral_mask!(BitVector(undef, h), h, idx)

#############
# Distances #
#############

export Distance
"""
    $(TYPEDEF)

Distance between two [`Sequence`](@ref)s.

# Implementation
The only required method for a distance `D<:Distance` to be usable for tree
contruction is

    distance(::D{T}, ::Sequence, ::Sequence) where T

It is also useful to implement a default constructor. For example, if
`D` is a discrete distance,

    D() = D{Int}()
"""
abstract type Distance{T} end

export distance
"""
    $(FUNCTIONNAME)(Dist, h1, h2)
    $(FUNCTIONNAME)(Dist, H)

Compute distances between sequences.

# Methods
$(METHODLIST)
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

export Hamming
"""
    $(TYPEDEF)

Hamming distance.
"""
struct Hamming{T} <: Distance{T} end

Hamming() = Hamming{Int64}()

function distance(::Hamming{T}, h1::Sequence, h2::Sequence) where T
    d = zero(T)
    nblocks = length(h1.data.chunks)

    @inbounds @simd for k ∈ 1:(nblocks-1)
        d += count_ones(h1.data.chunks[k] ⊻ h2.data.chunks[k])
    end

    lastchunk = last(h1.data.chunks) ⊻ last(h2.data.chunks)
    n_not0 = blocksize(h1) * nblocks - length(h1)
    lastchunk &= 0xffffffffffffffff >> n_not0
    d += count_ones(lastchunk)

    d
end

export LeftM
"""
    $(TYPEDEF)

Left marker distance.

This is only the discrete metric for the leftmost marker. Technically not a
metric for haplotypes, but widely used.
"""
struct LeftM{T} <: Distance{T} end

LeftM() = LeftM{Int64}()

distance(::LeftM{T}, h1::Sequence, h2::Sequence) where {T} =
    convert(T, first(h1) ⊻ first(h2))

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

Base.getindex(h::Sequence, i::Integer) = getindex(h.data, i)
Base.getindex(h::Sequence, I...) = Sequence(getindex(h.data, I...))

## Array.
@generated Base.IndexStyle(::Type{Sequence}) = Base.IndexLinear()

Base.similar(seq::Sequence, args...) = similar(seq.data, args...)

Base.setindex!(sequence::Sequence, x, i) = setindex!(sequence.data, x, i)

#          +----------------------------------------------------------+
#          |                         Plotting                         |
#          +----------------------------------------------------------+

export plot_sequence
"""
    $(SIGNATURES)

Unicode-graphic representation of an haplotype.

# Arguments
* `nbins` (`clamp(length(h), 1, 69)`): number of bins
* ̀`height` (`7`): number of rows
* `kwargs`: arguments passed directly to `UnicodePlots.heatmap`
"""
function plot_sequence(h::Sequence;
                       nbins = clamp(length(h), 1, 69),
                       height = 7,
                       kwargs...)
    counts_vec = map(idx -> count(h[idx]), index_chunks(1:length(h), n = nbins))
    counts = repeat(reshape(counts_vec, (1, length(counts_vec))), height)

    plt = heatmap(counts,
                  width = nbins,
                  border = :none,
                  margin = 0,
                  title = (string ∘ length)(h) * "-markers Sequence",
                  height = height,
                  colorbar = true,
                  zlabel = "#1",
                  colormap = default_colormap;
                  kwargs...)

    for k ∈ 1:height
        label!(plt, :l, k, "")
    end

    label!(plt, :bl, string(1), color = :white)
    label!(plt, :br, (string ∘ length)(h), color = :white)
    label!(plt, :b, string(length(h) ÷ 2), color = :white)

    plt
end
