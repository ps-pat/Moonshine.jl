import Base:
    empty,
    similar,
    ~, &, |, xor, >>>, <<,
    ==,
    show,
    isempty,
    string,
    convert

using Random: GLOBAL_RNG, AbstractRNG

"Efficient storage of marker data."
struct Sequence{T <: Unsigned}
    data::Vector{T}
    n::Int
end

isempty(seq::Sequence) = iszero(length(seq))

@generated empty(::Sequence{T}) where T = Sequence(T[], 0)

function string(seq::Sequence)
    isempty(seq) && return ""

    strings = bitstring.(seq.data)
    strings[1] = strings[1][posinblock(seq, 1):end]

    prod(strings)
end

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

@generated blocksize(::Type{Sequence{T}}) where T = 8sizeof(T)

@generated blocksize(::Sequence{T}) where T = blocksize(Sequence{T})

"""
    offsetidx(seq, idx)
    offsetidx(n, blocksize, idx)

Offset the index of a marker to account for the fact that the number
of markers is not necessarily a multiple of the size of a block.

"""
function offsetidx end

function offsetidx(n, blocksize, idx)
    ## Number of markers stored in the last (leftmost) block.
    lbas = n % blocksize

    ## Account for the offset created by a `nmarker` that is not a
    ## multiple of `blocksize`. Also `-1` because idx starts at 1...
    idx + (blocksize - lbas) % blocksize - 1
end

"""
    blockidx(seq, idx)
    blockidx(n, blocksize, idx)

Return the index (starting at 1) of the block containing marker `idx`
(starting at 1 as well).
"""
function blockidx end

"""
    posinblock(seq, idx)
    posinblock(n, blocksize, idx)

Return the index (starting at 1) of the `idx`th marker (starting at 1)
inside its block.
"""
function posinblock end

for (fun, op) ∈ Dict(:blockidx => :÷, :posinblock => :%)
    @eval function $fun(n, blocksize, idx)
        $op(offsetidx(n, blocksize, idx), blocksize) + 1
    end
end

"""
    blocklength(seq, idx)
    blocklength(n, blocksize, idx)

Return the actual number of markers stored in a block. `idx` is the
index of the marker of interest, not of the block itself.
"""
blocklength(n, blocksize, idx) =
    blocksize - (blockidx(n, blocksize, idx) <= 1 ?
    offsetidx(n, blocksize, 1) : zero(Int))


## "Sequence" versions
for fun ∈ [:offsetidx,
           :blockidx,
           :posinblock,
           :blocklength]
    @eval function $fun(seq::Sequence, idx)
        $fun(length(seq), blocksize(seq), idx)
    end
end

"""
    nblocks(n, blocksize)

Numberof blocks of size `blocksize` needed to store `n` markers.
"""
nblocks(n, blocksize) = div(n, blocksize, RoundUp)

function ~(seq::Sequence)
    isempty(seq) && return deepcopy(seq)

    data = (~).(seq.data)
    data[1] &= typemax(data[1]) >>> offsetidx(seq, 1)

    Sequence(data, length(seq))
end

for fun ∈ [:&, :|, :xor]
    @eval function $fun(seq1::Sequence, seq2::Sequence)
        data = $fun.(seq1.data, seq2.data)
        n = seq1.n

        Sequence(data, n)
    end
end

function <<(seq::Sequence{T}, shift) where T
    iszero(shift) && return deepcopy(seq)
    shift < 0 && return >>>(seq, -shift)
    shift >= length(seq) && return Sequence(false, length(seq))

    newseq = Sequence{T}(undef, length(seq))
    n = length(newseq.data)
    bs = blocksize(seq)

    intershift = shift ÷ bs
    intrashift = shift % bs

    @inbounds @simd for k ∈ range(1, n - intershift - 1)
        newseq.data[k] = (seq.data[k + intershift] << intrashift) |
            (seq.data[k + intershift + 1] << (intrashift - bs))
    end

    newseq.data[n - intershift] = (seq.data[n] << intrashift)

    @inbounds @simd for k ∈ range(n - intershift + 1, n)
        newseq.data[k] = 0
    end

    newseq
end

function >>>(seq::Sequence{T}, shift) where T
    iszero(shift) && return deepcopy(seq)
    shift < 0 && return <<(seq, -shift)
    shift >= length(seq) && return Sequence(false, length(seq))

    newseq = Sequence{T}(undef, length(seq))
    n = length(newseq.data)
    bs = blocksize(seq)

    intershift = shift ÷ bs
    intrashift = shift % bs

    @inbounds @simd for k ∈ 1:intershift
        newseq.data[k] = 0
    end

    newseq.data[intershift + 1] = seq.data[1] >>> intrashift

    @inbounds @simd for k ∈ range(intershift + 2, n)
        newseq.data[k] = (seq.data[k - intershift] >>> intrashift) |
            (seq.data[k - intershift - 1] >>> (intrashift - bs))
    end

    newseq
end

==(seq1::Sequence, seq2::Sequence) =
    (seq1.n == seq2.n) && (seq1.data == seq2.data)
==(seq::Sequence, str::AbstractString) = string(seq) == str
==(str::AbstractString, seq::Sequence) = seq == str

"""
    Sequence()
    Sequence(undef, n)
    Sequence(data, n)
    Sequence(rng = GLOBAL_RNG, minlength, maxlength = 0)

- `Sequence{T <: Unsigned}()`: empty sequence.
- `Sequence{T <: Unsigned}(undef, n)`: Uninitialized sequence of length `n`.
- `Sequence{T <: Unsigned}(data, n)`: sequence of `n` markers containing `data`.
- `Sequence{T <: Unsigned}(rng, minlength, maxlength)`: random sequence of
  length uniformly sampled in `minlength:maxlength`. If `iszero(maxlength)`,
  `maxlength = minlength`.

# Arguments
- `T`: DataType used to store data. If unspecified, defaults to `UInt`.
- `n`: Number of markers.
- `data`: Vector containing data.
- `rng::Random.AbstractRNG`: random number generator.
- `minLength, maxLength`: bounds for sequence length.
"""
function Sequence end

@generated Sequence{T}() where T = Sequence(T[], 0)
@generated Sequence() = Sequence{UInt}()

Sequence{T}(::UndefInitializer, n) where T <: Unsigned =
    Sequence(Vector{T}(undef, nblocks(n, 8sizeof(T))), n)
Sequence(::UndefInitializer, n) = Sequence{UInt}(undef, n)

function Sequence{T}(rng::AbstractRNG,
                     minlength::Integer,
                     maxlength::Integer = 0) where T <: Unsigned
    blocksize = 8sizeof(T)
    n = iszero(maxlength) ? minlength : rand(rng, minlength:maxlength)
    nb = nblocks(n, blocksize)
    actual_n = blocksize * nb

    data = rand(rng, T, nb)
    Sequence{T}((Sequence{T}(data, actual_n) >>> (actual_n - n)).data, n)
end
Sequence(rng::AbstractRNG, minlength::Integer, maxlength::Integer = 0) =
    Sequence{UInt}(rng, minlength, maxlength)


Sequence{T}(minlength::Integer, maxlength::Integer = 0) where T <: Unsigned =
    Sequence{T}(GLOBAL_RNG, minlength, maxlength)
Sequence(minlength::Integer, maxlength::Integer = 0) =
    Sequence{UInt}(minlength, maxlength)

function convert(::Type{Sequence{T}}, str::AbstractString) where T
    n = length(str)
    bs = 8sizeof(T)
    nb = nblocks(n, bs)

    data::Vector{UInt} =
        [chop(str, tail = Int(block .* 64)) |>
         σ -> last(σ, 64) |>
         σ -> parse(T, σ, base = 2)
         for block ∈ range(nb - 1, 0, step = -1)]

    n = isempty(data) ? 0 : n

    Sequence(data, n)
end

convert(::Type{Sequence}, str::AbstractString) = convert(Sequence{UInt}, str)

"""
    fillseq(derived, n, blocktype = UInt)

Create a sequence of length `n` filled exclusively with 0s or 1s.
"""
function fillseq(derived, n, blocktype = UInt)
    iszero(n) && return Sequence{blocktype}([], 0)

    blocksize = 8sizeof(blocktype)
    nb = nblocks(n, blocksize)
    data = derived ? fill(typemax(blocktype), nb) : zeros(blocktype, nb)
    derived && (data[1] >>>= nb * blocksize - n)

    Sequence(data, n)
end

"""
    andmask(seq, unmasked)
    andmask(n, unmasked)
    andmask(T, unmasked)

Construct an "and" mask for a sequence.
"""
function andmask end

function andmask(n::Integer, unmasked::UnitRange{Int})
    isempty(unmasked) && return fillseq(false, n)

    rshift = first(unmasked) - 1
    lshift = n - last(unmasked)

    seq = fillseq(true, n)

    (seq << lshift) & (seq >>> rshift)
end

## TODO: Implement for StepRange.

function andmask(n::Integer, unmasked::Integer)
    mask = fillseq(false, n)
    mask[unmasked] = true
    mask
end

## Sequence versions.
for par ∈ [:(UnitRange{<:Integer}), :Integer]
    @eval function andmask(seq::Sequence, unmasked::$par)
        andmask(length(seq), unmasked)
    end
end

andmask(T::Type{<:Unsigned}, unmasked::Integer) =
    one(T) << (8sizeof(T) - unmasked)

"""
    derivedpos(seq)
    ancestralpos(seq)

Derived/ancestral positions of a sequence.
"""
function derivedpos end,
function ancestralpos end

derivedpos(seq::Sequence) = findall(collect(seq))

ancestralpos(seq::Sequence) =
    setdiff(range(1, length = length(seq)), derivedpos(seq))

## Iteration.
Base.iterate(seq::Sequence, state = 1) =
    state > lastindex(seq) ? nothing : (seq[state], state + 1)

@generated Base.eltype(::Type{Sequence}) = Bool

@generated Base.eltype(::Sequence) = eltype(Sequence)

Base.length(seq::Sequence) = seq.n

Base.size(seq::Sequence) = (Base.length(seq),)

## Indexing.
Base.firstindex(seq::Sequence) = length(seq) > 0 ? 1 : 0

Base.lastindex(seq::Sequence) = length(seq)

Base.getindex(seq::Sequence{T}, i) where T =
    seq.data[blockidx(seq, i)] & andmask(T, posinblock(seq, i)) > 0

Base.getindex(seq::Sequence, rng::UnitRange) =
    map(idx -> getindex(seq, idx), rng)

## Array.
@generated Base.IndexStyle(::Type{Sequence}) = Base.IndexLinear()

Base.similar(seq::Sequence{T}, n = length(seq)) where T = Sequence(undef, n, T)

function Base.setindex!(seq::Sequence{T}, x, i) where T
    x_bool = convert(Bool, x)

    mask = andmask(T, posinblock(seq, i))
    _blockidx = blockidx(seq, i)
    if x_bool
        seq.data[_blockidx] |= mask
    else
        seq.data[_blockidx] &= ~mask
    end
end
