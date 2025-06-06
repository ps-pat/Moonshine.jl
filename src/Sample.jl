import Base: iterate,
             eltype,
             length,
             size,
             getindex,
             IndexStyle

using Graphs: Edge

export Sample
"""
    $(TYPEDEF)

Contain a sample of haplotypes and informations about them.

Implements [the iteration interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-iteration) and [the array interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array).

# Fields
$(TYPEDFIELDS)

# Constructors
!!! info
    Random constructor sample sequences via [msprime](https://tskit.dev/msprime/docs/stable/intro.html) using a [binary mutation model](https://tskit.dev/msprime/docs/stable/api.html#msprime.BinaryMutationModel).

$(METHODLIST)

where:
* `n`: number of sequences
* `rng`: random number generator
* `ts`: [Tree Sequence](https://tskit.dev/tskit/docs/latest/python-api.html#the-treesequence-class)
"""
struct Sample <: AbstractVector{Sequence}
    "Vector of haplotypes"
    H::Vector{Sequence}
    "Unscaled (per-locus) mutation rate"
    μ::Float64
    "Unscaled (per-locus) recombination rate"
    ρ::Float64
    "Effective population size"
    Ne::Float64
    "Sequence length"
    sequence_length::Float64
    "Marker's positions"
    positions::Vector{Float64}
    "Coefficients for positions' line"
    coefs::NTuple{2, Float64}
end

function Sample(H::AbstractVector{Sequence};
                μ = 1e-8, ρ = 0, Ne = 10_000,
                sequence_length = length(H) + 1,
                positions = 1:length(H))
    ## Compute coefficients.
    n = length(positions)
    sum_positions = sum(positions)

    m = (12 * (1:n)' * positions) / ((n - 1) * n * (n + 1)) -
        (6 * sum_positions) / ((n - 1) * n)
    b = sum_positions / n - (n + 1) * m / 2

    Sample(H, μ, ρ, Ne, sequence_length, positions, (b, m))
end

function Sample(ts::TreeSequence)
    n = zero(Int)
    μ = zero(Float64)
    ρ = zero(Float64)
    Ne = zero(Float64)
    sequence_length = one(Float64)

    ## Retreive relevant genetic parameters.
    for p ∈ provenances(ts)
        parameters = p.record[:parameters]

        n = get(parameters, :samples, n)
        μ = get(parameters, :rate, μ)
        ρ = get(parameters, :recombination_rate, ρ)
        Ne = get(parameters, :population_size, Ne)
        sequence_length = get(parameters, :sequence_length, sequence_length)
    end

    ## Retreive positions and reconstruct sequence.
    nmarkers = pyconvert(Int, ts.obj.num_mutations)
    positions = Vector{Float64}(undef, nmarkers)
    H = [zeros(Sequence, nmarkers) for _ ∈ 1:n]

    @inbounds for (mutation, variant) ∈ enumerate(variants(ts))
        pos = pyconvert(Float64, variant.site.position)
        positions[mutation] = pos

        genotypes = pyconvert(Vector{Int}, variant.genotypes)
        for (sam, marker) ∈ enumerate(genotypes)
            iszero(marker) && continue
            H[sam][mutation] = true
        end
    end

    Sample(H, μ = μ, ρ = ρ, Ne = Ne,
           sequence_length = sequence_length, positions = positions)
end

function Sample(rng::AbstractRNG, n, μ, ρ, Ne, sequence_length)
    sim_ancestry = msprime[].sim_ancestry
    sim_mutations = msprime[].sim_mutations
    BinaryMutationModel = msprime[].BinaryMutationModel

    seed = rand(rng, UInt32)
    seed_mutations = rand(rng, UInt32)

    ts = sim_ancestry(random_seed = seed,
                      samples = n, ploidy = 1,
                      sequence_length = sequence_length,
                      recombination_rate = ρ,
                      population_size = Ne,
                      record_provenance = true)

    mutated_ts = sim_mutations(ts,
                               random_seed = seed_mutations,
                               rate = μ,
                               discrete_genome = false,
                               model = BinaryMutationModel(false),
                               record_provenance = true)

    Sample(pyconvert(TreeSequence, mutated_ts))
end

###################
# Pretty printing #
###################

function show(io::IO, m::MIME"text/plain", sample::Sample)
    invoke(show, Tuple{IO, MIME"text/plain", AbstractVector}, io, m, sample)

    print(io, "\nlength = " * string(sample.sequence_length) *
              ", μ = " * string(sample.μ) *
              ", ρ = " * string(sample.ρ) *
              ", Ne = " * string(sample.Ne))
end

##################################
# Array and Iteration Interfaces #
##################################

for met ∈ (:length, :size, :firstindex, :lastindex)
    @eval $met(sample::Sample) = $met(sample.H)
end

iterate(sample::Sample, state = 1) =
    state > lastindex(sample) ? nothing : (sample[state], state + 1)

@generated eltype(::Type{Sample}) = Sequence

@generated eltype(::Sample) = eltype(Sample)

getindex(sample::Sample, I...) = getindex(sample.H, I...)

@generated IndexStyle(::Type{Sample}) = Base.IndexLinear()

#################
# Short Methods #
#################

positions(sample::Sample) = sample.positions

nmarkers(sample::Sample) = (length ∘ first)(sample.H)

nmarkers(sample::Sample, ω::Ω) = length(postoidx(sample, ω))

idxtopos(sample::Sample, idx) =
    iszero(idx) ? zero(Float64) : getindex(positions(sample), idx)

function _postoidx_approx(sample, pos)
    idx = floor(Int, (pos - first(sample.coefs)) / last(sample.coefs))
    clamp(idx, 1, nmarkers(sample))
end

function postoidx(sample::Sample, pos::Float64)
    lidx = (firstindex ∘ positions)(sample)
    ridx = (lastindex ∘ positions)(sample)

    pos <= (first ∘ positions)(sample) && return lidx
    pos >= (last ∘ positions)(sample) && return ridx

    ## Initial guess
    idx = _postoidx_approx(sample, pos)

    ## One iteration of binary search
    if idxtopos(sample, idx)< pos
        lidx = idx
    elseif idxtopos(sample, idx) > pos
        ridx = idx
    else
        return idx
    end

    ridx <= lidx + 1 && return ridx

    ## Interpolation-sequential search
    idx = lidx + floor(Int, ((pos - idxtopos(sample, lidx)) * (ridx - lidx)) /
                (idxtopos(sample, ridx) - idxtopos(sample, lidx)))
    if idxtopos(sample, idx) < pos
        idx += 1
        @inbounds while idxtopos(sample, idx) < pos
            idx += 1
        end
    elseif idxtopos(sample, idx) > pos
        idx -= 1
       @inbounds while idxtopos(sample, idx) > pos
            idx -= 1
        end

        if idxtopos(sample, idx) != pos
            idx += 1
        end
    end

    idx
end

function postoidx(sample::Sample, ω::Ω)
    lpos, rpos = endpoints(ω)

    lidx = postoidx(sample, lpos)
    ridx = postoidx(sample, rpos)

    if idxtopos(sample, ridx) >= rpos
        ridx -= 1
    end

    lidx:ridx
end

#          +----------------------------------------------------------+
#          |                          Masks                           |
#          +----------------------------------------------------------+

"""
    $(FUNCTIONNAME)(h, sample, idx; wipe = true)

Construct a mask for a range of markers.

If `wipe = true`, `h` is wiped beforehand.

See also [`wipe!`](@ref).

# Methods
$(METHODLIST)

--*Internal*--
"""
function mask!(h::AbstractVector{UInt64}, sample::Sample, idx::AbstractUnitRange;
               wipe = true)
    wipe && wipe!(h)

    chunks_idx = chunkidx(Sequence, idx.start), chunkidx(Sequence, idx.stop)

    ## Evrything is in reverse in line with the way data is stored in
    ## BitVectors.

    ## First chunk
    @inbounds let nleading0 = idxinchunk(Sequence, idx.start) - 1,
        ntrailing0 = blocksize(Sequence) - (first(chunks_idx) == last(chunks_idx) ?
            idxinchunk(Sequence, idx.stop) : blocksize(Sequence))

        mask = typemax(UInt)
        mask >>>= ntrailing0
        mask &= typemax(UInt) << nleading0

        h[first(chunks_idx)] |= mask
    end

    ## Early termination in case only one chunk needs to be modified
    first(chunks_idx) == last(chunks_idx) && return h

    ## Middle chunks
    @inbounds let nchunks = last(chunks_idx) - first(chunks_idx) - 1,
        extra_chunks = nchunks % simd_chunksize,
        lane = VecRange{simd_chunksize}(0)

        for k ∈ range(first(chunks_idx) + 1, last(chunks_idx) - 1 - extra_chunks, step = simd_chunksize)
            h[lane + k] |= typemax(UInt)
        end

        for k ∈ range(last(chunks_idx) - extra_chunks, last(chunks_idx) - 1)
            h[k] |= typemax(UInt)
        end
    end

    ## Last chunk
    @inbounds let ntrailing0 = blocksize(Sequence) - idxinchunk(Sequence, idx.stop)
        mask = typemax(UInt)
        mask >>>= ntrailing0

        h[last(chunks_idx)] |= mask
    end

    h
end

function ancestral_mask!(v::AbstractVector{UInt64}, sample::Sample, ω::Ω;
                         wipe = true)
    wipe && wipe!(v)

    idx = postoidx(sample, ω)
    (iszero ∘ length)(idx) && return v

    mask!(v, sample, idx, wipe = false)
end

function ancestral_mask!(η::Sequence, sample::Sample, ω::Ω; wipe = true)
    wipe && wipe!(η)
    η.data[postoidx(sample, ω)] .= true
    η
end

function ancestral_mask!(η, sample::Sample, ωs; wipe = true)
    wipe && wipe!(η)

    for ω ∈ ωs
        ancestral_mask!(η, sample, ω, wipe = false)
    end

    η
end

ancestral_mask(sample::Sample, x) =
    ancestral_mask!(Sequence(falses(nmarkers(sample))), sample, x, wipe = false)

function ancestral_mask!(η, sample::Sample, pos::AbstractFloat; wipe = true)
    wipe && wipe!(η)
    η[postoidx(sample, pos)] = true
    η
end

wipe!(h) = fill!(h, 0)

#          +----------------------------------------------------------+
#          |              Mutation & Recombination Rates              |
#          +----------------------------------------------------------+

export mut_rate
"""
    $(FUNCTIONNAME)(sample, scaled = true)

(Scaled) mutation rate.
"""
function mut_rate end

export rec_rate
"""
    $(FUNCTIONNAME)(sample, scaled = true)

(Scaled) recombination rate.
"""
function rec_rate end

for (f, symb) ∈ Dict(:mut_rate => :(:μ), :rec_rate => :(:ρ))
    @eval function $f(sample::Sample, scaled = true)
        ret = getfield(sample, $symb)
        scaled || return ret
        4 * ret * sample.Ne * sample.sequence_length
    end
end
