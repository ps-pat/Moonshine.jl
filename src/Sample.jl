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
* `H`: haplotypes
* `μ` (-1): mutation rate
* `ρ` (-1): recombination rate
* `Ne` (1): effective population size
* `positions`: position of markers
* `sequence_length` (`maximum(positions)`): sequence length
* `n`: number of sequences
* `rng`: random number generator
* `ts`: [Tree Sequence](https://tskit.dev/tskit/docs/latest/python-api.html#the-treesequence-class)
* `mat`: matrix of haplotypes (haplotypes -> columns)
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
                μ = -1, ρ = -1, Ne = 1,
                positions = 1:(length ∘ first)(H),
                sequence_length = maximum(positions))
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

function Sample(mat::BitMatrix, positions::AbstractVector{<:Real};
                μ = -1, ρ = -1, Ne = 1)
    nmarkers, n = size(mat)
    chunks = mat.chunks

    H = Vector{Sequence}(undef, n)
    for k ∈ eachindex(H)
        firstidx = (k - 1) * nmarkers + 1
        first_chunkidx = chunkidx(64, firstidx)
        first_idxinchunk = idxinchunk(64, firstidx)
        ntocopy = nmarkers + first_idxinchunk - 1

        ## TODO: lazy, optimize sometime
        data = BitVector(undef, ntocopy)
        unsafe_copyto!(data.chunks, 1, chunks, first_chunkidx, length(data.chunks))
        data <<= first_idxinchunk - 1
        resize!(data, nmarkers)
        H[k] = Sequence(data)
    end

    sequence_length = last(positions) - first(positions) + 1
    Sample(H, μ = μ, ρ = ρ, Ne = Ne,
           sequence_length = sequence_length, positions = positions)
end

###################
# Pretty printing #
###################

function show(io::IO, m::MIME"text/plain", sample::Sample)
    invoke(show, Tuple{IO, MIME"text/plain", AbstractVector}, io, m, sample)

    print(io, "\nsize = " * (string ∘ length)(sample) *
              ", length = " * string(sample.sequence_length) *
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

function ancestral_mask!(mask, sample::Sample, ω::Ω; wipe = true)
    wipe && wipe!(mask)

    idx = postoidx(sample, ω)
    (iszero ∘ length)(idx) && return mask

    ancestral_mask!(mask, nmarkers(sample), idx, wipe = false)
end

ancestral_mask(sample::Sample, ω::Ω) =
    ancestral_mask!(BitVector(undef, nmarkers(sample)), sample, ω)

function ancestral_mask!(mask, sample::Sample, ωs; wipe = true)
    wipe && wipe!(mask)

    for ω ∈ ωs
        ancestral_mask!(mask, sample, ω, wipe = false)
    end

    mask
end

ancestral_mask(sample::Sample, ωs) =
    ancestral_mask!(BitVector(undef, nmarkers(sample)), sample, ωs)

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
