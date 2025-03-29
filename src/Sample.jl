import Base: iterate,
             eltype,
             length,
             size,
             getindex,
             IndexStyle

using Graphs: Edge

export Ω
const Ω = Interval{:closed, :open, Float64}

export Sample
"""
    struct Sample

Contain a sample of haplotypes and informations about them.
"""
struct Sample <: AbstractVector{Sequence}
    H::Vector{Sequence}
    μ::Float64
    ρ::Float64
    Ne::Float64
    sequence_length::Float64
    positions::Vector{Float64}
    coefs::NTuple{2, Float64}
end

"""
    _validate_positions(positions, nmarkers)

Validate a marker of positions. Tries to fix it if invalid. Return a vector
of uniformly spaced positions if it cannot be fixed.
"""
function _validate_positions(positions, nmarkers)
    if length(positions) != nmarkers || !issorted(positions)
        @info((isempty(positions) ? "A" : "Invalid positions: a") *
              "ssuming equally spaced markers")
        positions = isone(nmarkers) ?
                    [0.0] : (collect ∘ range)(0, 1, length = nmarkers)
    end
    if minimum(positions) < 0
        @info("First position is < 0: shifting positions right")
        positions .-= minimum(positions)
    end
    if maximum(positions) > 1
        @info("Last position is > 0: scaling positions")
        positions ./= maximum(positions)
    end

    positions
end

function Sample(H::AbstractVector{Sequence};
                μ = 1e-8, ρ = 0, Ne = 10_000,
                sequence_length = length(H) + 1,
                positions = 1:length(H))
    # positions = _validate_positions(positions, (length ∘ first)(H))

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

function ancestral_mask!(v::AbstractVector{UInt64}, sample::Sample, ω::Ω;
                         wipe = true)
    wipe && _wipe!(v)

    idx = postoidx(sample, ω)
    iszero(length(idx)) && return v

    chunks_idx = chunkidx(Sequence, idx.start), chunkidx(Sequence, idx.stop)

    ## Evrything is in reverse in line with the way data is stored in
    ## BitVectors.

    ## First chunk
    let nleading0 = idxinchunk(Sequence, idx.start) - 1,
        ntrailing0 = blocksize(Sequence) - (first(chunks_idx) == last(chunks_idx) ?
            idxinchunk(Sequence, idx.stop) : blocksize(Sequence))

        mask = typemax(UInt)
        mask >>>= ntrailing0
        mask &= typemax(UInt) << nleading0

        v[first(chunks_idx)] |= mask
    end

    ## Early termination in case only one chunk needs to be modified
    first(chunks_idx) == last(chunks_idx) && return v

    ## Middle chunks
    @turbo for k ∈ range(first(chunks_idx) + 1, last(chunks_idx) - 1)
        v[k] |= typemax(UInt)
    end

    ## Last chunk
    let ntrailing0 = blocksize(Sequence) - idxinchunk(Sequence, idx.stop)
        mask = typemax(UInt)
        mask >>>= ntrailing0

        v[last(chunks_idx)] |= mask
    end

    v
end

function ancestral_mask!(η::Sequence, sample::Sample, ω::Ω; wipe = true)
    wipe && _wipe!(η)
    η.data[postoidx(sample, ω)] .= true
    η
end

function ancestral_mask!(η, sample::Sample, ωs; wipe = true)
    wipe && _wipe!(η)

    for ω ∈ ωs
        ancestral_mask!(η, sample, ω, wipe = false)
    end

    η
end

ancestral_mask(sample::Sample, x) =
    ancestral_mask!(Sequence(falses(nmarkers(sample))), sample, x, wipe = false)

function ancestral_mask!(η, sample::Sample, pos::AbstractFloat; wipe = true)
    wipe && _wipe!(η)
    η[postoidx(sample, pos)] = true
    η
end

_wipe!(η::Sequence) = η.data.chunks .⊻= η.data.chunks

_wipe!(h) = fill!(h, 0)

for (f, symb) ∈ Dict(:mut_rate => :(:μ), :rec_rate => :(:ρ))
    @eval function $f(sample::Sample, scaled = true)
        ret = getfield(sample, $symb)
        scaled || return ret
        4 * ret * sample.Ne * sample.sequence_length
    end
end
