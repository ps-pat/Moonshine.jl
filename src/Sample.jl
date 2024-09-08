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
                μ = 0, ρ = 0, Ne = 0, sequence_length = 1, positions = [])
    positions = _validate_positions(positions, (length ∘ first)(H))
    Sample(H, μ, ρ, Ne, sequence_length, positions)
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

    Sample(H, μ, ρ, Ne, sequence_length, positions)
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

idxtopos(sample::Sample, idx) =
    iszero(idx) ? zero(Float64) : getindex(positions(sample), idx)

function postoidx(sample::Sample, pos)
    ret = 1
    while idxtopos(sample, ret) < pos
        ret += 1
    end

    ret
end

function postoidx(sample::Sample, ω::Ω)
    lpos, rpos = endpoints(ω)

    lidx = postoidx(sample, lpos)

    ridx = nmarkers(sample)
    while ridx > 0 && idxtopos(sample, ridx) >= rpos
        ridx -= 1
    end

    lidx:ridx
end

function ancestral_mask!(η, sample::Sample, ω::Ω; wipe = true)
    wipe && _wipe!(η)
    η.data[postoidx(sample, ω)] .= true
    η
end

function ancestral_mask!(η, sample::Sample, ωs::AbstractSet{Ω}; wipe = true)
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

_wipe!(η) = η.data.chunks .⊻= η.data.chunks

for (f, symb) ∈ Dict(:mut_rate => :(:μ), :rec_rate => :(:ρ))
    @eval function $f(sample::Sample, scaled = true)
        ret = getfield(sample, $symb)
        scaled || return ret
        ret * sample.Ne * sample.sequence_length
    end
end
