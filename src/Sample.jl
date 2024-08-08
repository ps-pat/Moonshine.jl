import Base: iterate,
             eltype,
             length,
             size,
             getindex,
             IndexStyle


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

function Sample(H; μ = 0, ρ = 0, Ne = 0, sequence_length = 1, positions = [])
    positions = _validate_positions(positions, (length ∘ first)(H))
    Sample(H, μ, ρ, Ne, sequence_length, positions)
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

@generated IndexStyle(::Type{Sequence}) = Base.IndexLinear()

#################
# Short Methods #
#################

positions(sample::Sample) = sample.positions

nmarkers(sample::Sample) = (length ∘ first)(sample.H)

idxtopos(sample::Sample, idx) = getindex(positions(sample), idx)

function postoidx(sample::Sample, pos)
    @inbounds for (k, p) ∈ enumerate(positions(sample))
        p > pos && return k - 1
    end

    nmarkers(sample)
end

"""
    ancestral_mask!(η, x, ω; wipe = true)
    ancestral_mask(x, ω)

Mask non ancestral positions to 0. If `wipe = true`, all markers in `η` wil be
initialized at 0.
"""
function ancestral_mask! end,
function ancestral_mask end

function ancestral_mask!(η, sample::Sample, ω::Ω; wipe = true)
    lpos, rpos = endpoints(ω)

    lidx = 1
    while lpos > positions(sample)[lidx]
        lidx += 1
    end

    ridx = postoidx(sample, rpos)

    wipe && _wipe!(η)
    η.data[range(lidx, ridx)] .= true

    η
end

function ancestral_mask!(η, sample::Sample, ωs::Set{Ω}; wipe = true)
    wipe && _wipe!(η)

    for ω ∈ ωs
        ancestral_mask!(η, sample, ω, wipe = false)
    end

    η
end

ancestral_mask(sample::Sample, ω) =
    ancestral_mask!(Sequence(falses(nmarkers(sample))), sample, ω,
                    wipe = false)

ancestral_mask!(η, ωs, sample::Sample, x::Union{VertexType, Edge}; wipe = true) =
    ancestral_mask!(η, sample, ancestral_intervals!(ωs, sample, x), wipe = wipe)

ancestral_mask(sample::Sample, x::Union{VertexType, Edge}) =
    ancestral_mask!(Sequence(undef, nmarkers(sample)), Set{Ω}(), sample, x)

function ancestral_mask!(η, sample::Sample, x::AbstractFloat; wipe = true)
    wipe && _wipe!(η)
    η[postoidx(sample, x)] = true
    η
end

_wipe!(η) = η.data.chunks .⊻= η.data.chunks
