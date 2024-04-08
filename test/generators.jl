using Random: AbstractRNG

import JCheck: specialcases, generate

############
# Sequence #
############

specialcases(::Type{Sequence}) = [Sequence()]

generate(rng::AbstractRNG, ::Type{Sequence}, n) =
    [Sequence(rng, 1, 200) for _ ∈ 1:n]

# ## Fixed length sequences.
# struct SequenceN{N}
#     sequence::Sequence
# end

# specialcases(::Type{SequenceN}) = []

# generate(rng::AbstractRNG, ::Type{Sequence{N}}, n) where N =
#     [Sequence(rng, N) for _ ∈ 1:n]

########
# Tree #
########

specialcases(::Type{Tree}) = [Tree()]

function generate(rng::AbstractRNG, ::Type{Tree}, n)
    μ_loc = 1e-7
    seq_length = 10
    Ne = 1000

    ret = [Tree(rng, 1, 1, 100, 100, seq_length = seq_length, Ne = Ne, μ_loc = μ_loc)
           for _ ∈ 1:n]

    for k ∈ 1:n
        build!(rng, ret[k])
    end

    ret
end

# #######
# # Arg #
# #######

# ## ARG on 2 leaves, 2 markers, no recoalescence.
# struct Arg2{T}
#     arg::Arg{T}
# end

# specialcases(::Type{Arg2{T}}) where T = []

# function generate(rng::AbstractRNG, ::Type{Arg2{T}}, n) where T
#     μ_loc = 5e-7
#     seq_length = 1
#     Ne = 1000

#     ret = [Arg2(Arg([Sequence(rng, 1) for _ ∈ 1:2],
#                seq_length = seq_length,
#                effective_popsize = Ne,
#                μ_loc = μ_loc,
#                positions = [0]))
#            for _ ∈ 1:n]

#     for k ∈ eachindex(ret)
#         buildtree!(rng, ret[k].arg)
#     end

#     ret
# end
