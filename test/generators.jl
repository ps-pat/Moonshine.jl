using Random: AbstractRNG

import JCheck: specialcases, generate

## Sequences
specialcases(::Type{Sequence{T}}) where T = [Sequence()]

generate(rng::AbstractRNG, ::Type{Sequence{T}}, n)  where T =
    [Sequence(rng, 1, 200) for _ ∈ 1:n]
