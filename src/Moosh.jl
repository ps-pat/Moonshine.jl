module Moosh

using PrecompileTools

using StatsFuns: logtwo

using IntervalSets

## Some constants.
const VertexType = Int
const ∞ = Inf

export Ω
const Ω = Interval{:closed, :open, Float64}

include("Sequence.jl")
include("tools.jl")
include("CheapStack.jl")
include("Genealogy.jl")
include("Tree.jl")
include("Arg.jl")
include("genealogy_common.jl")
include("Copulas/Copulas.jl")
include("ArgDensity.jl")
include("Distributions.jl")
include("ImportanceSampling.jl")

##################
# Precompilation #
##################

precompile(build!, (Xoshiro, Tree))
precompile(fit!, (Xoshiro, CopulaFrechet{Bernoulli{Float64}}, Vector{Bool}, Vector{Sequence}, Type{Tree}))

end # module Moosh
