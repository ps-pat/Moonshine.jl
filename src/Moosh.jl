module Moosh

using PrecompileTools

using StatsFuns: logtwo, log2π

using IntervalSets

using LoopVectorization

using MultiFloats

## Some constants.
const VertexType = Int
const ∞ = Inf

export Ω
const Ω = Interval{:closed, :open, Float64}

include("Sequence.jl")
include("tools.jl")
include("Sample.jl")
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

@setup_workload begin
    rng = Xoshiro(42)
    n = 2
    m = 1
    @compile_workload begin
        H = [Sequence(rng, m) for _ ∈ 1:n]
        tree = Tree(H, positions = [0.], seq_length = 1., Ne = 1000, μloc = 1e-5)
        build!(rng, tree)
    end
end

end # module Moosh
