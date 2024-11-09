module Moosh

using PrecompileTools

using StatsFuns: logtwo, log2π

using IntervalSets
const AI = AbstractInterval

using LoopVectorization

using MultiFloats

using DataStructures: Stack, Queue

using Bumper

using AutoHashEquals

import Term.Progress

## Some constants.
const VertexType = Int
const ∞ = Inf

include("workarounds.jl")
include("Sequence.jl")
include("tools.jl")
include("foreign.jl")
include("CheapStack.jl")
include("Sample.jl")
include("Genealogy.jl")
include("Tree.jl")
include("Arg.jl")
include("genealogy_common.jl")
include("Copulas/Copulas.jl")
include("ArgDensity.jl")
include("Distributions.jl")
include("ImportanceSampling.jl")

function __init__()
    __init_msprime__()
end

##################
# Precompilation #
##################

@setup_workload begin
    rng = Xoshiro(42)
    n = 2
    m = 1
    H = Sample([Sequence(rng, m) for _ ∈ 1:n], 1e-8, 1e-8, 1000, 1000, [0.])
    @compile_workload begin
        tree = Tree(H)
        build!(rng, tree)
        arg = Arg(tree)
    end
end

end # module Moosh
