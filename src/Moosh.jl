module Moosh

using PrecompileTools

using StatsFuns: logtwo, log2π

using IntervalSets
const AI = AbstractInterval

using LoopVectorization

using MultiFloats

using DataStructures: Stack, Queue

using Bumper

using UnsafeArrays: UnsafeArray

using Base: unsafe_convert

using ProgressMeter

using FunctionWrappers: FunctionWrapper

## Some constants.
const VertexType = Int
const ∞ = Inf
const mmn_chunksize = 8
const mmn_chunktype = (eval ∘ Symbol)("UInt" * string(mmn_chunksize))
default_colormap = :Paired_3

include("workarounds.jl")
include("Sequence.jl")
include("tools.jl")
include("foreign.jl")
include("CheapStack.jl")
include("AncestralIntervals.jl")
const AIsType = AIs{Vector{Ω}}
include("Sample.jl")
include("Genealogy.jl")
include("Tree.jl")
include("Arg/Arg.jl")
include("genealogy_common.jl")
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
