module Moosh

using PrecompileTools

using StatsFuns: logtwo, log2π

using IntervalSets
const AI = AbstractInterval

using LoopVectorization

using DoubleFloats: Double64

using DataStructures: Stack, Queue

using Bumper

using UnsafeArrays: UnsafeArray

using Base: unsafe_convert

using ProgressMeter

using FunctionWrappers: FunctionWrapper

using CpuId: simdbytes

using SIMD: VecRange

using Distributions

#          +----------------------------------------------------------+
#          |                        Constants.                        |
#          +----------------------------------------------------------+

const VertexType = Int
const ∞ = Inf
const mmn_chunksize = 8
const mmn_chunktype = (eval ∘ Symbol)("UInt" * string(mmn_chunksize))
const simd_vecsize = simdbytes()

#          +----------------------------------------------------------+
#          |                     Global variables                     |
#          +----------------------------------------------------------+

"""
    default_colormap = :Paired_3

Color map used by UnicodePlots.jl. See
[ColorSchemes.jl](https://juliagraphics.github.io/ColorSchemes.jl/stable/).
"""
default_colormap = :Paired_3

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

end # module Moosh
