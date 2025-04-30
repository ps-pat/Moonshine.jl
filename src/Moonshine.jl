module Moonshine

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

using DocStringExtensions

#          +----------------------------------------------------------+
#          |                        Constants.                        |
#          +----------------------------------------------------------+

const VertexType = Int
const ∞ = Inf
const mmn_chunksize = 1
const mmn_chunktype = (eval ∘ Symbol)("UInt" * string(8mmn_chunksize))
const simd_vecsize = simdbytes()
const simd_chunksize = (simd_vecsize ÷ 8)

#          +----------------------------------------------------------+
#          |                     Global variables                     |
#          +----------------------------------------------------------+

"""
    default_colormap = :Paired_3

Color map used by UnicodePlots.jl. See
[ColorSchemes.jl](https://juliagraphics.github.io/ColorSchemes.jl/stable/).
"""
default_colormap = :Paired_3

#          +----------------------------------------------------------+
#          |                        Inclusions                        |
#          +----------------------------------------------------------+

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

#          +----------------------------------------------------------+
#          |                      Precompilation                      |
#          +----------------------------------------------------------+

precompile(Arg, (Xoshiro, Float64, Float64, Float64, Float64, Float64))
precompile(build!, (Xoshiro, Arg))

# ----------------------------------------------------------------------

function __init__()
    __init_msprime__()
end

end
