module Moonshine

using StatsFuns: logtwo, log2π

using IntervalSets

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

using Preferences

#          +----------------------------------------------------------+
#          |                        Constants.                        |
#          +----------------------------------------------------------+

const ∞ = Inf

"""
    const AI = AbstractInterval

Abstract types for intervals from
[`IntervalSets`](https://juliamath.github.io/IntervalSets.jl/stable/).

See also [`IntervalSets.AbstractInterval`](@extref) and [`Ω`](@ref).

--*Internal*--
"""
const AI = AbstractInterval

"""
    const VertexType = Int32
Type of the vertices. Any `genealogy` that is an instance of
[`AbstractGenealogy`](@ref) should be constructed so that `graph(genealogy)` is
an instance of `AbstractGraphs{VertexType}`
(see [`Graphs.AbstractGraph`](@extref).)

--*Internal*--
"""
const VertexType = @load_preference("VertexType", Int32)

"""
    const mmn_chunksize = 1
Size (in *bytes*) of a chunk of markers in [`next_inconsistent_idx`](@ref).

--*Internal*--
"""
const mmn_chunksize = @load_preference("mmn_chunksize", 1)

"""
    const mmn_chunktype = (eval ∘ Symbol)("UInt" * string(8mmn_chunksize))
Type of a chunk of markers in [`next_inconsistent_idx`](@ref).

--*Internal*--
"""
const mmn_chunktype = (eval ∘ Symbol)("UInt" * string(8mmn_chunksize))

"""
    const simd_vecsize = simdbytes()

Size (in *bytes*) of a SIMD registry on host machine.

--*Internal*--
"""
const simd_vecsize = simdbytes()

"""
    const simd_chunksize = simd_vecsize ÷ 8

Number of chunks of markers that fit in a SIMD registry.

--*Internal*--
"""
const simd_chunksize = simd_vecsize ÷ 8

"""
    const clat_batchsize = 100

Size of batches for recoalescence latitude sampling.

Customizable via Preferences.

--*Internal*--
"""
const clat_batchsize = @load_preference("clat_batchsize", 100)

"""
    clat_shortcut = 500

Maximum number of candidate recoalescence latitudes.

Customizable via Preferences.

--*Internal*--
"""
const clat_shortcut = @load_preference("clat_shortcut", 500)

#          +----------------------------------------------------------+
#          |                     Global variables                     |
#          +----------------------------------------------------------+

export default_colormap
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
include("ThreeTree.jl")
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
