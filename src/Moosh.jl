module Moosh

using Base: Fix1

include("Sequence.jl")

export Sequence

export fillseq

include("Arg.jl")

include("ArgDensity.jl")

include("ImportanceSampling.jl")

end # module Moosh
