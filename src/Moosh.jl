module Moosh

using Base: Fix1

include("Sequence.jl")

export Sequence

export fillseq

include("Arg.jl")

export Arg

include("ArgDensity.jl")

end # module Moosh
