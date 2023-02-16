module Moosh

## Sequence
include("Sequence.jl")

export Sequence
export fillseq

include("Arg.jl")

export Arg

export
    nleaves,
    leaves,
    isleaf,
    nmarkers,
    sequences,
    latitude,
    mrca,
    tmrca

export buildtree!

export argplot

end # module Moosh
