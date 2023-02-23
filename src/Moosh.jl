module Moosh

using Base: Fix1

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
    tmrca,
    children,
    parents

export
    buildtree!,
    first_inconsistent_position

export argplot

end # module Moosh
