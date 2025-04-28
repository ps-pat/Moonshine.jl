#!/usr/bin/env julia

using Moosh

import Aqua

let piracy_allowed = [
    Moosh._pointer] # See https://github.com/eschnett/SIMD.jl/issues/121

    @time Aqua.test_all(Moosh, piracies = (treat_as_own = piracy_allowed,))
end

using Test:
    @test,
    @testset,
    @inferred

using JCheck:
    Quickcheck,
    @add_predicate,
    @quickcheck

using RandomNumbers.PCG: PCGStateOneseq

using Graphs

include("generators.jl")
include("shrinkers.jl")

@time begin
    @testset begin
        include("Sequence.jl")
        # include("Tree.jl")
        # include("Arg.jl")
    end
end
