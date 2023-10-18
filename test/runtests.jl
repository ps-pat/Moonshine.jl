#!/usr/bin/env julia

using Moosh

import Aqua

Aqua.test_all(Moosh,
              project_extras = false,
              ambiguities = false,
              unbound_args = false,
              piracy = false,
              deps_compat = false)

using Test:
    @test,
    @testset,
    @inferred

using JCheck:
    Quickcheck,
    @add_predicate,
    @quickcheck

using RandomNumbers.PCG: PCGStateOneseq

include("generators.jl")
include("shrinkers.jl")

@time begin
    @testset begin
        include("Sequence.jl")
        include("Arg.jl")
    end
end
