#!/usr/bin/env julia

using Moosh

import Aqua

Aqua.test_all(Moosh,
              project_extras = false,
              ambiguities = false,
              unbound_args = false,
              piracies = false,
              deps_compat = false,
              persistent_tasks = false)

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
        include("Tree.jl")
        #include("Arg.jl")
    end
end
