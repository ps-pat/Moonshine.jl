#!/usr/bin/env julia

using Moosh

import Aqua

Aqua.test_all(Moosh,
              project_extras = false,
              ambiguities = false,
              unbound_args = false,
              piracy = false)

using Test:
    @test,
    @testset,
    @inferred

using JCheck:
    Quickcheck,
    @add_predicate,
    @quickcheck

include("generators.jl")
include("shrinkers.jl")

@time begin
    @testset begin
        include("Sequence.jl")
    end
end
