module Moosh

using Base: Fix1

using Chain: @chain

using PrecompileTools

include("Sequence.jl")

export Sequence

export fillseq

include("Arg.jl")

include("ArgDensity.jl")

include("ImportanceSampling.jl")

include("CheapStack.jl")

include("tools.jl")

## Precompilation.
@setup_workload begin
    using Random: Xoshiro

    @compile_workload begin
        rng = Xoshiro(42)

        test = Arg(rng, 10, 10, μ_loc = 5e-5, ρ_loc = 1e-8)
        buildtree!(rng, test)
    end
end

end # module Moosh
