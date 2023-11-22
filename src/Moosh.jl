module Moosh

using Base: Fix1

using PrecompileTools

include("Sequence.jl")

export Sequence

export fillseq

include("tools.jl")

include("CheapStack.jl")

include("Genealogy.jl")

# include("Arg.jl")

include("Tree.jl")

include("Copulas/Copulas.jl")

include("ArgDensity.jl")

## Precompilation.
@setup_workload begin
    using Random: Xoshiro

    @compile_workload begin
        rng = Xoshiro(42)

        tree = Tree(rng, 10, 10, μ_loc = 5e-5)
        build!(rng, tree)
    end
end

end # module Moosh
