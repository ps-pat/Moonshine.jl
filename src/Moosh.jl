module Moosh

using Base: Fix1

using PrecompileTools

using StatsFuns: logtwo

include("Sequence.jl")

include("tools.jl")

include("CheapStack.jl")

include("Genealogy.jl")

# include("Arg.jl")

include("Tree.jl")

include("Copulas/Copulas.jl")

include("ArgDensity.jl")

include("Distributions.jl")

include("ImportanceSampling.jl")

##################
# Precompilation #
##################

precompile(build!, (Xoshiro, Tree))
precompile(fit!, (Xoshiro, CopulaFrechet{Bernoulli{Float64}}, Vector{Bool}, Vector{Sequence}, Type{Tree}))

end # module Moosh
