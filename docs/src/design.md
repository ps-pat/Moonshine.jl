```@meta
CurrentModule = Moonshine
```

# Naming Conventions
## Arguments
### Genealogies
* `genealogy`: [`AbstractGenealogy`](@ref)
* `tree`: coalescent tree (usually of type [`Tree`](@ref))
* `arg`: ancestral recombination graph (usually of type [`Arg`](@ref))
  - Keep in mind that by definition, a coalescent tree is an ancestral recombination graph.

### Graphs
* `v`: vertex (usually of type [`VertexType`](@ref))
* `vs`: vertices
* `e`: edge (usually of type `Edge{VertexType}`, see [`Graphs.Edge`](@extref))
* `es`: edges

### Sequences
* `h`: haplotype (usually of type [`Sequence`](@ref), sometimes a `BitVector`
  or `Vector{UInt}`
* `idx`: index or indices of markers, i.e. element or subset or the range
  `1:n` where `n` is the number of markers
* `pos`: position or positions of markers, i.e. element or subset of
  ``\mathbb R``

### Intervals
* `ω`: right semi-open interval (usually of type [`Ω`](@ref))
* `ωs`: `ω` or collection of `ω`s (usually of type [`AIs`](@ref))

## Keyword Arguments
* `buffer`: `SlabBuffer` or `AllocBuffer` from
  [Bumper.jl](https://github.com/MasonProtter/Bumper.jl)

# Constructors
For a type `T`:
* `T()` is the empty constructor
* `T(rng, ...)` is a random constructor
