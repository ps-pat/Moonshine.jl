```@meta
CurrentModule = Moonshine
ShareDefaultModule = true
```

# Interactive visualization with Lorax

Standard graph visualization tools such as [GraphMakie](https://github.com/MakieOrg/GraphMakie.jl) (made available by [`plot_genealogy`](@ref)) are of limited use when dealing with even moderately large ancestries. For that purpose, we recommend [Lorax](https://lorax.ucsc.edu/documentation#introduction), a web platform specifically designed for interactive visualization of ancestral recombination graphs. While [a public instance is available](https://lorax.ucsc.edu), it has also been packaged with Moonshine since version 0.4.3 and can be used entirely locally via the [`lorax`](@ref) method. Below are three ways to do so.

## Visualize ancestries on the fly
The easiest way to use Lorax is to call [`lorax`](@ref) on an object that has a [`ts`](@ref) method defined (e.g., an [`Arg`](@ref)). Conversion to `TreeSequence` is handled automatically, and the resulting object is loaded into a new local Lorax instance.

```julia
arg = Arg(rng, 10, 1e-7, 1e-7, 10000, 1e6)
build!(rng, arg)
lorax(arg)
```

Lorax binds to port 3000 on localhost by default. This can be configured via the `host` and `port` keyword arguments.

## Visualize an ARG stored in a file
[`lorax`](@ref) accepts a path to a file compatible with Lorax ([details](https://lorax.ucsc.edu/documentation#supported-inputs)).

```julia
lorax("my-arg.trees")
lorax("my-compressed-arg.trees.tsz")
lorax("my-newick-encoded-arg.csv")
```

## Interactive session
[`lorax`](@ref) can be called without arguments:

```julia
lorax()
```

In that case, ARGs can be loaded from files interactively via the web interface, similar to [https://lorax.ucsc.edu](https://lorax.ucsc.edu).

# References
- Pratik Katte, Russell Corbett-Detig, Interactive exploration of biobank-scale ancestral recombination graphs with Lorax, *Bioinformatics*, Volume 42, Issue 7, July 2026, btag458, [https://doi.org/10.1093/bioinformatics/btag458](https://doi.org/10.1093/bioinformatics/btag458)
- Lorax: Interactive visualization of Ancestral Recombination Graphs. [https://github.com/pratikkatte/lorax](https://github.com/pratikkatte/lorax)
