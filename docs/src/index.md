```@meta
CurrentModule = Moonshine
```

# Moonshine.jl Documentation
!!! warning
    Under active development, breaking changes are likely, and documentation is
    a work in progress.

    Expect nothing. Live frugally on surprise. Report any issue
    [to the project's Github repository](https://github.com/ps-pat/Moonshine.jl/issues).

Moonshine.jl is a Julia framework for
[coalescent](https://en.wikipedia.org/wiki/Coalescent_theory) modelling oriented
towards [ancestral recombination graph
(ARG)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011110)
inference. It provides the [AbstractGenealogy](@ref) type which leverages
[Graph.jl](https://github.com/JuliaGraphs/Graphs.jl) for convenient
implementation of graph-theoretical models of molecular evolution.

!!! info "QuickStart"
    First things first, you are probably looking for the
    [QuickStart page](quickstart.md).

    Read it. It's awesome.

## Peer Review
As of today, Moonshine.jl has not undergone any kind of official peer review
process, but that should happen soon enough as a paper is currently in the
works. I'll keep this section updated for those of you who care about scientific
rigour.

On that matter, note that this package places a strong emphasis on repetability
via design choices that will be documented soon. Much of it is made possible via
[Melissa O'Neil](https://www.cs.hmc.edu/~oneill/index.html)'s
[PCG family of random number generator](https://www.pcg-random.org/). Her
website and, in particular, [her blog](https://www.pcg-random.org/blog/),
definitely deserve a read if you care about random number generation,
especially when done concurrently.

## Inference
Out of the box, types representing coalescent trees and ancestral recombination
graphs are provided: see [`Tree`](@ref) and [`Arg`](@ref) respectively.
Efficient inference algorithms are implemented: see [`build!`](@ref).

## Package Ecology
### msprime
Moonshine can infer trees and ARG's directly from
[`TreeSequence` objects](https://tskit.dev/tskit/docs/latest/python-api.html#the-treesequence-class).
For convenience, it is also possible to generate a sample of haplotypes directly
from [`msprime`](https://tskit.dev/msprime/docs/stable/intro.html) with a single
line of Julia code. Both these functionalities are documented in type
[`Sample`](@ref)'s documentation.

## Conventions
Documentation implements the following conventions:
* constructors are documented with types;
* non-exported symbols are documented and marked with --*Internal*--.
