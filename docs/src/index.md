```@meta
CurrentModule = Moonshine
```
!!! warning "Documentation in progres"
    This documentation is still a work in progress. Among other things, API reference is not complete yet.

`Moonshine` is a Julia framework for [coalescent](https://en.wikipedia.org/wiki/Coalescent_theory) modelling oriented towards [ancestral recombination graph (ARG)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011110) inference. It provides the [AbstractGenealogy](@ref) type which leverages [Graph.jl](https://github.com/JuliaGraphs/Graphs.jl) for convenient implementation of graph-theoretical phylogenetic models.

!!! info "QuickStart"
    If it is your first visit, chances are you are looking for the [Quick Start](@ref) page.

## Reference
`Moonshine` is the subject of the following open access publication:

[Fournier P and Larribe F (2026). Moonshine.jl: a Julia package for genome-scale model-based ancestral recombination graph inference. Front. Genet. 17:1753780. doi: 10.3389/fgene.2026.1753780](https://doi.org/10.3389/fgene.2026.1753780)

If you use `Moonshine` in your work, please cite this article. References in standard formats are provided with the [HTML version of the paper](https://doi.org/10.3389/fgene.2026.1753780) or on the [GitHub page of the package](https://github.com/ps-pat/Moonshine.jl) ([`CITATION`](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files) file).

## Inference
Out of the box, types representing coalescent trees and ancestral recombination
graphs are provided: see [`Tree`](@ref) and [`Arg`](@ref) respectively.
Efficient inference algorithms are implemented: see [`build!`](@ref).

## Package Ecology
### msprime
* `Moonshine` can infer trees and ARG's directly from [`TreeSequence` objects](https://tskit.dev/tskit/docs/latest/python-api.html#the-treesequence-class). For convenience, it is also possible to generate a sample of haplotypes directly from [`msprime`](https://tskit.dev/msprime/docs/stable/intro.html) with a single line of Julia code. Both these functionalities are documented in type [`Sample`](@ref)'s documentation.
* [`Arg`](@ref)s can be converted to `TreeSequence`s: see [Calling Moonshine From Python](@ref).

### OpenMendel
Genomic data can be imported from VCF formatted files using the [VCFTools](https://github.com/OpenMendel/VCFTools.jl) package: see [Working with VCF files](@ref).

## Conventions
Documentation implements the following conventions:
* constructors are documented with types;
* non-exported symbols are documented and marked with --*Internal*--.
