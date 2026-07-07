```@meta
CurrentModule = Moonshine
```

# AbstractGenealogy
```@contents
Pages = ["AbstractGenealogy.md"]
Depth = 2:4
```

## Types
```@docs
AbstractGenealogy
EdgesInterval
```

### Edge Iteration
```@docs
AbstractEGIter
AbstractEGIterTD
AbstractEGIterBU
EIterTD
EIterBU
```

## Functions
### Interface
`AbstractGenealogy` interface. Subtypes must implement every method unless
otherwise specified.

#### General
```@docs
dens
describe
graph
block_predicate
```

#### Vertices/Edges
```@docs
latitudes
latitude
leaves
ivertices
maxdads
maxchildren
```

#### Sequences
```@docs
positions
sequences
```

#### Ancestrality
```@docs
ancestral_intervals
ancestral_intervals!
mrca
```

#### Coalescences
```@docs
iscoalescence
```

#### Recombinations
```@docs
recombinations
isrecombination
nrecombinations
```

#### Plotting
```@docs
plot_layout
```

### General
#### Vertices/Edges
```@docs
nleaves
nivertices
isleaf
isroot
isivertex
branchlength
edges_interval
edgesmap
nlive!
ismutation_edge
```

#### Sequences
```@docs
sequence
idxtopos
postoidx
nmarkers
nmutations
```

#### Ancestrality
```@docs
dad
dads
child
children
ancestors
descendants
siblings!
siblings
sibling
tmrca
```

#### Plotting
```@docs
plot_genealogy
plot_latitudes
```
