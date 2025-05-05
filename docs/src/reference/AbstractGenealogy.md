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

## Functions
### Interface
`AbstractGenealogy` interface. Subtypes must implement every method unless
otherwise specified.

#### General
```@docs
dens
describe
graph
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
nlive
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
nmutations!
```

#### Ancestrality
```@docs
dad
dads
child
children
ancestors
descendants
sibling
siblings
ancestral_mask
ancestral_mask!
tmrca
```

#### Plotting
```@docs
graphplot(::AbstractGenealogy, ::Any)
plot_latitudes
```
