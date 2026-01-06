## [0.3.8] - 2026-01-06

### 🎨 Styling

- Formatting of `Project.toml`

### ⚙️ Miscellaneous Tasks

- Bump compat entries
- Bump dependencies

## [0.3.7] - 2025-11-28

### 🚀 Features

- *(Arg)* Improve derived recoalescence edge & latitude sampling

### 🐛 Bug Fixes

- *(Sample)* Revert to old default values for `Sample` constructors

### 📚 Documentation

- Remove `clat_shortcut`, replace `clat_batchsize` by `clat_gridsize`

## [0.3.6] - 2025-11-24

### 🐛 Bug Fixes

- Signatures of `mrca`/`tmrca`
- Sort `plot_tmrcas` grid

### ⚡ Performance

- *(Arg)* Implement wild recoalescence latitude sampling via time-scale transformation

## [0.3.5] - 2025-11-19

### 🚀 Features

- *(Arg)* Customizable batch size for recoalescence latitude sampling
- *(Arg)* Shortcut recoalescence latitude sampling

### 📚 Documentation

- Collapsible menus

## [0.3.4] - 2025-11-16

### 🚀 Features

- Use preferences for some constants
- *(Sample)* `BitMatrix` constructor

### 🐛 Bug Fixes

- `PythonCall` warning

### 📚 Documentation

- *(Readme)* Fix typo in Readme.md
- Bump Julia & `Documenter.jl`
- Add `CITATION.cff`
- *(Sample)* Improve documentation of `Sample`'s constructors.

### ⚡ Performance

- *(Arg)* Remove useless test in `sample_recombination_constrained!`

### ⚙️ Miscellaneous Tasks

- *(Sample)* Change `Sample` constructors default values for genetic parameters
- Bump dependencies

## [0.3.3] - 2025-10-20

### 🐛 Bug Fixes

- *(mmn)* Skip recombination edges
- *(ThreeTree)* Method ambiguity: `has_vertex(::ThreeTree, ...)`

### ⚡ Performance

- *(ThreeTree)* Remove useless root test from inneighbors
- *(mmn)* Optimize `block_predicate(::EdgesIterMMN, e)`
- *(Sequence)* Use `mod1` and (specialized) `fld1` for `idxinchunk`/`chunkidx`

## [0.3.2] - 2025-10-11

### 🐛 Bug Fixes

- `@inline` -> `@inbounds`
- Do not assume `fadjlist`/`badjlist` attributes in `AbstractGenealogy` subtypes

### 📚 Documentation

- *(ThreeTree)* Documentation for type `ThreeTree`

### ⚡ Performance

- *(ThreeTree)* Implement type `ThreeTree`
- *(Tree)* Encode `Tree` topologies with `ThreeTree`
- *(ThreeTree)* Add recombination/recoalescence methods for `ThreeTree`
- *(Arg)* Encode `Arg` topologies with `ThreeTree`

### ⚙️ Miscellaneous Tasks

- Bump Julia version to 1.12.0

## [0.3.1] - 2025-10-07

### 🐛 Bug Fixes

- *(ARG construction)* Increase size of edges & vertices stacks
- `@inline` -> `@inbounds`

## [0.3.0] - 2025-10-06

### 🚀 Features

- *(sampling)* Enlarge the support of recombination positions
- *(sampling)* Implement gene conversion
- Implement `isrecoalescence`

### 🐛 Bug Fixes

- [**breaking**] `mutationsidx!`/`mutation_edges!`/`mutation_edges`

### 🚜 Refactor

- [**breaking**] `ancestral_mask!`/`ancestral_mask`
- [**breaking**] `nmutations`
- *(graph traversal)* Create abstract types

### 📚 Documentation

- Start using git cliff

### ⚡ Performance

- *(mmn)* Exclude outedges of recombination vertices from mmn algorithm
- Early termination for `isdisjoint`
- *(mmn)* Implements bottom-up marginal graph traversal with early termination

