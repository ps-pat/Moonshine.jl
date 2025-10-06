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
