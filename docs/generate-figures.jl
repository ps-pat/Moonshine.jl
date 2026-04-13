using Moonshine

using Random: Xoshiro

using GraphMakie, CairoMakie

const assets = abspath("docs/src/assets")

rng = Xoshiro(42)
s = Sample(rng, 10, 1e-8, 1e-8, 1e4, 1e6)

tree = Tree(s)

rng_orig = copy(rng)
build!(rng, tree)
save(assets * "/plot_genealogy_tree1.png", plot_genealogy(tree), size = (800, 600))

rng = copy(rng_orig)
tree0 = Tree(s)
build!(rng, tree0, bias0 = Inf)
save(assets * "/plot_genealogy_tree2.png", plot_genealogy(tree0), size = (800, 600))

save(assets * "/plot_genealogy_tree3.png", plot_genealogy(tree, Ω(positions(tree)[1:2]...)), size = (800, 600))

tree_left = Tree(s)
build!(copy(rng_orig), tree_left, Dist = LeftM(), bias0 = Inf)
save(assets * "/plot_genealogy_tree_left.png", plot_genealogy(tree_left, Ω(positions(tree_left)[1:2]...)))

rng = Xoshiro(42)
tree_lazy = Tree(rng, 10, 1e-8, 1e-8, 1e4, 1e6)
build!(rng, tree_lazy)
