tree_qc = Quickcheck("Test Tree properties", serialize_fails = false)

## Validity.
@add_predicate(tree_qc,
               "Tree Valid",
               tree::Tree -> Moosh.validate(tree))

@add_predicate(tree_qc,
               "Mutation Edges Correct",
               tree::Tree -> begin
    me = mutation_edges(tree, Ω((first ∘ positions)(tree), Inf))

    for (marker, dedges) ∈ enumerate(me)
        for e ∈ edges(tree)
            m1 = sequence(tree, src(e))[marker]
            m2 = sequence(tree, dst(e))[marker]
            derived = m1 ⊻ m2

            in_dedges = e ∈ dedges
            (derived ⊻ in_dedges) && return false
        end
    end

    true
end)

@testset "Tree" begin
    @quickcheck tree_qc
end