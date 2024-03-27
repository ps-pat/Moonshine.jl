tree_qc = Quickcheck("Test Tree properties", serialize_fails = false)

## Validity.
@add_predicate(tree_qc,
               "Tree Valid",
               tree::Tree -> isvalid(tree))
