using PythonCall
using PythonCall: pyconvert_return, pyconvert_add_rule

import JSON

###########
# msprime #
###########

## TODO: finish implementation.
## TreeSequence conversion.
struct TreeSequence
    obj::Py
end

## Provenance
function provenance(ts::TreeSequence, id)
    p = ts.obj.provenance(id)
    (id = pyconvert(Int, p.id),
     timestamp = pyconvert(String, p.timestamp),
     record = JSON.parse(pyconvert(String, p.record), dicttype = Dict{Symbol, Any}))
end

function provenances(ts::TreeSequence)
    n = pyconvert(Int, ts.obj.num_provenances)
    map(id -> provenance(ts, id), 0:(n-1))
end

## Variants
function variants(ts::TreeSequence)
    ts.obj.variants()
end

## Conversion
pyconvert_TreeSequence(S, ts::Py) = (pyconvert_return ∘ TreeSequence)(ts)

const msprime = Ref{Py}()

function __init_msprime__()
    msprime[] = pyimport("msprime")
    pyconvert_add_rule("tskit.trees:TreeSequence", TreeSequence,
                       pyconvert_TreeSequence)
end