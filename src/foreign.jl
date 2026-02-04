using PythonCall
using PythonCall: pyconvert_return, pyconvert_add_rule
import PythonCall: Py

import JSON

## TreeSequence conversion.
struct TreeSequence
    obj::Py
end

ispy(::TreeSequence) = true
Py(ts::TreeSequence) = ts.obj

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

## Conversion
pyconvert_TreeSequence(::Any, ts::Py) = (pyconvert_return ∘ TreeSequence)(ts)

const msprime = Ref{Py}()
const tskit = Ref{Py}()

function __init_msprime__()
    msprime[] = pyimport("msprime")
    tskit[] = pyimport("tskit")

    pyconvert_add_rule("tskit.trees:TreeSequence", TreeSequence,
                       pyconvert_TreeSequence)
end
