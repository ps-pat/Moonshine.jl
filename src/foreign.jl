import JSON

_validate_python_class_macro(objname, classname, obj, class) =
    Bool(pybuiltins.isinstance(obj, class)) ||
    (throw ∘ ArgumentError)(
        "$objname must be an instance of $classname"
    )

macro _validate_python_class(obj, class)
    objname = string(obj)
    classname = string(class)
    eobj = esc(obj)
    eclass = esc(class)

    quote
        _validate_python_class_macro(
            $objname,
            $classname,
            $eobj,
            $eclass
        )
    end
end

const msprime = Ref{Py}()
const tskit = Ref{Py}()

function __init_msprime__()
    msprime[] = pyimport("msprime")
    tskit[] = pyimport("tskit")

    convert_provenance(::Any, provenance) = (
        id = pyconvert(Int, provenance.id),
        timestamp = pyconvert(String, provenance.timestamp),
        record = JSON.parse(
            pyconvert(String, provenance.record),
            dicttype = Dict{Symbol, Any}
        )
    )

    pyconvert_add_rule(
        "tskit.trees:Provenance",
        NamedTuple,
        convert_provenance
    )
end
