export lorax
"""
    $(FUNCTIONNAME)([x; host = "127.0.0.1", port = 3000])

Visualize ancestries interactively with Lorax.

`x` can either be an [`AbstractGenealogy`](@ref) or the path to a file in a
format supported by Lorax. If argument `x` is missing, Lorax is launched in
interactive mode.

`host` and `port` default to Lorax default values, see [https://lorax.ucsc.edu/documentation](https://lorax.ucsc.edu/documentation).

# Example

```@example
arg = Arg(Xoshiro(42), 10, 1e-7, 1e-7, 10000, 1e6)
build!(Xoshiro(666), arg)
lorax(arg)
```
"""
function lorax end

function lorax(path = nothing; host = "127.0.0.1", port = 3000)
    c = ["lorax"]
    isnothing(path) || push!(c, "--file", path)
    push!(c, "--host", host)
    push!(c, "--port", string(port))
    run(Cmd(c))
end

function lorax(arg::AbstractGenealogy; host = "127.0.0.1", port = 3000)
    path = tempname(suffix = ".trees")
    ts(arg).dump(path)
    lorax(path, host = host, port = port)
end
