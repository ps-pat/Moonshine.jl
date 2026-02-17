# Python & tskit
```@contents
Pages = ["python.md"]
Depth = 2:2
```

## Converting to tskit
An [Arg](@ref) can be converted to a `TreeSequence` using the [`ts`](@ref) method. From there, you have access to the complete [`TreeSequence` API](https://tskit.dev/tskit/docs/stable/python-api.html#treesequence-api).

Conversion of Python objects towards Julia is mainly done by the `pyconvert` methods ([details here](https://juliapy.github.io/PythonCall.jl/stable/conversion-to-julia/)).

## Installing Moonshine From Python
Moonshine can be installed from Python using the [juliapkg](https://github.com/JuliaPy/pyjuliapkg) package, which is available from pip. After installation, simply run:
```python
import juliapkg
juliapkg.add("Moonshine")
```

## Calling Moonshine From Python
To actually use Moonshine, you will need the [juliacall](https://juliapy.github.io/PythonCall.jl/stable/juliacall/) package. Once installed, you can load Moonshine using the `seval` method:
```python
import juliacall
juliacall.Main.seval("using Moonshine")
```
You should now have access to the entirety of Moonshine. For instance, you can construct a random arg tree:
```python
M = juliacall.Main.Moonshine

juliacall.Main.seval("using Random")
rng = juliacall.Main.Xoshiro(42)

arg = M.Arg(rng, 1000, 1e-8, 1e-8, 1e4, 1e6)
M.build_b(rng, arg)
```
Since Python does not allow the character `!` in symbol names, it is replaced with `_b` as shown in the last line. `arg` is a Julia object, but we can easily convert it to a Python `TreeSequence` with the aforementioned `ts` method:
```python
trees = M.ts(arg)
```
