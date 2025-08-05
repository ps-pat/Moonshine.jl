```@meta
CurrentModule = Moonshine
```

```@setup quickstart
using GraphMakie
using CairoMakie
```

# Quick Start
```@contents
Pages = ["quickstart.md"]
Depth = 2:2
```

The aim of this page is to present the package's main functionalities and guide
you in inferring your first ancestral recombination graph. Moonshine allows you
to do things like implement your own exotic model of ancestry and even crazy
stuff like analyzing ARGs as linear operators; tutorials on those more advanced
topics will follow soon.

As of writing these lines, Julia's dominance in the field of molecular biology
has not yet been achieved. We understand that many users might have little to
no experience with it. This guide is written with that in mind. The aim is for
an average Python user to have a satisfactory graph inference experience after
what we expect to be a short and somewhat enjoyable read. In particular, short
explanations of differences between Julia and more mainstream languages are
given when relevant: keep your eyes peeled for Julia👶 admonitions. On that
matter, the official documentation provides short and to the point descriptions
of noteworthy departures from
[Python](https://docs.julialang.org/en/v1/manual/noteworthy-differences/#Noteworthy-differences-from-Python),
[R](https://docs.julialang.org/en/v1/manual/noteworthy-differences/#Noteworthy-differences-from-R),
[C/C++](https://docs.julialang.org/en/v1/manual/noteworthy-differences/#Noteworthy-differences-from-C/C),
[MATLAB](https://docs.julialang.org/en/v1/manual/noteworthy-differences/#Noteworthy-differences-from-MATLAB) and [Common Lisp](https://docs.julialang.org/en/v1/manual/noteworthy-differences/#Noteworthy-differences-from-Common-Lisp).

[Feedback, especially from first-time users, is appreciated](https://codeberg.org/ptrk/Moonshine.jl/issues)

## Storing Haplotypes: the [`Sample`](@ref) Type
!!! note "Julia👶: Type vs Class"
    **{tl;dr} Types are classes (in the OOP sense)**

    In Julia, the concept of *type* is roughly equivalent to that of
    *class* in Python, C++ or MATLAB. Functionally, the expression "class `T`"
    can often be substituted for "type `T`". However, Julia's types
    underlie a generic function-based dispatch system and are more similar to
    R's S3 or Common Lisp's CLOS classes.

    [Official documentation on types](https://docs.julialang.org/en/v1/manual/types/)

Since we are *inferring* rather than *simulating* ARGs, we need some way to
store data about haplotypes of interest. In Moonshine, the [`Sample`](@ref)
type provides this functionality. Before getting into the nitty-gritty, a few
things to keep in mind about the data itself:
* only biallelic markers are supported;
* there is currently no support for ploidy other than 1;
* data should be phased (which is basically a consequence of the last restriction);
* the wild (ancestral) allele should be known and encoded as 0 for every marker.
These limitations, especially the one regarding ploidy, might change in the future.

A neat feature of Moonshine is its ability to transparently call
[`msprime`](https://tskit.dev/msprime/docs/stable/intro.html) to generate a
high quality sample. This should work out of the box: if you installed
Moonshine as described in the
[README](https://codeberg.org/ptrk/Moonshine.jl/src/branch/master/Readme.md),
you should have received a pre-packaged version of `msprime` at the same time.
You only need to:
1. import Moonshine;
2. import an RNG;
3. instantiate the RNG;
4. construct a sample.
This can be done interactively via [the Julia
REPL](https://docs.julialang.org/en/v1/stdlib/REPL/):
```@repl quickstart
using Moonshine
using Random: Xoshiro
rng = Xoshiro(42)
s = Sample(rng, 10, 1e-8, 1e-8, 1e4, 1e6)
```
The first line imports the `Moonshine` module and exposes its exported symbols.
The second line imports and exposes the `Xoshiro` symbol from the `Random`
module. The third line instantiates the `Xoshiro` type, creating an RNG seeded
at 42. As you might have guessed, the instance is stored in a variable called
`rng`. The last line creates the sample, calling
[`msprime.sim_ancestry`](https://tskit.dev/msprime/docs/stable/api.html#msprime.sim_ancestry)
under the hood. The meaning of each positional argument is detailed in
[`Sample`](@ref)'s documentation.

!!! note "Julia👶: Standard Library"
    `Random` is part of Julia's standard library, so you do not have to explicitly
    install it (although you still have to import it and its symbols). A
    complete list of the standard library’s modules is available in the
    [Official documentation](https://docs.julialang.org/en/v1/), section
    "Standard Library".

If you have a decent terminal with good Unicode support, you should get an
output similar to the one above. Otherwise, do yourself a favor and download a
modern terminal emulator
([list here](https://github.com/cdleon/awesome-terminals)).

If you're new to Julia, yes, we can be fast *and* look good. Sequences are
displayed as one-dimensional heatmaps for convenience. You can get a slightly
more detailed output using the [`plot_sequence`](@ref) method:
```@repl quickstart
plot_sequence(s[1])
```
The result is a little bit distorted (at least on my browser) but should look
alright on your terminal. Notice how we extracted the first (arrays
start at 1 over here) sequence? If you have a look at [`Sample`](@ref)'s
documentation, you will notice two things. First, the `Sample` type is a *subtype*
of *abstract* type `AbstractVector{Sequence}` (the
`<: AbstractVector{Sequence}` bit). This is similar to *inheritance*: any
argument of type `AbstractVector{Sequence}` (or a *supertype* thereof, for that
matter) is satisfiable by a `Sample` in a method call. Second, `Sample`
implements two *interfaces*: the array interface and the iteration interface.
This means we can basically treat them as arrays (like we did) and as iterables
(using a `for` loop or a higher order function such as
[`Base.argmax`](@extref)).

## Coalescent Trees: the [`Tree`](@ref) type
Now that we have some data, it's time to build our first coalescent tree. A
quick look at [`Tree`](@ref)'s documentation reveals that it can be constructed
from a sample:
```@repl quickstart
tree = Tree(s)
```
Conveniently, the REPL displays that `tree` is a coalescent tree with 10 leaves
(corresponding to the 10 sequences in `s`) and 563 markers. The tMRCA, which is
the latitude of the sample's "grand" MRCA, is 0. This is because there is a
distinction between *constructing* a tree, which merely means instantiating
`Tree`, and *building* it. The latter stands for sampling a *topology* (in a
graph-theoretical sense, vertices & edges) and *latitudes* (coalescence times).
Right now, the topology of `tree` is that of an edgeless graph with 10 isolated
vertices. Let's do something about it. We build a tree with the [`build!`](@ref)
method:
```@repl quickstart
rng_orig = copy(rng)
build!(rng, tree)
```

!!! note "Julia👶: Bang!"
    Notice the '!' character at the end of `build!`? This is called a "bang"
    so that the method would be pronounced "build-bang". An informal but widely
    adopted naming convention in Julia (and some other languages) is to end
    functions that mutate at least one argument with a bang. Technically, `build!`
    mutates `rng`, but we usually don't bother for such uses of random number
    generators. `tree`, however, is mutated, hence the bang! Conversely, you can
    generally assume bang-less functions not to mutate anything.

I added the first line for repeatability; it will come in handy later.
We can ignore it for now. Before explaining the building process any further,
we take a quick break to talk about visualization. First, a neat little
histogram of `tree`'s latitudes can be obtained via the
[`plot_latitudes`](@ref) method:
```@repl quickstart
plot_latitudes(tree)
```
While it's not going to win any beauty contest soon, it's quick, easy, and since
it's plain text, you can literally do a copy & paste and text it to someone
special. That's for the latitudes, but what about the topology? You can plot
that too, provided you installed two packages:
* [GraphMakie](https://github.com/MakieOrg/GraphMakie.jl);
* [GLMakie](https://docs.makie.org/stable/explanations/backends/glmakie) (or
  any other [Makie backend](https://docs.makie.org/v0.22/explanations/backends/backends)).
You can then plot `tree` via the [`plot_genealogy`](@ref) method:
```julia-repl
julia> using GraphMakie

julia> using GLMakie

julia> plot_genealogy(tree)
```
```@setup quickstart
save("assets/plot_genealogy_tree1.png", plot_genealogy(tree), size = (800, 600))
```
After some time, you should obtain the following plot:
![](assets/plot_genealogy_tree1.png)
This is all good and well, but you might wonder about the distribution of
`tree`. `tree` is a coalescent tree in the graph-theoretical sense, meaning that
it is a full binary tree. Since we are working conditional on a sample of
haplotypes however, it does not follow the "classical" distribution inherited
from Kingman's seminal paper. Moonshine's default distribution for coalescent
trees is related to three things:
1. the mutation rate;
2. the Hamming distance between sequences;
3. a bias parameter that favors coalescence between similar sequences.
Mutation events induce a distance between haplotypes, namely the number of
events necessary to turn one into another. For sequences of binary markers, this
is simply
[the Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance).
Moonshine assumes that mutation events are distributed according to a Poisson
process with user-defined mutation rate and sample coalescence events
accordingly. That being said, we might want to tweak the sampling distribution
in certain scenarios. One simple example is that of the *infinite site model* in
which each site is allowed to mutate at most once. This can be implemented via
a degenerate Poisson process, so to speak, where the distance between two
haplotypes is 0 if they are identical and infinite otherwise. In order to enable
the implementation of such models, distances can be biased via a user-defined
parameter. A bias value of ``b`` will add ``bd`` to a distance ``d``.

Enough talking, let's put that into practice. The mutation rate is defined by
the use as described in the [previous section](#Storing-Haplotypes:-the-[Sample](@ref)-Type)).
Bias can be tuned via the `bias0` *keyword argument* of `build!`. Valid values
go from `0` (no bias, default value) to `Inf` (strong bias). Let's build a tree
similar to the previous one but with infinite bias:
```@repl quickstart
rng = copy(rng_orig)
tree0 = Tree(s)
build!(rng, tree0, bias0 = Inf)
```
```@repl quickstart
plot_latitudes(tree0)
```
```@setup quickstart
save("assets/plot_genealogy_tree2.png", plot_genealogy(tree0), size = (800, 600))
```
```julia-repl
julia> plot_genealogy(tree0)
```
![](assets/plot_genealogy_tree2.png)
Notice the difference from the previous plot, even though we used the same sample
and RNG. Alright, one last tree example. A common requirement for spatial ARG
inference algorithms is building a tree consistent with the leftmost marker. To
accomplish that, we need to abandon Hamming distance, since it considers *all*
markers. Moonshine ships with the [`LeftM`](@ref) "distance", which is nothing
more than the discrete metric on the leftmost marker. First, let's have a look
at the situation on `tree`. We can color vertices according to the status of
their leftmost marker as follows:
```@setup quickstart
let fig = plot_genealogy(tree, Ω(positions(tree)[1:2]...))
    save("assets/plot_genealogy_tree3.png", fig, size = (800, 600))
end
```
```julia-repl
julia> plot_genealogy(tree, Ω(positions(tree)[1:2]...))
```
![](assets/plot_genealogy_tree3.png)

!!! note "Julia👶: Splat..."
    If you've never encountered it before in another language, the ellipsis
    ("...") is an operator that performs an operation known as "splatting". It
    allows using elements of a collection as positional arguments in a
    function call, rather than passing the collection itself as a single
    argument.

    For more details, see [`...`](@extref) in the official documentation.

In the method call above, [`Ω`](@ref) stands for a right semi-open interval:
`Ω(a, b)` is Moonshine's way of representing the interval ``[a, b)``. As you
might have guessed, `positions(tree)[1:2]` returns the positions of the first
two markers of `tree`. These are then splatted into `Ω` so that
`Ω(positions(tree)[1:2]...)` represents ``[p_1, p_2)`` where ``p_k`` is the
position of marker ``k``. Passing this interval as the second argument to
`plot_genealogy` tells it to only consider edges, vertices and markers that are
ancestral to it. In particular, vertices are colored according to the status
of the markers included in the interval (which is only the first marker in our
example): a vertex is red if all markers are derived, blue otherwise. This means
there are three *mutation edges* with respect to the first marker in `tree`:
``{17-1}``, ``{14-13}`` and ``{12-11}``. Let's build `tree` again,
this time using [`LeftM`](@ref) and an infinite bias, and see what happens:
```@repl quickstart
tree_left = Tree(s)
build!(copy(rng_orig), tree_left, Dist = LeftM(), bias0 = Inf)
```
```@setup quickstart
let fig = plot_genealogy(tree_left, Ω(positions(tree_left)[1:2]...))
    save("assets/plot_genealogy_tree_left.png", fig, size = (800, 600))
end
```
```julia-repl
julia> plot_genealogy(tree_left, Ω(positions(tree_left)[1:2]...))
```
![](assets/plot_genealogy_tree_left.png)
As expected, the leftmost marker mutates only once.

Before moving on to ARGs, I have to tell you about another handy constructor
for [`Tree`](@ref). Since multiple single-use samples are pretty common, it can
be somewhat cumbersome to explicitly construct them every time we want to build
a tree. To make our lives a little easier, there is a constructor that
transparently calls [`Sample`](@ref)'s random constructor and constructs a new
tree from the result. Long story short, we could have built `tree` more
succinctly:
```@repl quickstart
rng = Xoshiro(42)
tree_lazy = Tree(rng, 10, 1e-8, 1e-8, 1e4, 1e6)
build!(rng, tree_lazy)
```

## Ancestral Recombination Graphs: the [`Arg`](@ref) Type
Instances of [`Tree`](@ref) are built on top of instances of [`Sample`](@ref).
The design behind this is conceptual: coalescent trees contain all the
information present in the associated sample and more, namely a topology and
coalescence times. Of course, as far as we are concerned, this is merely a
convenient fiction since we have no direct knowledge of the tree: we literally
inferred it from the sample! But don't think about that too much. After all, we
are mathematicians (of sorts), so let's assume our trees are actual ancestries.
We can further improve on them by sampling recombination events. Hence, it
should come as no surprise that instances of [`Arg`](@ref), Moonshine's
tailor-made type to represent ARGs, are built on top of instances of
[`Tree`](@ref). You should be equally unsurprised when you learn that they can
be constructed in much the same way as [`Tree`](@ref) were
```@repl quickstart
arg_tree = Arg(tree)
```
and you should literally fall asleep from boredom when I tell you how to build
them:
```@repl quickstart
rng = Xoshiro(42)
build!(rng, arg_tree)
```
Just to make sure your lack of surprise is total, let me show you how to
plot its vertices' latitudes:
```@repl quickstart
plot_latitudes(arg_tree)
```
I won't `plot_genealogy` it, though, since, as you might have noticed, this is
quite a demanding process even for a tree with 19 vertices.

Moonshine's ARG building algorithm is of the spatial (as opposed to temporal)
kind: it iteratively adds recombination and recoalescence events to a tree until
every marker in the initial sample mutates at most once in its history. The
resulting graph is said to be *consistent* with the sample. As you have just
witnessed, unlike similar algorithms, the initial tree is not required to be
consistent with the leftmost marker. Moonshine really doesn't care what you
throw at it, as long as it is a valid ancestral recombination graph (or
coalescent tree, of which they are special cases). This opens the door to
things like MCMC sampling, which might be implemented in the future. For now,
let’s illustrate Moonshine’s functionalities further with something more
substantial:
```@repl quickstart
rng = Xoshiro(42)

@time begin
    arg = Arg(rng, 3000, 1e-8, 1e-8, 1e4, 1e6)
    build!(rng, arg)
end

plot_latitudes(arg)
```

!!! note "Julia👶: Macros"
    Julia has a very rich macro system, similar to Common Lisp's. Macros can be
    told appart from from functions by the '@' prefix in their names. If you
    are not familiar with Lisp, just keep in mind that their arguments do not need
    to be valid Julia expressions.

    [Official documentation on macros](https://docs.julialang.org/en/v1/manual/metaprogramming/#man-macros)

!!! note "Julia👶: begin...end"
    Denote a block of code. It was used above because the method of the
    [`Base.@time`](@extref) macro we are interested in only accepts a single
    argument.

    See [`begin`](@extref) in the official documentation.

A couple of things to unpack here. First, notice the usage of the `@time` macro.
In addition to the total execution time, it also informs us about the number
of memory allocations performed, the total memory usage (which is greater than
the peak memory usage), and the percentage of execution time dedicated to
garbage collection. The next couple of lines are similar to the ones
displayed for trees and tell us about the number of haplotypes ("leaves") and
markers associated with `arg`. We can easily plot the distribution of
breakpoints (recombination positions):
```@repl quickstart
plot_breakpoints(arg)
```

Another interesting and easily plottable feature are tMRCAs:
```@repl quickstart
plot_tmrcas(arg, npoints = 300, noprogress = true)
```
As you may have noticed, computing tMRCAs is a *slow* process. The `npoints`
argument allows trading precision for speed, as it limits computation to a grid
of approximately 200 points. By default, the tMRCA is evaluated at every
breakpoint, which is very slow. You may also have noticed that computation is
done concurrently; this helps speed things up a little. Finally, a progress
meter should be displayed. I disabled it via the `noprogress = true` parameter
because my output is static, but there is probably no reason for you to do so.

!!! note "Julia👶: Multi-Threading"
    Julia is built with multi-threading in mind. For instance, you can get the
    number of active threads via [`Base.Threads.nthreads`](@extref). To set
    the number of threads, just invoke Julia with the `-t N` switch, where `N` is
    the desired number.

---

The algorithm used to build `arg` is exact in the sense that each new event is
sampled conditional on the entire graph, or more precisely the subgraph that is
already built at the time of sampling. Of course, since we are also working
conditional on a set of haplotypes, we have to impose additional restrictions
on the distributions compared to samplers such as msprime. Moonshine is exact
in the sense that there is no SMC-type trickery involved. That being said, we
emphatically have nothing against Markovian approximations of recombination. In
fact, Moonshine is able to perform aformentioned trickery. [`build!`](@ref)
accepts a keyword argument named `winwidth` which, as its name suggests,
controls the window of positions considered when sampling an event. It is
infinite by default, leading to "exact" sampling. Setting it to 0 leads to
"approximate" SMC-type sampling with potentially significant speed gains. Let's
do a quick comparison.
```@repl quickstart
rng = Xoshiro(42)
s = Sample(rng, 1000, 1e-8, 1e-8, 1e4, 1e7);

tree_exact = Tree(s)
build!(rng, tree_exact)
tree_smc = deepcopy(tree_exact)

arg_exact = Arg(tree_exact)
arg_smc = Arg(tree_smc)

@time build!(rng, arg_exact, noprogress = true)
@time build!(rng, arg_smc, winwidth = 0, noprogress = true)
```
As you can see, there are significant gains. Memory usage is a little higher
for approximate sampling due to the increased number of recombination events,
ultimately leading to, well, more stuff to store in memory. We can actually verify
that:
```@repl quickstart
nrecombinations(arg_exact)
nrecombinations(arg_smc)
```
Of course, window widths between 0 and ``\infty`` lead to "proportional" levels
of approximation.

---

[`Tree`](@ref), [`Arg`](@ref) and any properly designed subtype of
[`AbstractGenealogy`](@ref) also implement the
[`Graphs.AbstractGraph`](@extref) interface. This ensures complete
compatibility with [Graphs.jl](https://juliagraphs.org/Graphs.jl/stable/) and,
by extension, [Julia's whole graph-theoretical
ecosystem](https://juliagraphs.org/). As an example, let's pretend we are very
interested in computing maximum flow between the first two vertices of `arg`.
No need to resort to complex code: this can be done in a couple of lines by
invoking [`GraphsFlows.maximum_flow`](@extref) from
[GraphFlows.jl](https://juliagraphs.org/GraphsFlows.jl/dev/):
```@repl quickstart
using Graphs, GraphsFlows, SparseArrays
capmat = spzeros(Float64, nv(arg), nv(arg)); # 0-initialized sparse matrix
for e ∈ edges(arg) # fill capacity matrix
    ibl = (inv ∘ branchlength)(arg, e) # capacity = 1/length of branch
    capmat[src(e), dst(e)] = capmat[dst(e), src(e)] = ibl
end
flow, flowmat = maximum_flow(arg, 1, 2, capmat)
droptol!(flowmat, (eps ∘ eltype)(flowmat)) # optional cleanup
```

!!! note "Julia👶: function ∘ composition"
    **{tl;dr} `(f ∘ g)(x...)` = `f(g(x...))`**

    `∘` is obtained in the REPL via LaTeX-type syntax and tab completion: just
    type `\circ<TAB>`. Learn all about it in
    [the official documentation](https://docs.julialang.org/en/v1/manual/functions/#Function-composition-and-piping)
In addition, Moonshine implements a variety of specialized, coalescent
theory-related methods. These are described in the [AbstractGenealogy](@ref)
section. Methods specific to coalescent trees and ancestral recombination
graphs are also available and documented in [Tree & ARG](@ref). The API is not
fully documented yet, but will be improved shortly.
That's about it for now. Go on and have fun coalescing!
