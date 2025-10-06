using Graphs

import Graphs: add_vertices!, add_edge!, rem_edge!

using Random

using StatsBase: samplepair, ProbabilityWeights, fit, Histogram

using StaticArrays: @SVector

using UnicodePlots: heatmap, label!, stairs

##################
# Arg Definition #
##################

export Arg
"""
    $(TYPEDEF)

Ancestral recombination graph.

See also [`Tree`](@ref).

# Fields
$(TYPEDFIELDS)

# Constructors
!!! info
    Random constructor calls [`Sample`](@ref)'s random constructor.

!!! warning
    These **do not** actually build the arg. For that, see
    [`build!(rng, arg)`](@ref).

$(METHODLIST)

# Arguments
* `tree`: coalescent [`Tree`](@ref)
* other arguments are identical to [`Sample`](@ref)
"""
struct Arg <: AbstractGenealogy
    "Graph's topology"
    graph::SimpleDiGraph{VertexType}
    "Vertices' latitudes"
    latitudes::Vector{Float64}
    "∩-mask for ancestral intervals"
    recombination_mask::Vector{AIsType}
    "Arg's grand MRCA"
    mrca::Base.RefValue{VertexType}
    "Vertices' haplotypes"
    sequences::Vector{Sequence}
    "Edges' ancestral intervals"
    ancestral_intervals::Dict{Edge{VertexType}, AIsType}
    "Associated [`Sample`](@ref)"
    sample::Sample
    "Log-value of the associated pdf"
    logdensity::Base.RefValue{Double64}
end

function Arg(tree::Tree)
    ancestral_intervals = Dict{Edge{VertexType}, AIsType}()
    for e ∈ edges(tree)
        ancestral_intervals[e] = AIsType([Ω(0, ∞)])
    end

    Arg(graph(tree),
        latitudes(tree),
        Vector{AIsType}(undef, 0),
        Ref(mrca(tree)),
        sequences(tree),
        ancestral_intervals,
        sam(tree),
        tree.logdensity)
end

function Arg(rng::AbstractRNG, n, μ, ρ, Ne, sequence_length)
    tree = Tree(rng, n, μ, ρ, Ne, sequence_length)
    build!(rng, tree)
    Arg(tree)
end

###############################
# AbstractGenealogy Interface #
###############################

nleaves(arg::Arg) = length(sam(arg).H)

describe(::Arg, long = true) = long ? "Ancestral Recombination Graph" : "ARG"

isrecombination(::Arg, v, n) = iseven(v) && v >= 2n

isrecombination(arg::Arg, v) = isrecombination(arg, v, nleaves(arg))

function iscoalescence(arg::Arg, v)
    n = nleaves(arg)

    v <= n && return false
    v < 2n && return true
    isodd(v) && return true

    false
end

function recombinations(arg::Arg)
    isempty(arg.recombination_mask) && return StepRange{Int, Int}(0, 1, 0)

    start = 2nleaves(arg)
    step = 2
    stop = nv(arg)

    StepRange{Int, Int}(start, step, stop)
end

nrecombinations(arg::Arg) = ne(arg) - nv(arg) + 1

mrca(arg::Arg) = arg.mrca[]

mrca(arg, ωs) = mrca(arg, leaves(arg), ωs)

@generated maxdads(::Type{Arg}) = 2

@generated maxchildren(::Type{Arg}) = 2

let funs = (:dads, :children),
    typesandtests = ((:Real, in), (:AI, !isdisjoint), (:(AIs{<:AbstractVector{<:AI}}), !isdisjoint))

    for fun ∈ funs
        edge = fun == :dads ? :(Edge(u => v)) : :(Edge(v => u))
        isrecmodifier = fun == :dads ? :! : :identity

        for (Argtype, testfun) ∈ typesandtests
            @eval @inbounds function $fun(arg, v, ωs::$Argtype)
                neig = $fun(arg, v)
                isempty(neig) && return view(neig, 0x01:0x00)

                u = first(neig)
                ai1 = ancestral_intervals(arg, $edge)

                if $isrecmodifier(isrecombination)(arg, v)
                    idx = UInt8($testfun(ωs, ai1))
                    return view(neig, 0x01:idx)
                end

                u = last(neig)
                ai2 = ancestral_intervals(arg, $edge)
                idx = 0x06
                idx ⊻= 0x03 * $testfun(ωs, ai1)
                idx ⊻= 0x0c * $testfun(ωs, ai2)

                view(neig, range(idx & 0x03, idx >> 0x02))
            end
        end
    end
end

###########################
# AbstractGraph Interface #
###########################

function add_vertices!(arg::Arg, H, lats)
    append!(latitudes(arg), lats)
    append!(sequences(arg), H)
    add_vertices!(graph(arg), length(H))
end

function add_edge!(arg::Arg, e, ints::AIsType)
    arg.ancestral_intervals[e] = ints
    add_edge!(graph(arg), e)
end

function rem_edge!(arg::Arg, e)
    delete!(arg.ancestral_intervals, e)
    rem_edge!(graph(arg), e)
end

function plot_layout(arg::Arg)
    initxs = rand(1:nleaves(arg), nv(arg))
    ys = vcat(zeros(nleaves(arg)), latitudes(arg))
    Spring(
        initialpos = (collect ∘ zip)(initxs, ys),
        # initialpos = vcat(
            # [(v, 0) for v ∈ leaves(arg)],
            # [(1, 10) for _ ∈ ivertices(arg)]
        # ),
        pin = [(false, true) for _ ∈ vertices(arg)])
end

########################
# Ancestrality Methods #
########################

function ancestral_intervals!(ωs, arg::Arg, e::Edge; wipe = true)
    wipe && empty!(ωs)

    union!(ωs, ancestral_intervals(arg, e))

    ωs
end

ancestral_intervals(arg::Arg, e::Edge) = arg.ancestral_intervals[e]

function ancestral_intervals!(ωs, arg::Arg, v::VertexType; wipe = true)
    wipe && empty!(ωs)
    isleaf(arg, v) && return push!(ωs, Ω(0, ∞))

    for child ∈ children(arg, v)
        ancestral_intervals!(ωs, arg, Edge(v => child), wipe = false)
    end

    ωs
end

ancestral_intervals(arg::Arg, v::VertexType) = ancestral_intervals!(AIsType(), arg, v)

_recidx(arg, v) = (v - 2(nleaves(arg) - 1)) ÷ 2

function _coalidx(arg, v)
    n = nleaves(arg)
    if v > 2n
        v -= _recidx(arg, v)
    end
    v - n
end

export recombination_mask
"""
    $(SIGNATURES)

    Intervals of the partition of [0, ∞) induced by a recombination event for
    which an edge is ancestral.
"""
function recombination_mask(arg, e::Edge)
    s, d = src(e), dst(e)

    inc = s > otherdad(arg, s, d)
    arg.recombination_mask[2_recidx(arg, d) - 1 + inc]
end

export breakpoints
"""
    $(SIGNATURES)

Recombination events' positions
"""
breakpoints(arg::Arg) =
    Iterators.map(rightendpoint, @view arg.recombination_mask[1:2:end])

export plot_breakpoints
"""
    $(SIGNATURES)

Heatmap of recombination events' positions.

See also [`breakpoints`](@ref)

# Keywords
* `nbins` (`clamp(nrecombinations(arg) ÷ 100, 1, 69)`): number of bins. Its
  default maximum value of 69 produces a nice 80 columns wide plot.
* `height` (`7`): height of the plot.
* `kwargs...`: additional keywords arguments are passed directly to
  [`UnicodePlots.histogram`](https://github.com/JuliaPlots/UnicodePlots.jl).
"""
function plot_breakpoints(arg;
                          nbins = clamp(nrecombinations(arg) ÷ 100, 1, 69),
                          height = 7,
                          kwargs...)
    bins = range(0, sam(arg).sequence_length, length = nbins + 1)
    bps = (collect ∘ breakpoints)(arg)
    h = fit(Histogram, bps, bins)
    counts = repeat(reshape(h.weights, (1, nbins)), 10)

    plt = heatmap(counts,
                  width = nbins,
                  border = :none,
                  margin = 0,
                  title = "Recombinations' Positions",
                  height = height,
                  colorbar = true,
                  zlabel = "#ρ",
                  colormap = default_colormap;
                  kwargs...)

    for k ∈ 1:height
        label!(plt, :l, k, "")
    end

    mid(x) = getindex(x, length(x) ÷ 2)
    label!(plt, :bl, (string ∘ Int ∘ round ∘ first ∘ positions)(arg), color = :white)
    label!(plt, :br, (string ∘ Int ∘ round ∘ last ∘ positions)(arg), color = :white)
    label!(plt, :b, (string ∘ Int ∘ round ∘ mid ∘ positions)(arg), color = :white)

    plt
end

export plot_tmrcas
"""
    $(SIGNATURES)

Staircase plot of time to the most recent common ancestor. Additional keywords
arguments are passed directly to [`UnicodePlots.histogram`](@ref).

This function actually computes the tmrca of each recombinatoin breakpoint,
which is **very** demanding, much more than simulating the graph. It is
multithreaded in order to reduce computation time.

See also [`tmrca`](@ref)

# Keywords
* `npoints` (`nrecombinations(arg) + 1`): approximate number of points at which
  to compute the TMRCA. By default, every recombination breakpoint is
  considered. Use a smaller number of points to reduce computing time.
* `width` (`76`): width of the plot. Its default maximum value of 76 produces
  a 80 columns wide plot.
* `noprogress` (`false`): hide progress meter
* `kwargs...`: additional keywords arguments are passed directly to
  [`UnicodePlots.histogram`](https://github.com/JuliaPlots/UnicodePlots.jl).
"""
function plot_tmrcas(arg::Arg;
                     width = 76,
                     npoints = nrecombinations(arg) + 1,
                     noprogress = false,
                     kwargs...)
    lastidx = nrecombinations(arg) + 1
    stride = max(1, div(lastidx, npoints - 1, RoundDown) - 1)
    idx = StepRange(1, stride, lastidx)

    grid = [0.; collect(breakpoints(arg))[idx]]
    times = similar(grid, Float64)

    @showprogress enabled = !noprogress Threads.@threads :greedy for k ∈ eachindex(grid)
        times[k] = tmrca(arg, leaves(arg), grid[k])
    end

    plt = stairs(grid, times,
                 width = width,
                 margin = 0,
                 title = "TMRCAs",
                 xlabel = "Position";
                 kwargs...)

    plt
end

#          +----------------------------------------------------------+
#          |                      Other methods                       |
#          +----------------------------------------------------------+

"""
    $(FUNCTIONNAME)(arg, s, d)
    $(FUNCTIONNAME)(arg, e)

Return the parent of `d` that is not `s` for a recombination vertex `d`. If `d`
is not a recombination vertex, returns `s`. Can also take an edge as argument.

# Methods
$(METHODLIST)

--*Internal*--
"""
function otherdad end

function otherdad(arg, s, d)
    ret = s
    for dad ∈ dads(arg, d)
        ret == dad && continue
        ret = dad
        break
    end

    ret
end

otherdad(arg, e) = otherdad(arg, src(e), dst(e))

include("Mmn.jl")
include("Recombination.jl")
include("Algebra.jl")

function validate(arg::Arg; check_mutations = true)
    flag = true
    n = nleaves(arg)

    ## General properties ##
    if check_mutations && nmutations(arg) != nmarkers(arg)
        @info "Number of mutations not equal to the number of markers"
        flag = false
    end

    ## Leaves ##
    for v ∈ leaves(arg)
        if (!iszero ∘ length)(children(arg, v))
            @info "Leaf with children" v
            flag = false
        end

        if (!isone ∘ length)(dads(arg, v))
            @info "Leaf with number of parents not equal to 1" v
            flag = false
        end

        if ancestral_intervals(arg, v) != AIsType([Ω(0, ∞)])
            @info "Ancestral interval of leaf's parental edge not equal to [0, ∞)" v
            flag = false
        end
    end

    ## Non-root coalescence vertices ##
    for v ∈ Iterators.flatten((range(n + 1, 2n - 1), range(2n + 1, nv(arg), step = 2)))
        v = VertexType(v)
        v == mrca(arg) && continue

        if length(children(arg, v)) != 2
            @info "Coalescence vertex with number of children not equal to 2" v
            flag = false
        end

        if (!isone ∘ length)(dads(arg, v))
            @info "Coalescence vertex with number of parents not equal to 1" v
            flag = false
        end

        ai_children = mapreduce(x -> ancestral_intervals(arg, Edge(v => x)), ∪, children(arg, v))
        if ancestral_intervals(arg, Edge(dad(arg, v) => v)) != ai_children
            msg = "Coalescence vertex whose parental edge's ancestral interval" *
                " is not equal to the union of the ancestral intervals of its" *
                " children's edge."
            @info msg v
            flag = false
        end
    end

    ## Recombination vertices ##
    for v ∈ range(2n, nv(arg), step = 2)
        if (!isone ∘ length)(children(arg, v))
            @info "Recombination vertex with number of children not equal to 1" v
            flag = false
        end

        if length(dads(arg, v)) != 2
            @info "Recombination vertex with number of parents not equal to 2" v
            flag = false
        end

        # ai_child = ancestral_intervals(arg, Edge(v => child(arg, v)))
        # ai_right = ancestral_intervals(arg, Edge(rightdad(arg, v) => v))
        # ai_left = ancestral_intervals(arg, Edge(leftdad(arg ,v) => v))
        # if ai_child != ai_left ∪ ai_right
        #     msg = "Recombination vertex whose children edge's ancestral" *
        #         " interval is not equal to the union of it parental edges's" *
        #         " ancestral intervals."
        #     @info msg v ai_child ai_left ai_right
        #     flag = false
        # end
    end

    ## Internal vertices ##
    for v ∈ ivertices(arg)
        v = VertexType(v)
        h = sequence(arg, v) & ancestral_mask(arg, v)

        ref = mapreduce(&, children(arg, v)) do _child
            sequence(arg, _child) | .~ancestral_mask(arg, Edge(v => _child))
        end
        ref &= ancestral_mask(arg, v)

        if h != ref
            @info "Inconsistent sequence" v sequence(arg, v) ref
            flag = false
        end

        if latitude(arg, v) < maximum(u -> latitude(arg, u), children(arg, v))
            @info "Inconsistent latitudes" v
            flag = false
        end
    end

    ## Edges ##
    # for e ∈ edges(arg)
    #     if isempty(ancestral_intervals(arg, e))
    #         @info "Empty ancestral interval" e
    #         flag = false
    #     end
    # end

    flag
end
