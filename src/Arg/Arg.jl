using Graphs

import Graphs: add_vertices!, add_edge!, rem_edge!

using Random

using StatsBase: samplepair, ProbabilityWeights, fit, Histogram

using SparseArrays

using Combinatorics: combinations

using Distributions

using StaticArrays: @SVector

using UnicodePlots: heatmap, label!

##################
# Arg Definition #
##################

export Arg
struct Arg <: AbstractGenealogy
    graph::SimpleDiGraph{VertexType}
    latitudes::Vector{Float64}
    recombination_mask::Vector{AIsType}
    mrca::Base.RefValue{VertexType}
    sequences::Vector{Sequence}
    ancestral_intervals::Dict{Edge{VertexType}, AIsType}
    sample::Sample
    logprob::Base.RefValue{Float64x2}
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
        Ref(prob(tree, logscale = true)))
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

rem_edge!(arg::Arg, e) = rem_edge!(graph(arg), e)

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

function ancestral_intervals!(ωs, arg::Arg, e::Edge; wipe = true, simplify = true)
    wipe && empty!(ωs)

    union!(ωs, ancestral_intervals(arg, e), simplify = simplify)

    ωs
end

ancestral_intervals(arg::Arg, e::Edge) =
    get!(() -> AIsType([Ω(0, ∞)]), arg.ancestral_intervals, e)

function ancestral_intervals!(ωs, arg::Arg, v::VertexType; wipe = true)
    wipe && empty!(ωs)
    isleaf(arg, v) && return push!(ωs, Ω(0, ∞))

    for child ∈ children(arg, v)
        ancestral_intervals!(ωs, arg, Edge(v => child), wipe = false, simplify = false)
    end

    simplify!(ωs)
end

ancestral_intervals(arg::Arg, v::VertexType) = ancestral_intervals!(AIsType(), arg, v)

ancestral_mask!(η, arg::Arg, e::Edge{VertexType}; wipe = true) =
    ancestral_mask!(η, sam(arg), ancestral_intervals(arg, e), wipe = wipe)

function ancestral_mask!(h, arg::Arg, v::VertexType;
                         buffer = default_buffer(), wipe = true)
    ## Compute number of Ωs
    len = sum(c -> (length ∘ ancestral_intervals)(arg, Edge(v => c)),
              children(arg, v))

    @no_escape buffer begin
        ωs = @alloc(Ω, len)
        k = 1
        for c ∈ children(arg, v)
            for ω ∈ ancestral_intervals(arg, Edge(v, c))
                ωs[k] = ω
                k += 1
            end
        end
        ancestral_mask!(h, sam(arg), ωs, wipe = wipe)
        nothing # Workaround, see Bumper.jl's issue #49
    end

    h
end

ancestral_mask(arg::Arg, x::Union{VertexType, Edge{VertexType}}) =
    ancestral_mask!(Sequence(falses(nmarkers(arg))), arg, x)

recidx(arg, v) = (v - 2(nleaves(arg) - 1)) ÷ 2

function ancestral_mask(e::Edge, arg)
    s, d = src(e), dst(e)

    inc = s > otherdad(arg, s, d)
    arg.recombination_mask[2recidx(arg, d) - 1 + inc]
end

export breakpoints
"""
    breakpoint(arg)

Recombination events' positions
"""
breakpoints(arg::Arg) =
    Iterators.map(rightendpoint, @view arg.recombination_mask[1:2:end])

export plot_breakpoints
"""
    plot_breakpoints(arg; kwargs...)

Heatmap of recombination events' positions. Additional keywords arguments are
passed directly to [`UnicodePlots.histogram`](@ref).

See also [`breakpoints`](@ref)
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

#          +----------------------------------------------------------+
#          |                      Other methods                       |
#          +----------------------------------------------------------+

"""
    otherdad(arg, s, d)
    otherdad(arg, e)

Return the parent of `d` that is not `s` for a recombination vertex `d`. If `d`
is not a recombination vertex, returns `s`. Can also take an edge as argument.
"""
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

#          +----------------------------------------------------------+
#          |                      Thévenin Tree                       |
#          +----------------------------------------------------------+

export thevenin!, thevenin
"""
    thevenin!(rng, tree, arg)
    thevenin(rng, arg)

Sample a Thévenin tree conditional on an ARG.
"""
function thevenin! end, function thevenin end

function _thevenin_impedance_helper!(tree, arg, v1_leaves, v2, C, Z2,
    edgesmap, estack, vqueue, visited)
    r = nrecombinations(arg)
    v2_leaves = descendants(tree, v2) ∩ leaves(tree)
    if isempty(v2_leaves)
        push!(v2_leaves, v2)
    end

    p = length(v1_leaves) + length(v2_leaves)
    C_view = @view C[:, range(1, r + p - 1)]
    impedance!(arg, v1_leaves, v2_leaves, C_view, Z2,
        edgesmap = edgesmap, estack = estack,
        vqueue = vqueue, visited = visited)
end

function thevenin!(tree, arg::Arg; ϵ = 1e-5)
    n = nleaves(arg)
    r = nrecombinations(arg)

    ## Initialize tree ##
    empty!(graph(tree).fadjlist)
    graph(tree).ne = 0
    add_vertices!(graph(tree), 2n - 1)

    resize!(latitudes(tree), n - 1)

    resize!(sequences(tree), 2n - 1)
    sequences(tree)[1:n] .= sequences(arg, leaves(arg))

    tree.logprob[] = 0

    # ## Build tree ##
    edgesid = edgesmap(arg)
    estack = Stack{Edge{VertexType}}(ceil(Int, log(nv(arg))))
    vqueue = Queue{VertexType}(ceil(Int, log(nv(arg))))
    visited = Set{VertexType}()

    ## To limit allocations, we make `C` large enough to handle the worst case
    ## scenario. We can pass an appropriate `view` to `impedance!`.
    C = _impedance_C(arg, n - 1, edgesid)
    Z2 = _impedance_Z(arg, edgesid, true)

    live = collect(leaves(tree))

    ## Distance matrix. Used for clustering
    dists = impedance_matrix(arg, estack, edgesid)
    for k ∈ 1:n
        dists[k, k] = ∞
    end

    ## Descendant leaves of live vertices.
    dleaves = Vector{Set{VertexType}}(undef, n)
    for k ∈ eachindex(dleaves)
        dleaves[k] = Set(k)
    end

    ## Impedance between live vertices and their descendant leaves
    live_impedances = zeros(Float64, n)

    for nlive ∈ range(length(live), 2, step = -1)
        ## Find the next coalescing pair ##
        ii, jj = Tuple(findfirst(!isinf, dists))
        @inbounds for ii_it ∈ 1:(nlive-1)
            for jj_it ∈ (ii_it+1):nlive
                dists[ii_it, jj_it] <= dists[ii, jj] || continue
                ii, jj = ii_it, jj_it
            end
        end

        ## Create coalescence event ##
        i, j = live[ii], live[jj]
        k = 2n + 1 - nlive

        add_edge!(tree, k, i)
        add_edge!(tree, k, j)

        newlat = 0.5 * (dists[ii, jj] +
            latitude(tree, i) - live_impedances[ii] +
            latitude(tree, j) - live_impedances[jj])
        Δlat = newlat - latitude(tree, k - 1)
        if Δlat < ϵ
            newlat = latitude(tree, k - 1) + ϵ
        end
        latitudes(tree)[k - n] = newlat

        ## Compute descendant leaves ##
        union!(dleaves[ii], dleaves[jj])

        ## Update `live_impedances` ##
        live_impedances[ii] = inv(
            inv(live_impedances[ii] + branchlength(tree, k, i)) +
            inv(live_impedances[jj] + branchlength(tree, k, j)))

        ## Update live vertices and impedance matrix ##
        ## Replace i by k.
        live[ii] = k

        ## Make the distance from any vertex to j infinite.
        @inbounds for ii_it ∈ 1:(jj-1)
            dists.data[jj, ii_it] = ∞
        end
        @inbounds for ii_it ∈ (jj+1):n
            dists.data[ii_it, jj] = ∞
        end

        ## Compute distances from k.
        leaves_ii = dleaves[ii]
        n_ii = length(leaves_ii)

        @inbounds for jj_it ∈ 1:(ii-1)
            isinf(dists.data[ii, jj_it]) && continue
            leaves_jj_it = dleaves[jj_it]
            n_jj_it = length(leaves_jj_it)

            C_view = view(C, :, 1:(r + n_ii + n_jj_it - 1))
            dists.data[ii, jj_it] = impedance!(arg, leaves_ii, leaves_jj_it, C_view, Z2,
                                               edgesmap = edgesid, estack = estack,
                                               vqueue = vqueue, visited = visited)
        end

        @inbounds for jj_it ∈ (ii+1):n
            isinf(dists.data[jj_it, ii]) && continue
            leaves_jj_it = dleaves[jj_it]
            n_jj_it = length(leaves_jj_it)

            C_view = view(C, :, 1:(r + n_ii + n_jj_it - 1))
            dists.data[jj_it, ii] = impedance!(arg, leaves_jj_it, leaves_ii, C_view, Z2,
                                               edgesmap = edgesid, estack = estack,
                                               vqueue = vqueue, visited = visited)
        end
    end

    tree
end

thevenin(arg::Arg) = thevenin!(Tree(sam(arg)), arg)

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
        ref = mapreduce(&, children(arg, v)) do _child
            sequence(arg, _child) | ~ancestral_mask(arg, Edge(v => _child))
        end
        ref &= ancestral_mask(arg, v)

        if sequence(arg, v) != ref
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
