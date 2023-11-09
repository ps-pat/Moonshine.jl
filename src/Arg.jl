using Graphs:
    SimpleDiGraph,
    SimpleEdge,
    AbstractSimpleGraph,
    degree,
    src, dst,
    Edge

import Graphs:
    edges,
    edgetype,
    has_edge,
    has_vertex,
    inneighbors,
    ne,
    nv,
    outneighbors,
    vertices,
    is_directed,
    add_vertex!,
    add_edge!,
    rem_edge!

import Base:
    eltype,
    show,
    union, intersect, in
using Base: @invoke

using Random:
    AbstractRNG,
    GLOBAL_RNG,
    shuffle!

using GraphMakie:
    graphplot,
    register_interaction!,
    deregister_interaction!,
    NodeDrag,
    NodeHoverHighlight
using LayeredLayouts
using GeometryBasics: Point

using Distributions:
    Bernoulli,
    Exponential,
    Uniform,
    DiscreteUniform,
    Categorical,
    truncated,
    logccdf,
    logpdf

using IntervalSets

using StaticArrays:
    SVector, MVector,
    SA

const VertexType = Int
const EdgeType = SimpleEdge{VertexType}
const Ω = Interval{:closed, :open, Float64}
const ∞ = Inf

Ω(x) = Ω(x, ∞)

export Ω

######################
# ArgCore definition #
######################

struct ArgCore{T}
    graph::SimpleDiGraph{VertexType}
    latitudes::Vector{Float64}
    sequences::Vector{Sequence{T}}
    nleaves::Int
    ancestral_interval::Dict{EdgeType, Set{Ω}}

    positions::Vector{Float64}
    recombinations_position::Vector{Float64}
    seq_length::Float64
    eff_popsize::Float64
    μ_loc::Float64
    ρ_loc::Float64
end

## TODO: Check performance impact of sizehint!.
function ArgCore{T}(leaves::AbstractVector{Sequence{T}};
                    positions = [],
                    seq_length = 1.0,
                    effective_popsize = 1.0,
                    μ_loc = 0.0,
                    ρ_loc = 0.0) where T
    n = length(leaves)
    latitudes = sizehint!(Float64[], 10n)

    sequences = sizehint!(similar(leaves, Sequence{T}, n), 10n)
    for (k, leaf) ∈ enumerate(leaves)
        sequences[k] = leaf
    end

    nmarkers = (length ∘ first)(leaves)
    if length(positions) != nmarkers || !issorted(positions)
        @info((isempty(positions) ? "A" : "Invalid positions: a") *
            "ssuming equally spaced markers")
        positions = isone(nmarkers) ?
            [0.0] : (collect ∘ range)(0, 1, length = nmarkers)
    end
    if minimum(positions) < 0
        @info("First position is < 0: shifting positions right")
        positions .-= minimum(positions)
    end
    if maximum(positions) > 1
        @info("Last position is > 0: scaling positions")
        positions ./= maximum(positions)
    end

    ArgCore{T}(SimpleDiGraph(n), latitudes, sequences, n, Dict(),
               positions, [], seq_length, effective_popsize, μ_loc, ρ_loc)
end

function ArgCore{T}(rng::AbstractRNG,
                    nmin::Integer, minlength::Integer,
                    nmax::Integer = 0, maxlength::Integer = 0;
                    genpars...) where T
    n = iszero(nmax) ? nmin : rand(rng, nmin:nmax)
    nmarkers = iszero(maxlength) ? minlength : rand(rng, minlength:maxlength)

    ArgCore{T}([Sequence{T}(rng, nmarkers) for _ ∈ 1:n]; genpars...)
end

function ArgCore{T}(nmin::Integer, minlength::Integer,
                    nmax::Integer = 0, maxlength::Integer = 0;
                    genpars...) where T
    ArgCore{T}(GLOBAL_RNG, nmin, minlength, nmax, maxlength; genpars...)
end

##################
# Arg definition #
##################

mutable struct Arg{T} <: AbstractSimpleGraph{VertexType}
    core::ArgCore{T}
    logprob::BigFloat
end
export Arg

Arg(leaves::AbstractVector{Sequence{T}}; genpars...) where T =
    Arg(ArgCore{T}(leaves; genpars...), zero(BigFloat))

function Arg{T}(rng::AbstractRNG,
                nmin::Integer, minlength::Integer,
                nmax::Integer = 0, maxlength::Integer = 0;
                genpars...) where T
    Arg(ArgCore{T}(rng, nmin, minlength, nmax, maxlength; genpars...),
        zero(BigFloat))
end
function Arg(rng::AbstractRNG,
             nmin::Integer, minlength::Integer,
             nmax::Integer = 0, maxlength::Integer = 0;
             genpars...)
    Arg{UInt}(rng, nmin, minlength, nmax, maxlength; genpars...)
end

function Arg{T}(nmin::Integer, minlength::Integer,
                nmax::Integer = 0, maxlength::Integer = 0;
                genpars...) where T
    Arg{T}(GLOBAL_RNG, nmin, minlength, nmax, maxlength; genpars...)
end

function Arg(nmin::Integer, minlength::Integer,
             nmax::Integer = 0, maxlength::Integer = 0;
             genpars...)
    Arg{UInt}(nmin, minlength, nmax, maxlength; genpars...)
end

## AbstractGraph Interface.
## See https://juliagraphs.org/Graphs.jl/stable/ecosystem/interface/
for fun ∈ [:edges, :vertices, :ne, :nv]
    @eval function $fun(arg::Arg)
        $fun(arg.core.graph)
    end
end

for (fun, ret) ∈ Dict(:eltype => VertexType,
                      :edgetype => EdgeType,
                      :is_directed => :true)
    @eval begin
        @generated $fun(::Type{Arg}) = $ret
        @generated $fun(::Arg) = $fun(Arg)
    end
end

let arg_edge = Expr(:(::), :e, EdgeType),
    arg_vertex = Expr(:(::), :v, VertexType)
    for (fun, a) ∈ Dict(:has_edge => arg_edge,
                        :has_vertex => arg_vertex,
                        :inneighbors => arg_vertex,
                        :outneighbors => arg_vertex)
        varname = first(a.args)
        @eval function $fun(arg::Arg, $a)
            $fun(arg.core.graph, $varname)
        end
    end
end

function add_vertex!(arg::Arg, seq = Sequence(), lat = ∞)
    push!(arg.core.latitudes, lat)
    push!(arg.core.sequences, seq)
    add_vertex!(arg.core.graph)
end

function add_edge!(arg::Arg, e, ints::Set{Ω})
    arg.core.ancestral_interval[e] = ints
    add_edge!(arg.core.graph, e)
end

add_edge!(arg::Arg, e, int::Ω) =
    add_edge!(arg, e, int ∩ ancestral_intervals(arg, dst(e)))

add_edge!(arg::Arg, e) = add_edge!(arg, e, Ω(0, ∞))

function rem_edge!(arg::Arg, e)
    delete!(arg.core.ancestral_interval, e)
    rem_edge!(arg.core.graph, e)
end

## Simple methods.
for field ∈ [:(:nleaves), :(:sequences), :(:latitudes),
             :(:seq_length), :(:eff_popsize), :(:positions),
             :(:recombinations_position)]
    fun_name = eval(field)
    @eval begin
        export $fun_name

        function $fun_name(arg)
            getfield(arg.core, $field)

        end
    end
end

nrecombinations(arg) = length(arg.core.recombinations_position)

leaves(arg) = Base.OneTo(arg.core.nleaves)
isleaf(arg, v) = v <= nleaves(arg)

function ivertices(arg)
    n = nleaves(arg)
    range(n + 1, length = n - 1)
end
nivertices(arg) = nleaves(arg) - 1

"""
    nmarkers(arg)

Number of markers of the sequences in an ARG.
"""
nmarkers(arg) = (length ∘ first)(arg.core.sequences)
export nmarkers

sequences(arg, v::VertexType) = sequences(arg)[v]
sequences(arg, e::EdgeType) = [sequences(arg, src(e)), sequences(arg, dst(e))]

sequences(arg, vs::AbstractArray{VertexType}) =
    Iterators.map(v -> sequences(arg, v), vs)

latitude(arg, v) =
    isleaf(arg, v) ? zero(Float64) : arg.core.latitudes[v - nleaves(arg)]
export latitude

edge_length(arg, e) = latitude(arg, src(e)) - latitude(arg, dst(e))

export mrca, tmrca
mrca(arg) = isempty(arg.core.latitudes) ?
    zero(Int) : argmax(arg.core.latitudes) + nleaves(arg)
tmrca(arg) = isempty(arg.core.latitudes) ?
    zero(Float64) : maximum(arg.core.latitudes)

"""
    children(arg, v, int = Ω(0, ∞))
    parents(arg, v, int = Ω(0, ∞))
    siblings(arg, v, int = Ω(0, ∞))

Return the children/parents/siblings of a vertex bearing ancestral material for
a given interval.
"""
function children end,
function parents end,
function siblings end
export children, parents, siblings

children(arg, v) = outneighbors(arg, v)

parents(arg, v) = inneighbors(arg, v)

function siblings(arg, v, int = Ω(0, ∞))
    ret = mapreduce(p -> children(arg, p, int), ∪, parents(arg, v, int),
                    init = eltype(arg)[])

    v_idx = findfirst(==(v), ret)
    isnothing(v_idx) && return ret

    @inbounds ret[setdiff(eachindex(ret), v_idx)]
end

for (fun, edge_fun) ∈ Dict(:children => (x, y) -> Edge(x, y),
                           :parents => (x, y) -> Edge(y, x))
    @eval $fun(arg, v, int) =
        filter(x -> !isempty(ancestral_intervals(arg, $edge_fun(v, x)) ∩ int),
               $fun(arg, v))
end

"""
    idxtopos(arg, idx)

Return the position of the marker given its index.
"""
idxtopos(arg, idx) = arg.core.positions[idx]

"""
    postoidx(arg, pos)

Return the largest marker's index that is at a position lesser than the one
given.
"""
function postoidx(arg, pos)
    @inbounds for (k, p) ∈ enumerate(arg.core.positions)
        p >= pos && return k - 1
    end

    nmarkers(arg)
end

"""
    ancestral_intervals(arg, e)
    ancestral_intervals(arg, σ, δ)
    ancestral_intervals(arg, v)

Compute the interval for which an edge or a vertex is ancestral.
"""
function ancestral_intervals end
export ancestral_intervals

ancestral_intervals(arg, e::EdgeType) = arg.core.ancestral_interval[e]
ancestral_intervals(arg, σ, δ) = ancestral_intervals(arg, Edge(σ, δ))

function ancestral_intervals(arg, v::VertexType)
    isleaf(arg, v) && return Set([Ω(0, ∞)])


    mapreduce(child -> ancestral_intervals(arg, v, child), ∪, children(arg, v))
end

@generated blocksize(::Arg{T}) where T = blocksize(Sequence{T})

export nmutations
nmutations(arg, e) = count_ones(xor(sequences(arg, e)...))

nmutations(arg) = mapreduce(e -> nmutations(arg, e), +, edges(arg),
                            init = zero(Int))

function branchlength_tree(arg)
    lats = view(latitudes(arg), range(1, length = nleaves(arg) - 1))
    sum(lats, init = last(lats))
end

export branchlength
function branchlength(arg)
    ν = nrecombinations(arg)
    lats = view(latitudes(arg), range(nleaves(arg), length = 2ν))
    rec_branchlength = mapreduce(vs -> last(vs) - first(vs), +,
                                 Iterators.partition(lats, 2), init = 0.0)

    branchlength_tree(arg) + rec_branchlength
end

for (fun, par) ∈ Dict(:mut_rate => :(:μ_loc), :rec_rate => :(:ρ_loc))
    @eval begin
        export $fun

        function $fun(arg, scaled = true)
            rate_loc = getfield(arg.core, $par)
            mul = scaled ? 4eff_popsize(arg) : one(rate_loc)

            mul * rate_loc
        end
    end
end

export subarg
function subarg(arg, int)
    newarg = deepcopy(arg)

    for (e, int_set) ∈ pairs(newarg.core.ancestral_interval)
        isempty(int ∩ int_set) || continue

        rem_edge!(newarg, e)
    end

    newarg
end

"""
    maxdepth(arg, v, pos, depth = 0)

Compute the depth of a vertex, that is the number of egdes between it and the
arg's mrca.
"""
function maxdepth(arg, v, int, depth = 0)
    _parents = parents(arg, v, int)
    isempty(_parents) && return depth

    mapreduce(parent -> maxdepth(arg, parent, int, depth + 1), max, _parents)
end


##############################
# Plotting & pretty printing #
##############################

ArgLayout(int) = function(arg)
    ## Force leaves to be in last layer.
    lastlayer = maximum(v -> maxdepth(arg, v, int), leaves(arg)) + 1

    xs, ys, _ = solve_positions(Zarate(),
                                arg.core.graph,
                                force_layer = Pair.(leaves(arg), lastlayer))

    ## Rotate by -π/2.
    Point.(zip(ys, -xs))
end

function argplot end
export argplot

function argplot(arg, int::Ω;
                 arrow_show = true,
                 wild_color = :blue,
                 derived_color = :red,
                 attributes...)
    newarg = deepcopy(arg)
    n = nv(newarg)

    ## Remove non ancestral edges.
    for (e, eint) ∈ newarg.core.ancestral_interval
        isempty(int ∩ eint) || continue
        rem_edge!(newarg, e)
    end

    nedges = ne(newarg)

    edgecolor = fill(:gray, nedges)

    mask = ancestral_mask(arg, int)
    nodecolor = map(sequences(newarg)) do σ
        any(σ & mask) ? derived_color : wild_color
    end

    node_labels = string.(1:n)

    node_sizes = fill(25.0, n)

    p = graphplot(newarg,
                  layout = ArgLayout(int),
                  nlabels = node_labels,
                  nlabels_distance = 10,
                  node_size = node_sizes,
                  node_color = nodecolor,
                  edge_color = edgecolor,
                  edge_width = fill(3, nedges),
                  elabels = nothing,
                  arrow_show = arrow_show,
                  attributes...)


    deregister_interaction!(p.axis, :rectanglezoom)
    register_interaction!(p.axis, :nodedrag, NodeDrag(p.plot))
    register_interaction!(p.axis, :nodehover, NodeHoverHighlight(p.plot))

    p
end

function show(io::IO, ::MIME"text/plain", arg::Arg)
    println(io, "Ancestral recombination graph:")
    println(io, nleaves(arg), " leaves, ",
            nmarkers(arg), " markers")
    print(io, "tMRCA: ", tmrca(arg))
end

function show(io::IO, arg::Arg)
    print(io, "ARG(")
    print(io, nleaves(arg))
    print(io, " leaves, ")
    print(io, nmarkers(arg))
    print(io, " markers)")
end

#####################
# Tree construction #
#####################

function coalesce!(rng::AbstractRNG, arg, vertices, nlive)
    ## Select coalescing pair.
    arg.logprob += (log ∘ inv ∘ binomial)(length(vertices), 2)
    shuffle!(rng, vertices)

    ## Sample a latitude for the coalescence event.
    Δdist = Exponential(inv(nlive - 1))
    Δ = rand(rng, Δdist)
    arg.logprob += logpdf(Δdist, Δ)
    newlat = latitude(arg, nv(arg)) + Δ

    ## Perform the coalescence.
    parent = coalesce!(arg, pop!(vertices), pop!(vertices), newlat)

    ## Update live vertices.
    push!(vertices, parent)

    @debug("$(first(_children)) and $(last(_children)) \
           coalesced into $parent at $(latitude(arg, parent))")
end

function coalesce!(arg::Arg, v1, v2, lat)
    newseq = sequences(arg, v1) & sequences(arg, v2)

    add_vertex!(arg, newseq, lat) ||
        @error "Could not add a vertex to ARG"
    parent = nv(arg)

    for child ∈ SA[v1, v2]
        e = Edge(parent, child)
        add_edge!(arg, e)
        arg.core.ancestral_interval[e] = Set([Ω(0, ∞)])
    end

    parent
end

"""
    buildtree!(rng = GLOBAL_RNG, arg, idx = 1)

Build a marginal tree consistent with a specific index.
"""
function buildtree! end
export buildtree!

function buildtree!(rng::AbstractRNG, arg::Arg, idx = 1)
    nv(arg) > nleaves(arg) && @warn "ARG contains non-leaf vertices"

    derived = findall(σ -> getindex(σ, idx), sequences(arg))
    wild = setdiff(range(1, length = length(sequences(arg))), derived)

    nlive = nleaves(arg)
    nlive_derived, nlive_wild = length(derived), length(wild)

    ## While there are more than 1 live derived vertex, we build two
    ## distinct subtrees.
    while nlive_derived > 1
        ## Determine the type of coalescence.
        wild_weight = nlive_wild * (nlive_wild - 1)
        derived_weight = nlive_derived * (nlive_derived - 1)
        derived_prob = derived_weight / (wild_weight + derived_weight)

        coalescence_dist = Bernoulli(derived_prob)
        derived_coalescence = rand(rng, coalescence_dist)
        cvertices = derived_coalescence ? derived : wild
        arg.logprob += logpdf(coalescence_dist, derived_coalescence)

        ## Perform the coalescence.
        coalesce!(rng, arg, cvertices, nlive)

        ## Update the number of live nodes.
        nlive -= 1
        nlive_derived -= derived_coalescence
        nlive_wild -= !derived_coalescence
    end

    cvertices = [derived; wild]
    while nlive > 1
        coalesce!(rng, arg, cvertices, nlive)
        nlive -= 1
    end

    arg
end
buildtree!(arg::Arg, idx = 1) = buildtree!(GLOBAL_RNG, arg, idx)

export quotient_leaves
function quotient_leaves(arg::Arg{T}) where T
    RetEltype = Vector{eltype(arg)}
    d = Dict{Sequence{T}, RetEltype}()

    for (v, σ) ∈ enumerate(sequences(arg, leaves(arg)))
        haskey(d, σ) || setindex!(d, RetEltype(), σ)
        push!(d[σ], v)
    end

    (collect ∘ values)(d)
end

#######
# MMN #
#######

"""
    mmn_blocks_idx(arg, e, int = Ω(0, ∞))

Compute the indices of blocks to be used in a pass of the mmn algorithm.
"""
function mmn_blocks_idx end

function mmn_blocks_idx(arg::Arg{T}, e::EdgeType, int = Ω(0, ∞)) where T
    target_intervals = ancestral_intervals(arg, e) ∩ int

    n = nmarkers(arg)
    bs = blocksize(arg)

    (sort ∘ mapreduce)(union, target_intervals, init = Int[]) do i
        _endpoints = endpoints(i)

        postoidx_arg = Fix1(postoidx, arg)
        lbound = (postoidx_arg ∘ first)(_endpoints)
        rbound = (postoidx_arg ∘ prevfloat ∘ last)(_endpoints)

        blockidx(n, bs, lbound), blockidx(n, bs, rbound)
    end
end

mmn_blocks_idx(arg, σ::VertexType, δ::VertexType, int = Ω(0, ∞)) =
    mmn_blocks_idx(arg, EdgeType(σ, δ), int)

function first_inconsistent_position end
export first_inconsistent_position

function first_inconsistent_position(arg, int::Ω)
    n = nmarkers(arg)
    _mrca = mrca(arg)
    (iszero(_mrca) || iszero(n)) && return ∞, EdgeType[]
    bs = blocksize(arg)

    ## Vertices to xor are stored in acc.
    acc = Set(_mrca)

    ## Contains mutation edges.
    mutation_edges = Dict{Int, Vector{EdgeType}}()

    while !isempty(acc)
        vertex = pop!(acc)

        for e ∈ map(Fix1(EdgeType, vertex), children(arg, vertex))
            blocks_idx = mmn_blocks_idx(arg, e, int)

            ## Non ancestral edge.
            isempty(blocks_idx) && continue

            ## Add child to acc.
            c = dst(e)
            isleaf(arg, c) || push!(acc, c)

            xored_sequences =
                (sequences(arg, src(e)).data[blocks_idx] .⊻
                sequences(arg, dst(e)).data[blocks_idx]) .&
                ancestral_mask(arg, e).data[blocks_idx]
            for (block_idx, block) ∈ zip(blocks_idx, xored_sequences)
                first_set_bit = leading_zeros(block) + 1
                if first_set_bit <= bs
                    idx = actualpos(n, bs, block_idx, first_set_bit)
                    haskey(mutation_edges, idx) || (mutation_edges[idx] = EdgeType[])
                    push!(mutation_edges[idx], e)
                    break
                end
            end
        end
    end

    ## Compute the first inconsistent position.
    inconsistent_idx =
        minimum(p -> (length ∘ last)(p) > 1 ? first(p) : typemax(Int),
                mutation_edges)

    idxtopos(arg, inconsistent_idx), Set(mutation_edges[inconsistent_idx])
end

first_inconsistent_position(arg, lbound::Real, ubound::Real = ∞) =
    first_inconsistent_position(arg, Ω(lbound, ubound))

####################
# ARG construction #
####################

"""
    recombine!(arg, redge, cedge, breakpoint, rlat, clat)
    recombine!([rng = GLOBAL_RNG], arg, redge, other_mutation_edges, breakpoint, α = ∞)

Create recombination/recoalescence events on an ARG. The coalescence edge is
returned.
"""
function recombine! end
export recombine!

function recombine!(arg::Arg, redge, cedge, breakpoint, rlat, clat)
    @debug "Recombination $(nrecombinations(arg) + 1)" redge = redge cedge = cedge

    ## Recombination and recoalescence vertices.
    rvertex, cvertex = nv(arg) .+ (1:2)

    ## Record ancestral intervals
    old_edges = SA[redge, cedge]
    ints_redge, ints_cedge = Fix1(ancestral_intervals, arg).(old_edges)

    ## Remove recombination and recoalescence edges.
    Fix1(rem_edge!, arg).(old_edges)

    ## Compute sequences.
    rseq = deepcopy(sequences(arg, dst(redge)))
    newseqs = SVector{2, eltype(sequences(arg))}(
        rseq,
        rseq & sequences(arg, dst(cedge)))

    ## Add recombination and recoalescences vertices.
    for (seq, lat) ∈ zip(newseqs, SA[rlat, clat])
        add_vertex!(arg, seq, lat)
    end

    ## New edges. `dst(cedge)` should be equal to `mrca(arg)` if the
    ## recoalescence happens above `arg` current MRCA. Additionaly,
    ## `src(cedge)` should be lesser than or equal to 0.
    add_edge!(arg, Edge(rvertex, dst(redge)), ints_redge) # Recombination vertex out
    add_edge!(arg, Edge(src(redge), rvertex), Ω(0, breakpoint)) # Old tree
    add_edge!(arg, Edge(cvertex, rvertex), Ω(breakpoint, ∞)) # New tree
    add_edge!(arg, Edge(cvertex, dst(cedge)), ints_cedge) # Old tree
    src(cedge) > 0 && add_edge!(arg, Edge(src(cedge), cvertex))

    push!(arg.core.recombinations_position, breakpoint)

    update_upstream!(arg, rvertex)

    cedge
end

function recombine!(rng::AbstractRNG, arg, redge,
                    other_mutation_edges, breakpoint; α = ∞)
    ## Sample a recombination latitude uniformly on the recombination edge.
    rlat_dist = Uniform(latitude(arg, dst(redge)), latitude(arg, src(redge)))
    rlat = rand(rng, rlat_dist)
    arg.logprob += logccdf(rlat_dist, rlat)

    ## Sample an exponential recoalescence latitude.
    derived_recombination =
        sequences(arg, dst(redge))[postoidx(arg, breakpoint) + 1]

    clat_ubound = derived_recombination ?
        maximum((Fix1(latitude, arg) ∘ src).(other_mutation_edges)) :
        Inf

    Δlat_dist = truncated(Exponential(nblive(arg, rlat, breakpoint)),
                          upper = clat_ubound - rlat)
    Δlat = rand(rng, Δlat_dist)
    arg.logprob += logccdf(Δlat_dist, Δlat)
    clat = rlat + Δlat

    if derived_recombination
        possible_cedges = possible_cedges_derived(arg, clat, breakpoint,
                                                  dst.(other_mutation_edges),
                                                  α)
        if isempty(possible_cedges)
            @info "No possible coalescence edges" rlat = rlat clat = clat
            return arg
        end
    end

    arg.logprob -= log(length(possible_cedges))
    recombine!(arg, redge, rand(rng, possible_cedges), breakpoint, rlat, clat)
end

recombine!(arg::Arg, redge, other_mutation_edges, breakpoint; α = ∞) =
    recombine!(GLOBAL_RNG, arg, redge, other_mutation_edges, breakpoint, α = α)

"""
    nblive(arg, latitude, pos)

Compute the number of live vertices.
"""
function nblive(arg, lat, pos)
    lat > tmrca(arg) && return one(Int)

    (length ∘ filter)(arg.core.ancestral_interval) do π
        e, int = π

        pos ∈ int || return false
        latitude(arg, dst(e)) <= lat <= latitude(arg, src(e)) || return false

        true
    end
end
export nblive

"""
    possible_cedges_derived(arg, clat, breakpoint, seeds = [mrca(arg)], α = ∞)

Compute possible coalescence edges for a derived recoalescence. Two conditions
must be met for an edge `e` to be admissible:
1. `e` must be ancestral for a position in [`breakpoint - α`, `breakpoint`)
2. `clat` ∈ [`dst(e)`, `src(e)`].

The search for coalescence edges goes downward from the vertices in `seeds`.
"""
function possible_cedges_derived(arg, clat, breakpoint,
                                 seeds = [mrca(arg)], α = ∞)
    idx = postoidx(arg, breakpoint) + 1
    cedges = Set{EdgeType}()
    ancestral_int = Ω(clamp(breakpoint - α, 0, breakpoint), breakpoint)

    acc = Set(seeds)
    while !isempty(acc)
        v = pop!(acc)

        if sequences(arg, v)[idx]
            ## Vertex is derived.
            for parent ∈ parents(arg, v, ancestral_int)
                if latitude(arg, v) ≤ clat ≤ latitude(arg, parent)
                    push!(cedges, Edge(parent, v))
                elseif latitude(arg, v) > clat
                    for child ∈ children(arg, v)
                        push!(acc, child)
                    end
                end
            end
        else
            ## Vertex is wild.
            for child ∈ children(arg, v)
                push!(acc, child)
            end
        end
    end
    cedges
end

# function possible_cedges_wild(arg, clat, breakpoint, α = ∞)
#     acc = Set(mrca(arg))
# end

## TODO: allow wild recombinations.
"""
    make_consistent!([rng = GLOBAL_RNG], arg, pos, mutation_edges; α = ∞)

Modifies an ARG with a sequence of recombination/recoalescence events in order
to make it consistent with a given position.
"""
function make_consistent! end
export make_consistent!

function make_consistent!(rng, arg, pos, mutation_edges; α = ∞)
    while length(mutation_edges) > 1
        ## Sample recombination breakpoint.
        bp_lbound = isempty(recombinations_position(arg)) ? 0 :
            max(idxtopos(arg, max(1, postoidx(arg, pos) - 1)),
                last(recombinations_position(arg)))
        Δbp_dist = truncated(Exponential(rec_rate(arg)), upper = pos - bp_lbound)
        Δbp = rand(rng, Δbp_dist)
        arg.logprob += logccdf(Δbp_dist, Δbp)
        breakpoint = bp_lbound + Δbp

        ## Sample a recombination edge.
        redge = sample_edge(rng, arg, collect(mutation_edges))
        delete!(mutation_edges, redge)

        ## Recombination/recoalescence
        cedge = recombine!(rng, arg, redge, mutation_edges, breakpoint, α = α)
        delete!(mutation_edges, cedge)

        ## Update the set of mutation edges.
        push!(mutation_edges, upstream_mutation_edge(arg, pos, breakpoint))
    end
end

function sample_edge(rng, arg, edges)
    edges_weights = map(Fix1(edge_length, arg), edges)
    redge_dist = Categorical(edges_weights / sum(edges_weights))
    redge_idx = rand(rng, redge_dist)
    arg.logprob += logpdf(redge_dist, redge_idx)

    edges[redge_idx]
end

function upstream_mutation_edge(arg, pos, breakpoint)
    idx = postoidx(arg, pos)
    child = nv(arg)
    parent = (first ∘ parents)(arg, child, Ω(breakpoint, ∞))

    while sequences(arg, parent)[idx]
        child = parent
        parent = (first ∘ parents)(arg, child, Ω(breakpoint, ∞))
    end

    Edge(parent, child)
end

"""
    ancestral_mask(arg, x)

Mask non ancestral positions to 0.
"""
ancestral_mask(arg, x::Set{Ω}) = mapreduce(|, x) do int
    lpos, rpos = endpoints(int)

    lidx = 1
    while lpos > arg.core.positions[lidx]
        lidx += 1
    end

    ridx = postoidx(arg, rpos)

    andmask(nmarkers(arg), UnitRange(lidx, ridx))
end

ancestral_mask(arg, x::Ω) = ancestral_mask(arg, Set([x]))

ancestral_mask(arg, x::EdgeType) =
    ancestral_mask(arg, ancestral_intervals(arg, x))

update_sequence!(arg, v) =
    arg.core.sequences[v] = mapreduce(&, children(arg, v)) do child
        mask = ancestral_mask(arg, Edge(v, child))
        sequences(arg, child) | ~mask
    end

update_intervals!(arg, e, ints_ubound) =
    arg.core.ancestral_interval[e] =
    (simplify! ∘ mapreduce)(Fix1(∩, ints_ubound), ∪,
                            arg.core.ancestral_interval[e])

function update_upstream!(arg, v)
    vstack = Set{eltype(arg)}(parents(arg, v))
    while !isempty(vstack)
        v = pop!(vstack)

        old_sequence = sequences(arg, v)
        if old_sequence != update_sequence!(arg, v)
            _parents = parents(arg, v)
            isempty(_parents) || push!(vstack, _parents...)
        end

        for parent ∈ parents(arg, v)
            e = Edge(parent, v)
            ints_ubound = ancestral_intervals(arg, v)

            old_ints = ancestral_intervals(arg, e)
            old_ints == update_intervals!(arg, e, ints_ubound) ||
                push!(vstack, parent)
        end
    end

    arg
end

export highest_under
"""
    highest_under(arg, lat)

Find the set S of vertices of `arg` such that:
1. Every vertices of S has a latitude <= to `lat`;
2. Every leaves of `arg` has an ancestor in S;
3. The cardinality of S is minimal.
"""
function highest_under(arg, lat)
    ret = CheapStack(eltype(arg), nleaves(arg))

    acc = CheapStack(eltype(arg), nleaves(arg))
    push!(acc, mrca(arg))
    while !isempty(acc)
        v = pop!(acc)
        if latitude(arg, v) <= lat
            push!(ret, v)
            continue
        end
        push!(acc, children(arg, v)...)
    end

    ret
end
