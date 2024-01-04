using Graphs: SimpleEdge,
              AbstractSimpleGraph,
              diameter,
              indegree, outdegree

import Graphs: edges, vertices, ne, nv,
               eltype, edgetype, is_directed,
               has_edge, has_vertex, inneighbors, outneighbors

import GraphMakie: graphplot

using LayeredLayouts

using IntervalSets

using StaticArrays: SA

using GeometryBasics: Point

using Bumper

const VertexType = Int
const EdgeType = SimpleEdge{VertexType}
const ∞ = Inf

export Ω
const Ω = Interval{:closed,:open,Float64}
Ω(x) = Ω(x, ∞)

abstract type AbstractGenealogy <: AbstractSimpleGraph{VertexType} end

#############
# Interface #
#############

export graph
"""
    graph(genealogy)

Return the underlying graph of a genealogy.
"""
function graph end

"""
    describe(genealogy, long = true)

Return a string containing a long or short description of a genealogy.
Used internally by pretty printing functions.
"""
function describe end

export latitudes
"""
    latitudes(genealogy[, ivs])

Return the latitudes of (a subset of) the internal vertices of a genealogy.

See also [`latitude`](@ref) to get the latitude of a single vertex.
"""
function latitudes end

export latitude
"""
    latitude(genealogy, v)

Latitude of a vertex.

See also [`latitudes`](@ref) to get the latitudes of internal vertices stored
in a genealogy.
"""
function latitude end

export leaves, ivertices
"""
    leaves(genealogy)
    ivertices(genealogy_vertex)

An interable containing the leaves/internal vertices of a genealogy.

See also [`nleaves`](@ref) and [`nivertices`](@ref) for the number of leaves
and internal vertices.
"""
function leaves end,
function ivertices end

export sequences
"""
    sequences(genealogy[, e])
    sequences(genealogy[, vs])

Sequences of a genealogy. If an edge is specified, return the sequences
associated with the vertices incident to that edge. If an iterable of
vertices is specified, return the sequences associated with these vertices.

See also [`sequence`](@ref) to get the sequence associated with a vertex.

# Implementation

Custom types only need to implement `sequences(::T)`.
"""
function sequences end

sequences(genealogy, vs) = Iterators.map(v -> sequence(genealogy, v), vs)

sequences(genealogy, e::EdgeType) = sequences(genealogy, (src(e), dst(e)))

export mrca
"""
    mrca(genealogy[, vs])

Most recent common ancestor of a set of vertices. If omited, returns the mrca
of the whole genealogy.

See also [`tmrca`](@ref) for the time to the most recent common ancestor.

# Implementation

Custom type only need to implement `mrca(::T)`.
"""
function mrca end

# function mrca(genealogy, vs)
#     μ = mrca(genealogy)

#     possible_μ = CheapStack{VertexType}(nv(genealogy))

#     ## Mandatory to avoid dynamic dispatch.
#     for child ∈ children(genealogy, μ)
#         push!(possible_μ, child)
#     end

#     while !isempty(possible_μ)
#         v = pop!(possible_μ)
#         vs ⊆ descendants(genealogy, v) || continue

#         μ = v
#         empty!(possible_μ)

#         ## Mandatory to avoid dynamic dispatch.
#         for child ∈ children(genealogy, μ)
#             push!(possible_μ, child)
#         end
#     end

#     ## Type inference needs a little bit of help here for some reason...
#     something(μ)
# end

function mrca(genealogy, vs)
    μ = mrca(genealogy)

    @no_escape begin
        possible_μ = @alloc(eltype(genealogy), nv(genealogy))
        ptr = firstindex(possible_μ)

        for child ∈ children(genealogy, μ)
            possible_μ[ptr] = child
            ptr += 1
        end

        while ptr > firstindex(possible_μ)
            v = possible_μ[ptr - 1]
            ptr -= 1
            vs ⊆ descendants(genealogy, v) || continue

            μ = v
            ptr = firstindex(possible_μ)

            for child ∈ children(genealogy, μ)
                possible_μ[ptr] = child
                ptr += 1
            end
        end
    end

    μ
end

export distance
"""
    distance(genealogy, v1, v2)

Distance between two vertices.

# Implementation,

Only mandatory if copulas are to be fitted on the genealogy.
See [`loglikelihood`](@ref).
"""
function distance end

export prob
"""
    prob(genealogy; logscale = false)

Probability of the genealogy.

# Implementation,

Only mandatory if copulas are to be fitted on the genealogy.
See [`loglikelihood`](@ref).
"""
function prob end

export positions
"""
    positions(genealogy)

Positions of the markers.
"""
function positions end

#############
# Utilities #
#############

"""
    validate_positions(positions, nmarkers)

Validate a marker of positions. Tries to fix it if invalid. Return a vector
of uniformly spaced positions if it cannot be fixed.
"""
function validate_positions(positions, nmarkers)
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

    positions
end

"""
    idxtopos(genealogy, idx)

Return the position of the marker given its index.
"""
idxtopos(genealogy, idx) = getindex(positions(genealogy), idx)

"""
    postoidx(genealogy, pos)

Return the largest marker's index that is at a position lesser than the one
given.
"""
function postoidx(genealogy, pos)
    @inbounds for (k, p) ∈ enumerate(positions(genealogy))
        p >= pos && return k - 1
    end

    nmarkers(genealogy)
end

"""
    ancestral_mask(genealogy, x)

Mask non ancestral positions to 0.
"""
function ancestral_mask(genealogy, ω::Ω)
    lpos, rpos = endpoints(ω)

    lidx = 1
    while lpos > positions(genealogy)[lidx]
        lidx += 1
    end

    ridx = postoidx(genealogy, rpos)

    mask = falses(nmarkers(genealogy))
    mask[range(lidx, ridx)] .= true

    Sequence(mask)
end

ancestral_mask(genealogy, ωs::Set{Ω}) =
    mapreduce(ω -> ancestral_mask(genealogy, ω), |, ωs)

ancestral_mask(genealogy, x::EdgeType) =
    ancestral_mask(genealogy, ancestral_intervals(genealogy, x))

###################
# Pretty printing #
###################

function show(io::IO, ::MIME"text/plain", genealogy::AbstractGenealogy)
    println(io, describe(genealogy) * ":")
    println(io, nleaves(genealogy), " leaves, ",
            nmarkers(genealogy), " markers")
    print(io, "tMRCA: ", tmrca(genealogy))
end

function show(io::IO, genealogy::AbstractGenealogy)
    print(io, describe(genealogy, false) * "(")
    print(io, nleaves(genealogy))
    print(io, " leaves, ")
    print(io, nmarkers(genealogy))
    print(io, " markers)")
end

###########################
# AbstractGraph Interface #
###########################

## See https://juliagraphs.org/Graphs.jl/stable/ecosystem/interface/

for fun ∈ [:edges, :vertices, :ne, :nv]
    @eval function $fun(genealogy::AbstractGenealogy)
        $fun(graph(genealogy))
    end
end

for (fun, ret) ∈ Dict(:eltype => VertexType,
                      :edgetype => EdgeType,
                      :is_directed => :true)
    @eval begin
        @generated $fun(::Type{<:AbstractGenealogy}) = $ret
        @generated $fun(::AbstractGenealogy) = $ret
    end
end

let genealogy_edge = Expr(:(::), :e, EdgeType),
    genealogy_vertex = Expr(:(::), :v, VertexType)

    for (fun, a) ∈ Dict(:has_edge => genealogy_edge,
                        :has_vertex => genealogy_vertex,
                        :inneighbors => genealogy_vertex,
                        :outneighbors => genealogy_vertex)
        varname = first(a.args)
        @eval function $fun(genealogy::AbstractGenealogy, $a)
            $fun(graph(genealogy), $varname)
        end
    end
end

############
# Plotting #
############

"""
    maxdepth(genealogy, v)

Compute the depth of a vertex, that is the number of egdes between it and the
genealogy's mrca.
"""
function maxdepth(genealogy, v, depth = 0)
    _dads = dads(genealogy, v)
    isempty(_dads) && return depth

    mapreduce(d -> maxdepth(genealogy, d, depth + 1), max, _dads)
end

"""
    GenLayout(leaveslayer)

Layout function designed to plot genealogies. Leaves are all placed on the
same layer.
"""
GenLayout(leaveslayer, leaves) = function (genealogy_graph)
    xs, ys, _ = solve_positions(Zarate(),
                                genealogy_graph,
                                force_layer = Pair.(leaves, leaveslayer))

    ## Rotate by -π/2.
    Point.(zip(ys, -xs))
end

function graphplot(genealogy::AbstractGenealogy;
                   arrow_show = false,
                   edge_color = :gray,
                   edge_width = 3,
                   node_color = :black,
                   layout = nothing,
                   attributes...)
    if isnothing(layout)
        lastlayer = maximum(v -> maxdepth(genealogy, v), leaves(genealogy)) + 1
        layout = GenLayout(lastlayer, leaves(genealogy))
    end

    graphplot(graph(genealogy),
              layout = layout,
              ilabels = string.(range(1, nv(genealogy))),
              ilabels_color = :white,
              edge_color = edge_color,
              edge_width = edge_width,
              arrow_show = arrow_show,
              node_color = node_color,
              attributes...)
end

##################
# Common Methods #
##################

export isleaf, isroot, isivertex
"""
    isleaf(genealogy, v)
    isroot(genealogy, v)
    isivertex(genealogy, v)

True if `v` is a leaf/root in the given genealogy.
"""
function isleaf end,
function isroot end,
function isivertex end

isleaf(genealogy, v) = (iszero ∘ outdegree)(genealogy, v)

isroot(genealogy, v) = (iszero ∘ indegree)(genealogy, v)

isivertex(genealogy, v) = !isleaf(genealogy, v)

export nleaves, nivertices
"""
    nleaves(genealogy)
    nivertices(genealogy)

Number of leaves/internal vertices in a genealogy.

See also [`leaves`](@ref) and [`ivertices`](@ref) for an iterable over the
leaves/vertices.
"""
function nleaves end,
function nivertices end

nleaves(genealogy) = (length ∘ leaves)(genealogy)

nivertices(genealogy) = nv(genealogy) - nleaves(genealogy)

export sequence
"""
    sequence(genealogy, v)

Sequence of a genealogy associatex with a vertex.

See also [`sequences`](@ref) to get all the sequences of a genealogy.
"""
sequence(genealogy, v::VertexType) = getindex(sequences(genealogy), v)

export nmarkers
"""
    nmarkers(genealogy)

Number of markers in the sequences of a genealogy.
"""
nmarkers(genealogy) = (length ∘ first)(sequences(genealogy))

export branchlength
"""
    branchlength(genealogy[, e])

Total branch length of a genealogy. If an edge is specified, returns the length
of that edge.
"""
function branchlength end

branchlength(genealogy, e) = latitude(genealogy, src(e)) - latitude(genealogy, dst(e))

branchlength(genealogy) = mapreduce(e -> branchlength(genealogy, e), +, edges(genealogy))

export tmrca
function tmrca(genealogy)
    _mrca = mrca(genealogy)
    iszero(_mrca) && return zero(Float64)

    latitude(genealogy, _mrca)
end

tmrca(genealogy, vs) = latitude(genealogy, mrca(genealogy, vs))

export dads, children, siblings
"""
    dads(tree, v)
    children(tree, v)
    siblings(tree, v)

Parents/children/siblings of a vertex.
"""
function dads end,
function children end,
function siblings end

for (fun, graph_method) ∈ Dict(:dads => :inneighbors,
                               :children => :outneighbors)
    @eval function $fun(genealogy, v)
        $graph_method(genealogy, v)
    end
end

function siblings(genealogy, v)
    res = mapreduce(v -> children(genealogy, v), vcat, dads(genealogy, v))
    filter(!=(v), res)
end

export descendants
"""
    descendants(genealogy, v)

Transitive closure of the "children of" relation.
"""
function descendants(genealogy, v)
    descendants = CheapStack{VertexType}(nv(genealogy) - 1)
    _children = CheapStack{VertexType}(nv(genealogy) - 1)

    ## Mandatory to avoid dynamic dispatch.
    for child ∈ children(genealogy, v)
        push!(_children, child)
    end

    while !isempty(_children)
        v = pop!(_children)

        ## Mandatory to avoid dynamic dispatch.
        for child ∈ children(genealogy, v)
            push!(_children, child)
        end
        push!(descendants, v)
    end

    vec(descendants)
end

export nmutations
"""
    nmutations(genealogy[, e])

Number of mutation on a genealogy. If an edge is specified, return only the
number of mutations on that edge.
"""
function nmutations end

nmutations(genealogy, e) = count(xor(sequences(genealogy, e)...))

function nmutations(genealogy)
    mapreduce(e -> nmutations(genealogy, e), +, edges(genealogy),
              init = zero(Int))
end

## Edge functions.
for fun ∈ [:branchlength, :nmutations]
    @eval function $fun(genealogy, σ, δ)
        $fun(genealogy, Edge(σ, δ))
    end
end

export coalesce!
"""
    coalesce!(genealogy, v1, v2, lat)

Make two vertices coalesce at a specified latitude. Each marker in the
sequence associated with the newly created parent vertex is the conjunction (&)
of the corresponding marker in the children's sequences.
"""
function coalesce!(genealogy, v1, v2, lat)
    newseq = sequence(genealogy, v1) & sequence(genealogy, v2)

    add_vertex!(genealogy, newseq, lat) ||
        @error "Could not add vertex to genealogy"
    _dad = nv(genealogy)

    for child ∈ SA[v1, v2]
        e = Edge(_dad, child)
        add_edge!(genealogy, e)
    end

    _dad
end
