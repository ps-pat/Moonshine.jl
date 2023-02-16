using Graphs:
    SimpleDiGraph,
    SimpleEdge,
    AbstractSimpleGraph,
    degree

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
    show

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

using Distributions:
    Bernoulli,
    Exponential,
    Uniform,
    DiscreteUniform,
    Categorical,
    truncated,
    logccdf,
    logpdf

const VertexType = Int
const EdgeType = SimpleEdge{VertexType}

######################
# ArgCore definition #
######################

struct ArgCore{T}
    graph::SimpleDiGraph{VertexType}
    latitudes::Vector{Float64}
    sequences::Vector{Sequence{T}}
    nleaves::Int
end

function ArgCore(leaves::AbstractArray{Sequence{T}}) where T
    n = length(leaves)
    latitudes = sizehint!(zeros(Float64, n - 1), 10n)
    sequences = sizehint!(leaves, 10n)

    ArgCore{T}(SimpleDiGraph(n), latitudes, sequences, n)
end

## TODO: Check performance impact of sizehint!.
function ArgCore{T}(rng::AbstractRNG,
                    nmin::Integer, minlength::Integer,
                    nmax::Integer = 0, maxlength::Integer = 0) where T
    n = iszero(nmax) ? nmin : rand(rng, nmin:nmax)
    nmarkers = iszero(maxlength) ? minlength : rand(rng, minlength:maxlength)

    ArgCore([Sequence{T}(rng, nmarkers) for _ ∈ 1:n])
end

function ArgCore{T}(nmin::Integer, minlength::Integer,
                    nmax::Integer = 0, maxlength::Integer = 0) where T
    ArgCore{T}(GLOBAL_RNG, nmin, minlength, nmax, maxlength)
end

#############################
# End of ArgCore definition #
#############################

##################
# Arg definition #
##################

mutable struct Arg{T} <: AbstractSimpleGraph{VertexType}
    core::ArgCore{T}
    nrecombinations::Int
    logprob::BigFloat
end

Arg(seqs::AbstractArray{Sequence{T}}) where T =
    Arg(ArgCore{T}(SimpleDiGraph(length(seqs)), Float64[], seqs), 0, zero(BigFloat))

function Arg{T}(rng::AbstractRNG,
                nmin::Integer, minlength::Integer,
                nmax::Integer = 0, maxlength::Integer = 0) where T
    Arg(ArgCore{T}(rng, nmin, minlength, nmax, maxlength), 0, zero(BigFloat))
end
function Arg(rng::AbstractRNG,
             nmin::Integer, minlength::Integer,
             nmax::Integer = 0, maxlength::Integer = 0)
    Arg{UInt}(rng, nmin, minlength, nmax, maxlength)
end

function Arg{T}(nmin::Integer, minlength::Integer,
                nmax::Integer = 0, maxlength::Integer = 0) where T
    Arg{T}(GLOBAL_RNG, nmin, minlength, nmax, maxlength)
end

function Arg(nmin::Integer, minlength::Integer,
             nmax::Integer = 0, maxlength::Integer = 0)
    Arg{UInt}(nmin, minlength, nmax, maxlength)
end

## AbstractGraph Interface.
## See https://juliagraphs.org/Graphs.jl/stable/ecosystem/interface/
for fun ∈ [:edges, :vertices, :ne, :nv, :add_vertex!]
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
                        :add_edge! => arg_edge,
                        :rem_edge! => arg_edge,
                        :has_vertex => arg_vertex,
                        :inneighbors => arg_vertex,
                        :outneighbors => arg_vertex)
        varname = first(a.args)
        @eval function $fun(arg::Arg, $a)
            $fun(arg.core.graph, $varname)
        end
    end
end

## Simple methods.
nleaves(arg) = arg.core.nleaves
leaves(arg) = Base.OneTo(arg.core.nleaves)
isleaf(arg, v) = v ∈ leaves(arg)

nmarkers(arg) = (length ∘ first)(arg.core.sequences)
sequences(arg) = arg.core.sequences
latitude(arg, v) =
    isleaf(v) ? zero(Float64) : arg.core.latitudes[v - nleaves(arg)]

mrca(arg) = argmax(arg.core.latitudes)
tmrca(arg) = maximum(arg.core.latitudes)

children(arg, v) = outneighbors(arg, v)
parents(arg, v) = inneighbors(arg, v)

#########################
# End of Arg definition #
#########################

##############################
# Plotting & pretty printing #
##############################

function argplot(arg, x;
                 arrow_show = true,
                 wild_color = :blue,
                 derived_color = :red,
                 attributes...)
    n = nv(arg)
    nedges = ne(arg)

    edgecolor = fill(:gray, nedges)

    nodecolor = map(sequences(arg)) do vertex
        vertex[x] ? derived_color : wild_color
    end

    node_labels = string.(1:n)

    node_sizes = fill(25.0, n)

    p = graphplot(arg,
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

#####################################
# End of plotting & pretty printing #
#####################################

####################
# ARG construction #
####################

function coalesce!(rng, arg, vertices, nlive)
    ## Select coalescing pair.
    arg.logprob += (log ∘ inv ∘ binomial)(length(vertices), 2)
    shuffle!(rng, vertices)

    ## Sample a latitude for the coalescence event.
    Δdist = Exponential(inv(nlive))
    Δ = rand(rng, Δdist)
    arg.logprob += logccdf(Δdist, Δ)

    ## Perform the coalescence.
    add_vertex!(arg) || @error "Could not add a vertex to ARG"
    parent = nv(arg)
    _children = Vector{eltype(arg)}(undef, 2)
    for child ∈ 1:2
        _children[child] = pop!(vertices)
        add_edge!(arg, parent, _children[child])
    end

    push!(vertices, parent)

    push!(sequences(arg),
          sequences(arg)[first(_children)] &
              sequences(arg)[last(_children)])


    @info("$(first(_children)) and $(last(_children)) \
           coalesced into $parent")
end

"""
    buildtree!(rng = GLOBAL_RNG, arg, idx = 1)

Build a marginal tree consistent with a specific index.
"""
function buildtree!(rng::AbstractRNG, arg::Arg, idx = 1)
    nv(arg) > nleaves(arg) && @warn "ARG contains non-leaf vertices"

    ## Group wild/derived sequences together to make the algorithm simpler.
    sort!(sequences(arg), lt = (x, y) -> x[idx] <= y[idx])

    nlive = nleaves(arg)
    nlive_wild = map(s -> s[idx], sequences(arg)) .|> iszero |> sum
    nlive_derived = nlive - nlive_wild

    wild = (collect ∘ range)(1, length = nlive_wild)
    derived = (collect ∘ range)(nlive_wild + 1, length = nlive_derived)

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

###########################
# End of ARG construction #
###########################
