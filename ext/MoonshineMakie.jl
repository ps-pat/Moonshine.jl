module MoonshineMakie

using Moonshine
import Moonshine: plot_genealogy
using Moonshine: ancestral_mask

using Graphs

using GraphMakie: graphplot

const ∞ = Inf

function plot_genealogy(
    genealogy::AbstractGenealogy,
    ω;
    wild_color = :blue,
    derived_color = :red,
    arrow_show = false,
    edge_color = :gray,
    edge_width = 3,
    layout = Moonshine.plot_layout(genealogy),
    attributes...
)
    vlabels = string.(range(1, nv(genealogy)))

    ## Color of the vertices.
    mask = fill(ancestral_mask(genealogy, ω), nv(genealogy))
    node_color =
        ifelse.(any.(sequences(genealogy) .& mask), derived_color, wild_color)

    ## Hide non ancestral edges and vertices ##
    ewidth = Dict{Edge{Moonshine.VertexType}, Int}()
    for e ∈ edges(genealogy)
        ewidth[e] = isdisjoint(ancestral_intervals(genealogy, e), ω) ? 0 : edge_width
    end

    vsize = Dict{Moonshine.VertexType, Any}()
    default_vsize = 30
    for v ∈ leaves(genealogy)
        vsize[v] = default_vsize
    end

    for v ∈ ivertices(genealogy)
        if isdisjoint(ancestral_intervals(genealogy, v), ω)
            vsize[v] = 0
            vlabels[v] = ""
        else
            vsize[v] = default_vsize
        end
    end

    graphplot(
        graph(genealogy),
        layout = layout,
        node_size = vsize,
        ilabels = vlabels,
        ilabels_color = :white,
        edge_color = edge_color,
        edge_width = ewidth,
        arrow_show = arrow_show,
        node_color = node_color,
        attributes...
    )
end

plot_genealogy(genealogy::AbstractGenealogy; attributes...) =
    plot_genealogy(genealogy, Ω(0, ∞); attributes...)
end
