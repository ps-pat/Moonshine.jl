module MoonshineMakie

using Moonshine
import Moonshine: plot_genealogy

using Graphs

using DataStructures: DefaultDict

using GraphMakie: graphplot

const ∞ = Inf

function plot_genealogy(genealogy::AbstractGenealogy, ω;
                        wild_color = :blue,
                        derived_color = :red,
                        arrow_show = false,
                        edge_color = :gray,
                        edge_width = 3,
                        layout = Moonshine.plot_layout(genealogy),
                        attributes...)
    vlabels = string.(range(1, nv(genealogy)))

    ## Color of the vertices.
    mask = fill(ancestral_mask(genealogy, ω), nv(genealogy))
    node_color = ifelse.(any.(sequences(genealogy) .& mask),
                         derived_color, wild_color)

    ## Hide non ancestral edges and vertices ##
    ewidth = DefaultDict{Edge{Moonshine.VertexType}, Int}(edge_width)
    for e ∈ edges(genealogy)
        isdisjoint(ancestral_intervals(genealogy, e), ω) || continue
        ewidth[e] = 0
    end

    vsize = DefaultDict{Moonshine.VertexType, Any}(30)
    for v ∈ ivertices(genealogy)
        isdisjoint(ancestral_intervals(genealogy, v), ω) || continue
        vsize[v] = 0
        vlabels[v] = ""
    end

    graphplot(graph(genealogy),
              layout = layout,
              node_size = vsize,
              ilabels = vlabels,
              ilabels_color = :white,
              edge_color = edge_color,
              edge_width = ewidth,
              arrow_show = arrow_show,
              node_color = node_color,
              attributes...)
end

plot_genealogy(genealogy::AbstractGenealogy; attributes...) =
    plot_genealogy(genealogy, Ω(0, ∞); attributes...)
end
