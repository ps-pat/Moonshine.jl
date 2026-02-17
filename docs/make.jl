push!(LOAD_PATH,"../src/")

using Moonshine

using Documenter

using DocumenterInterLinks: InterLinks

links = InterLinks(
    "Graphs" => "https://juliagraphs.org/Graphs.jl/stable/objects.inv",
    "GraphMakie" => "https://graph.makie.org/stable/objects.inv",
    "GraphsFlows" => "https://juliagraphs.org/GraphsFlows.jl/dev/objects.inv",
    "Random" => "https://docs.julialang.org/en/v1/objects.inv",
    "Base" => "https://docs.julialang.org/en/v1/objects.inv",
    "IntervalSets" => "https://juliamath.github.io/IntervalSets.jl/stable/objects.inv")

writter = Documenter.HTMLWriter.HTML(
    assets = ["assets/custom.css"],
    canonical = "https://moonshine.patrickfournier.ca/stable",
    size_threshold_warn = 200 * 1024,
    size_threshold = nothing)

makedocs(sitename = "Moonshine.jl",
         format = Documenter.HTML(;collapselevel = 1),
         doctest = true,
         plugins = [links],
         repo = Remotes.GitHub("ps-pat", "Moonshine.jl"),
         pages = [
         "Home" => "index.md",
         "Design & conventions" => "design.md",
         "Guides" => [
             "QuickStart" => "quickstart.md",
             "Python & tskit" => "python.md"
         ],
         "Reference" => [
             "AbstractGenealogy" => "reference/AbstractGenealogy.md",
             "Tree & ARG" => "reference/Tree_Arg.md",
             "Ancestral Intervals" => "reference/AncestralIntervals.md",
             "Sequence" => "reference/Sequence.md",
             "Sample" => "reference/Sample.md",
             "CheapStack" => "reference/CheapStack.md",
             "ThreeTree" => "reference/ThreeTree.md",
             "Global Values" => "reference/globals.md",
             "Utilities" => "reference/utilities.md"
         ]
     ]
 )

if "publish" ∈ ARGS
    deploydocs(repo = "github.com/ps-pat/Moonshine.jl.git")
end
