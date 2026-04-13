push!(LOAD_PATH,"../src/")

using Moonshine

using Documenter

using DocumenterInterLinks: InterLinks

using Git: git

const draft = "draft" ∈ ARGS

## Determine which version we are building documentation for
current_version = try
    gitout = (read ∘ git)(["describe", "--tags", "--abbrev=0", "--exact-match"])
    version = mapreduce(Char, *, gitout[1:end - 1])
    match(r"v[0-9]+\.[0-9]+\.[0-9]+", version).match
catch
    "dev"
end

links = InterLinks(
    "Graphs" => "https://juliagraphs.org/Graphs.jl/stable/objects.inv",
    "GraphMakie" => "https://graph.makie.org/stable/objects.inv",
    "GraphsFlows" => "https://juliagraphs.org/GraphsFlows.jl/dev/objects.inv",
    "Random" => "https://docs.julialang.org/en/v1/objects.inv",
    "Base" => "https://docs.julialang.org/en/v1/objects.inv",
    "IntervalSets" => "https://juliamath.github.io/IntervalSets.jl/stable/objects.inv",
    "PythonCall" => "https://juliapy.github.io/PythonCall.jl/stable/objects.inv")

writter = Documenter.HTMLWriter.HTML(
    assets = ["assets/custom.css"],
    canonical = "https://patrickfournier.ca/software/documentation/moonshine/stable",
    edit_link = "master",
    highlights = ["python", "python-repl"],
    size_threshold_warn = 200 * 1024,
    size_threshold = nothing)

makedocs(build = "build/$current_version",
         draft = draft,
         sitename = "Moonshine.jl",
         format = Documenter.HTML(;collapselevel = 1),
         doctest = true,
         plugins = [links],
         repo = Remotes.GitHub("ps-pat", "Moonshine.jl"),
         pages = [
         "Home" => "index.md",
         "Design & conventions" => "design.md",
         "Guides" => [
             "QuickStart" => "quickstart.md",
             "Python & tskit" => "python.md",
             "Working with VCF files" => "vcf.md"
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
