using ExploratoryDataAnalysis
using Documenter

DocMeta.setdocmeta!(ExploratoryDataAnalysis, :DocTestSetup, :(using ExploratoryDataAnalysis); recursive=true)

makedocs(;
    modules=[ExploratoryDataAnalysis],
    authors="Daymond Ling",
    repo="https://github.com/DaymondLing/ExploratoryDataAnalysis.jl/blob/{commit}{path}#{line}",
    sitename="ExploratoryDataAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://DaymondLing.github.io/ExploratoryDataAnalysis.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/DaymondLing/ExploratoryDataAnalysis.jl",
)
