using EpithelialDynamics1D
using Documenter

DocMeta.setdocmeta!(EpithelialDynamics1D, :DocTestSetup, :(using EpithelialDynamics1D); recursive=true)

makedocs(;
    modules=[EpithelialDynamics1D],
    authors="Daniel VandenHeuvel <danj.vandenheuvel@gmail.com>",
    repo="https://github.com/DanielVandH/EpithelialDynamics1D.jl/blob/{commit}{path}#{line}",
    sitename="EpithelialDynamics1D.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://DanielVandH.github.io/EpithelialDynamics1D.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/DanielVandH/EpithelialDynamics1D.jl",
    devbranch="main",
)
