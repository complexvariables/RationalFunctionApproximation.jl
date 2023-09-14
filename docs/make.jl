using RationalFunctionApproximation
using Documenter

DocMeta.setdocmeta!(RationalFunctionApproximation, :DocTestSetup, :(using RationalFunctionApproximation); recursive=true)

makedocs(;
    modules=[RationalFunctionApproximation],
    authors="Toby Driscoll <driscoll@udel.edu> and contributors",
    repo="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/{commit}{path}#{line}",
    sitename="RationalFunctionApproximation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://complexvariables.github.io/RationalFunctionApproximation.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/complexvariables/RationalFunctionApproximation.jl",
    devbranch="main",
)
