import Pkg; Pkg.update();
using RationalFunctionApproximation
using Documenter, DocumenterVitepress

DocMeta.setdocmeta!(RationalFunctionApproximation, :DocTestSetup, :(using RationalFunctionApproximation); recursive=true)

makedocs(;
    modules=[RationalFunctionApproximation],
    authors="Toby Driscoll <driscoll@udel.edu> and contributors",
    #repo="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/{commit}{path}#{line}",
    repo=Remotes.GitHub("complexvariables", "RationalFunctionApproximation.jl"),
    sitename="RationalFunctionApproximation.jl",
    doctest=false,
    format=DocumenterVitepress.MarkdownVitepress(;
        repo = "https://github.com/complexvariables/RationalFunctionApproximation.jl",
    #format=Documenter.HTML(;
     #   prettyurls=get(ENV, "CI", "false") == "true",
      #  canonical="https://complexvariables.github.io/RationalFunctionApproximation.jl",
       # edit_link="main",
        #repolink="https://github.com/complexvariables/RationalFunctionApproximation.jl",
        #assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "Algorithms" => "algorithms.md",
        "Domains" => "domains.md",
        "Discrete data" => "discrete.md",
        "Minimax" => "minimax.md",
        "Usage from Python" => "python.md",
        "Functions" => "functions.md",
    ],
)

deploydocs(;
    repo="github.com/complexvariables/RationalFunctionApproximation.jl",
    target="build",
    devbranch="main",
    branch="gh-pages",
    push_preview=true
)
