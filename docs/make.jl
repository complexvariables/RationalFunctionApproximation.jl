using RationalFunctionApproximation
using Documenter, DocumenterVitepress
using DocumenterCitations
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))
include(joinpath(@__DIR__, "citations.jl"))

DocMeta.setdocmeta!(RationalFunctionApproximation, :DocTestSetup, :(using RationalFunctionApproximation); recursive=true)

makedocs(;
    modules=[RationalFunctionApproximation],
    authors="Toby Driscoll <driscoll@udel.edu>",
    repo=Remotes.GitHub("complexvariables", "RationalFunctionApproximation.jl"),
    sitename="RationalFunctionApproximation.jl",
    doctest=false,
    plugins=[bib],
    format=DocumenterVitepress.MarkdownVitepress(;
        repo = "https://github.com/complexvariables/RationalFunctionApproximation.jl",
        devbranch = "main",
    ),
    pages=[
        "Introduction" => "index.md",
        "Installation" => "install.md",
        "Algorithms" => "algorithms.md",
        "Domains" => "domains.md",
        "Discrete data" => "discrete.md",
        "Minimax" => "minimax.md",
        "Usage from Python" => "python.md",
        "Functions" => "functions.md",
        "References" => "references.md",
    ],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/complexvariables/RationalFunctionApproximation.jl",
    target = joinpath(@__DIR__, "build"),
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
