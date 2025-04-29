using qAlgebra
using Documenter

DocMeta.setdocmeta!(qAlgebra, :DocTestSetup, :(using qAlgebra); recursive=true)

makedocs(;
    modules=[qAlgebra],
    authors="Michael Schilling",
    sitename="qAlgebra.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API Reference" => [
            "Operator Spaces" => "api/qSpace.md",
            "qExpressions" => "api/qExpressions.md",
            "Defaults" => "api/defaults.md"
        ]
    ],
)
