using qAlgebra
using Documenter


makedocs(;
    modules=[qAlgebra],
    sitename="qAlgebra.jl Documentation",
    pages=[
        "Home" => "index.md",
        "API Reference" => [
            "Operator Spaces" => "api/statespace.md",
            "qExpressions" => "api/qExpressions.md",
            "Defaults" => "api/defaults.md"
        ]
    ],
)
