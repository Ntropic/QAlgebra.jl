using qAlgebra
using Documenter


makedocs(;
    modules=[qAlgebra, qAlgebra.FFunctions],
    sitename="qAlgebra.jl Documentation",
    pages=[
        "Home" => "index.md",
        "API Reference" => [
            "Default Values" => "api/defaults.md",
            "Operator Spaces" => "api/qSpace.md",
            "qExpressions" => "api/qExpressions.md",
            "FFunctions" => "api/FFunctions.md", 
            "Defaults" => "api/defaults.md"
        ]
    ],
)
