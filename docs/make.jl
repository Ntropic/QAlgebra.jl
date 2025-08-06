using QAlgebra
using Documenter


makedocs(;
    modules=[QAlgebra, QAlgebra.CFunctions],
    sitename="QAlgebra.jl Documentation",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages=[
        "Home" => "index.md",
        "API Reference" => [
            "Default Values" => "api/defaults.md",
            "Operator Spaces" => "api/QSpace.md",
            "QExpressions" => "api/QExpressions.md",
            "CFunctions" => "api/CFunctions.md", 
            "Other" => "api/other.md"
        ]
    ],
)
