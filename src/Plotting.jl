module Plotting

export plot_ensemble_samples

"""
    plot_ensemble_samples(args...; kwargs...)

Placeholder emitted when the plotting extension is unavailable. Load `CairoMakie`
and `using QAlgebra.Plotting` to activate plotting support.
"""
function plot_ensemble_samples(args...; kwargs...)
    error("""
    QAlgebra.Plotting requires CairoMakie.
    Add CairoMakie to your environment and `using CairoMakie, QAlgebra.Plotting`
    to enable `plot_ensemble_samples`.
    """)
end

end # module
