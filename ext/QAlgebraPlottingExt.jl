module QAlgebraPlottingExt

using CairoMakie
using LaTeXStrings: LaTeXString
import QAlgebra
using QAlgebra: QSpace, get_parameter_group, get_ensemble
import QAlgebra.Plotting
using QAlgebra.Sampler: pdf

const _DEFAULT_MARKER = :circle
const _DEFAULT_MARKERSIZE = 10
const _DEFAULT_MARKERCOLOR = :dodgerblue
const _DEFAULT_JOINT_RESOLUTION = 256

@inline _to_symbol(name::Symbol) = name
@inline _to_symbol(name::AbstractString) = Symbol(name)

_as_namedtuple(kwargs::NamedTuple) = kwargs
_as_namedtuple(kwargs::Base.Pairs) = (; kwargs...)
_as_namedtuple(::Nothing) = NamedTuple()
function _as_namedtuple(kwargs)
    error("Expected NamedTuple or keyword pairs, got $(typeof(kwargs)).")
end

function _resolve_aspect(aspect)
    aspect === nothing && return nothing
    if aspect isa Real
        return AxisAspect(aspect)
    elseif aspect isa Symbol
        if aspect in (:equal, :data)
            return DataAspect()
        else
            error("Unsupported symbolic aspect $(aspect). Use :equal, :data, or a numeric ratio.")
        end
    elseif aspect isa AxisAspect || aspect isa DataAspect
        return aspect
    else
        error("Unsupported aspect type $(typeof(aspect)).")
    end
end

@inline function _latex_label(qspace::QSpace, sym::Symbol)
    group_idx = get_parameter_group(qspace, sym)
    label = qspace.param_info.params_latex[group_idx]
    return LaTeXString(label)
end

@inline function _joint_line_defaults()
    return (; color=:black, linewidth=2)
end

@inline function _joint_contour_defaults()
    return (; colormap=:viridis, transparency=0.35, levels=12)
end

@inline function _joint_heatmap_defaults()
    return (; colormap=:viridis, transparency=0.35, interpolate=true)
end

@inline function _joint_contour_line_defaults()
    return (; colormap=:viridis, linewidth=2)
end

"""
    plot_ensemble_samples(qspace, ensemble_key; params=nothing, marker=:circle,
                           markersize=10, color=:dodgerblue, add_joint_probability=false,
                           joint_resolution=256, joint_mode=:heatmap, joint_kwargs=NamedTuple(), aspect=nothing,
                           axis_kwargs=NamedTuple(), kwargs...)

Visualise the sample positions registered on an ensemble inside `qspace`. The
`ensemble_key` can be either the ensemble's outer subspace symbol or one of its
inner keys. By default all sampled parameter groups (those backed by
`QDistribution`s) are plotted; supply `params` with a subset of parameter
symbols to restrict the axes. One- and two-dimensional plots are supported. If
`add_joint_probability` is true, the marginal/product distribution is overlayed as a line
(`ndims == 1`) or filled heatmap (`ndims == 2`). Additional keyword
arguments are forwarded to `Makie.scatter!` so you can customise markers (for example `(; color=:transparent, strokecolor=:black, strokewidth=1.5)` draws hollow circles).
Pass `joint_kwargs` to tweak the overlay styling (for example `(; transparency=0.35)` to control opacity or `(; colormap=:grays)` for greyscale). For 2D overlays a continuous heatmap is used by default; swap to filled contours with `joint_mode=:contourf`, or contour lines only with `joint_mode=:contour`. The legacy `alpha` keyword is translated automatically for Makie.

The returned `Figure` contains a single `Axis` with LaTeX-formatted labels that
match the requested parameters.
"""
function Plotting.plot_ensemble_samples(qspace::QSpace,
                                        ensemble_key::Union{Symbol,AbstractString};
                                        params::Union{Nothing,AbstractVector{<:Union{Symbol,AbstractString}}}=nothing,
                                        marker=_DEFAULT_MARKER,
                                        markersize::Real=_DEFAULT_MARKERSIZE,
                                        color=_DEFAULT_MARKERCOLOR,
                                        add_joint_probability::Union{Nothing,Bool}=nothing,
                                        joint_resolution::Int=_DEFAULT_JOINT_RESOLUTION,
                                        joint_mode::Symbol=:heatmap,
                                        joint_kwargs::Union{NamedTuple,Base.Pairs,Nothing}=NamedTuple(),
                                        aspect=nothing,
                                        axis_kwargs::Union{NamedTuple,Base.Pairs,Nothing}=NamedTuple(),
                                        kwargs...)
    sym = _to_symbol(ensemble_key)
    ensemble = get_ensemble(qspace, sym)
    sample = ensemble.sampler
    sample === nothing && error("Ensemble $(sym) does not have samples attached. Build samples before plotting.")

    group_lookup = Dict(sample.group_symbols[i] => i for i in eachindex(sample.group_symbols))
    selected_syms = params === nothing ? copy(sample.group_symbols) : _to_symbol.(params)
    isempty(selected_syms) && error("At least one parameter symbol is required for plotting.")

    selected_cols = Int[]
    for s in selected_syms
        idx = get(group_lookup, s, nothing)
        idx === nothing && error("Parameter symbol $(s) is not part of the sampling groups for ensemble $(sym).")
        push!(selected_cols, idx)
    end

    ndims = length(selected_cols)
    ndims > 2 && error("plot_ensemble_samples supports up to 2 parameters; got $(ndims).")

    coords = sample.samples[:, selected_cols]
    fig = Figure()
    aspect_setting = _resolve_aspect(aspect)
    axis_kw = _as_namedtuple(axis_kwargs)
    axis_kw = aspect_setting === nothing ? axis_kw : merge(axis_kw, (; aspect=aspect_setting))
    ax = Axis(fig[1, 1]; axis_kw...)

    scatter_kwargs = _as_namedtuple(kwargs)
    legacy_plot_joint = hasproperty(scatter_kwargs, :plot_joint) ? scatter_kwargs.plot_joint : nothing
    if legacy_plot_joint !== nothing
        if add_joint_probability !== nothing && add_joint_probability != legacy_plot_joint
            error("Both `add_joint_probability=$(add_joint_probability)` and legacy `plot_joint=$(legacy_plot_joint)` were provided. Use only one keyword or ensure they match.")
        end
        scatter_kwargs = (; (p for p in pairs(scatter_kwargs) if p.first != :plot_joint)...)
    end
    add_joint_flag = add_joint_probability === nothing ? (legacy_plot_joint === nothing ? false : legacy_plot_joint) : add_joint_probability

    joint_mode ∈ (:contourf, :heatmap, :contour) ||
        error("joint_mode must be one of :contourf, :heatmap, or :contour, got $(joint_mode).")
    overlay_mode = ndims == 2 ? joint_mode : :contourf

    joint_kw = _as_namedtuple(joint_kwargs)
    if :alpha in propertynames(joint_kw)
        alpha_val = joint_kw.alpha
        joint_kw = (; (p for p in pairs(joint_kw) if p.first != :alpha)..., transparency=alpha_val)
    end
    if overlay_mode === :heatmap && :levels in propertynames(joint_kw)
        joint_kw = (; (p for p in pairs(joint_kw) if p.first != :levels)...)
    end

    latex_labels = [_latex_label(qspace, s) for s in selected_syms]

    if ndims == 1
        xs = coords[:]
        dist = sample.distributions[selected_cols[1]]
        ax.xlabel = latex_labels[1]
        if add_joint_flag
            grid = range(dist.minimum, dist.maximum; length=joint_resolution)
            values = pdf.(Ref(dist), grid)
            line_kwargs = merge(_joint_line_defaults(), joint_kw)
            lines!(ax, grid, values; line_kwargs...)
            ax.ylabel = LaTeXString("p(" * string(latex_labels[1]) * ")")
        else
            ylims!(ax, -0.5, 0.5)
            ax.ylabel = LaTeXString(" ")
            ax.yticksvisible = false
            ax.yticklabelsvisible = false
            ax.ygridvisible = false
            ax.yminorgridvisible = false
            ax.yspinesvisible = false
            hlines!(ax, [0.0]; color=:gray70, linestyle=:dash, linewidth=1)
        end
        ys = add_joint_flag ? pdf.(Ref(dist), xs) : zeros(length(xs))
        scatter!(ax, xs, ys; marker=marker, markersize=markersize, color=color, scatter_kwargs...)
    else
        xs = coords[:, 1]
        ys = coords[:, 2]
        ax.xlabel = latex_labels[1]
        ax.ylabel = latex_labels[2]
        if add_joint_flag
            dist_x = sample.distributions[selected_cols[1]]
            dist_y = sample.distributions[selected_cols[2]]
            gx = range(dist_x.minimum, dist_x.maximum; length=joint_resolution)
            gy = range(dist_y.minimum, dist_y.maximum; length=joint_resolution)
            Z = [pdf(dist_x, x) * pdf(dist_y, y) for x in gx, y in gy]
            if overlay_mode === :heatmap
                heatmap_kwargs = merge(_joint_heatmap_defaults(), joint_kw)
                heatmap!(ax, gx, gy, Z; heatmap_kwargs...)
            elseif overlay_mode === :contourf
                contour_kwargs = merge(_joint_contour_defaults(), joint_kw)
                contourf!(ax, gx, gy, Z; contour_kwargs...)
            else
                contour_line_kwargs = merge(_joint_contour_line_defaults(), joint_kw)
                contour!(ax, gx, gy, Z; contour_line_kwargs...)
            end
        end
        scatter!(ax, xs, ys; marker=marker, markersize=markersize, color=color, scatter_kwargs...)
    end

    return fig
end

end # module
