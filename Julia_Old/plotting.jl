using CairoMakie
using LaTeXStrings
using Colors, ColorSchemes
using FixedPointNumbers
using Base.Threads
include("postprocessing.jl")
include("combinatorics.jl")
include("print_terms.jl")
include("fourier.jl")

macro usethreads(multithreaded, expr::Expr)
    ex = quote
        if $multithreaded
            Threads.@threads $expr
        else
            $expr
        end
    end
    esc(ex)
end
###########################################################################################################################
########### Plot Pulses ###################################################################################################
###########################################################################################################################

function plot_pulse(val_dict::Dict, pulse_param::Pulse_Param_Struct; n::Int=500, base_unit::String="Hz", factor=10^6, xlabel::String="", ylabel::String="", ops::Union{Symbol, Vector{Symbol}}=:all, fig::Union{Figure, Nothing}=nothing, ax_n::Int=1)
    # Plot a control pulse 
    # Args:
    #     val_dict: Dict: The dictionary containing \\beta and \\kappa values
    #     T: Float64: The total time of the pulse
    # Optional Args:
    #     n: Int: The number of points to plot (default 500)
    #     label: String: The label for the pulse (default "Control Pulse")
    #     base_unit: String: The base unit of the pulse (default "Hz")
    #     factor: Float64: The factor to multiply the pulse by (default 10^6)
    #     xlabel: String: The x-axis label (default "")
    #     ylabel: String: The y-axis label (default "")
    #     dolatex: Bool: Whether to use LaTeX for the labels (default true)
    #     kwargs: Any: Any other keyword arguments to pass to the plot function
    # Returns:
    #     fig: Figure: The figure object
    if ops == :all
        ops = [:real, :imag, :abs]
    end
    if isa(ops, Symbol)
        ops = [ops]
    end
    # check each entry of opt 
    for op in ops
        if !(op in [:abs, :real, :imag])
            error("ops should be :abs, :real, :imag or :all or a vector of these symbols")
        end
    end 
    pulse_fun = val_dict["\\beta"]
    kappa = val_dict["\\kappa"]
    T = pulse_param.T
    t_orig = collect(range(0, T, length=n))
    beta_vals = [pulse_fun(t) for t in t_orig]
    pulse_vals = sqrt.(kappa) * beta_vals
    ts, x_unit = transform_to_closest_unit(t_orig, "s", 1 / factor; dolatex=true)
    real_pulse_vals = real.(pulse_vals)
    imag_pulse_vals = imag.(pulse_vals)
    real_pulse_vals, imag_pulse_vals, y_unit = transform_to_closest_unit(real_pulse_vals, imag_pulse_vals, base_unit, factor; dolatex=true)
    if length(xlabel) == 0
        xlabel = L"$t$ [%$x_unit]"
    end
    if length(ylabel) == 0
        ylabel = L"$\sqrt{\kappa}\beta(t)$ [%$y_unit]"
    end
    if fig == nothing
        fig = Figure()
    end
    ax = Axis(fig[1, ax_n], xlabel=xlabel, ylabel=ylabel)
    counter = 0
    for op in ops
        counter += 1
        if op == :real
            lines!(ax, ts, real_pulse_vals, label=L"\text{Re}(\sqrt{\kappa}\beta(t))", color=Cycled(counter))
        elseif op == :imag
            lines!(ax, ts, imag_pulse_vals, label=L"\text{Im}(\sqrt{\kappa}\beta(t))", color=Cycled(counter))
        elseif op == :abs
            # area from - abs to abs
            alpha = 0.25
            abs_vals = abs.(real_pulse_vals+ 1im*imag_pulse_vals)
            color = Makie.wong_colors()[counter]
            color2 = (color, alpha)
            band!(ax, ts, -abs_vals, abs_vals, label=L"|\sqrt{\kappa}\beta(t)|", color=color2)
            #lines!(ax, ts, abs.(real_pulse_vals+ 1im*imag_pulse_vals), label=L"|\sqrt{\kappa}\beta(t)|", color=Cycled(counter))
            #lines!(ax, ts, -abs.(real_pulse_vals+ 1im*imag_pulse_vals), color = Cycled(counter))
        end
    end
    axislegend(ax, position=:rt)
    xlims!(ax, ts[1], ts[end])
    return fig
end
# function that automizes this plot 
function plot_fourier(val_dict::Dict, pulse_param::Pulse_Param_Struct; num_stds=4.0, n=2^10, reltol=1e-3, type::Symbol=:positive, ops::Union{Symbol, Vector{Symbol}}=:all, factor::Float64=1e6, base_unit::String="Hz", fig::Union{Figure, Nothing}=nothing, ax_n::Int=1)
    # ops can be :all, :abs, :real, :imag or a vector of these symbols 
    if ops == :all
        ops = [:real, :imag, :abs]
    end
    if isa(ops, Symbol)
        ops = [ops]
    end
    # check each entry of opt 
    for op in ops
        if !(op in [:abs, :real, :imag])
            error("ops should be :abs, :real, :imag or :all or a vector of these symbols")
        end
    end 
    fs, F_funs = fourier_in_interval(val_dict, pulse_param; num_stds=num_stds, num_points=n, reltol=reltol, type=type)
    fs, fs_unit = transform_to_closest_unit(fs, base_unit, factor, dolatex=true)
    real_F_funs, imag_F_funs, F_unit = transform_to_closest_unit(real.(F_funs), imag.(F_funs), base_unit, factor, dolatex=true)
    # generate figure 
    if fig == nothing
        fig = Figure()
    end
    ax = Axis(fig[1, ax_n], xlabel=L"f \text{ [%$fs_unit]}" , ylabel=L"\sqrt{\kappa}\beta(f) \text{ [%$F_unit]}")
    # add axe
    counter = 0
    for op in ops
        counter += 1
        if op == :real
            lines!(ax, fs, real_F_funs, label=L"\text{Re}(\sqrt{\kappa}\beta(f)) \text{ [%$F_unit]}", color=Cycled(counter))
        elseif op == :imag
            lines!(ax, fs, imag_F_funs, label=L"\text{Im}(\sqrt{\kappa}\beta(f)) \text{ [%$F_unit]}", color=Cycled(counter))
        elseif op == :abs
            alpha = 0.25
            abs_vals = abs.(real_F_funs + 1im*imag_F_funs)
            color = Makie.wong_colors()[counter]
            color2 = (color, alpha)
            band!(ax, fs, -abs_vals, abs_vals, label=L"|\sqrt{\kappa}\beta(f)| \text{ [%$F_unit]}", color=color2)
            #lines!(ax, fs, abs.(F_funs), label=L"|F(f)|", color=Cycled(counter))
            #lines!(ax, fs, -abs.(F_funs), color = Cycled(counter))
        end
    end
    # x axis tight 
    xlims!(ax, minimum(fs), maximum(fs))
    # add legend
    axislegend(ax)
    return fig
end
function plot_pulse_fourier(val_dict::Dict, pulse_param::Pulse_Param_Struct; num_stds=4.0, n::Int=1024, reltol=1e-3, base_unit::String="Hz", factor=10^6, type::Symbol=:positive, ops::Union{Symbol, Vector{Symbol}}=:all)
    #wider figure
    fig = Figure(size=(1200, 450))
    fig = plot_pulse(val_dict, pulse_param; n=n, base_unit=base_unit, factor=factor, fig=fig, ax_n=1)
    fig = plot_fourier(val_dict, pulse_param; num_stds=num_stds, n=n, reltol=reltol, type=type, ops=ops, factor=factor, base_unit=base_unit, fig=fig, ax_n=2)
    return fig
end

###########################################################################################################################
########### Plot Weighted Samples and Sample Locations in 2D ##############################################################
###########################################################################################################################

function make_scatter_weight_legend(f, markersizes, markerstr; strokewidth=1.0)
    #colorstr = ["-", "+"]
    colorstr = [L"\text{negative}", L"\text{positive}"]

    colors = [Cycled(2), Cycled(1)]

    group_size = [MarkerElement(marker=:circle, color=:lightgray, markersize=ms, strokewidth=strokewidth) for ms in markersizes]

    group_color = [MarkerElement(marker=:circle, color=color, markersize=markersizes[end], strokewidth=strokewidth) for color in colors]

    legend = Legend(f,
        [group_size, group_color],
        [markerstr, colorstr],
        [L"Weight: $\frac{w_i}{w_\text{max}}$", L"Sign: $\text{sgn}(w_i)$"], tellheight=true)

    legend.orientation = :horizontal
    legend.tellheight = true
    legend.tellwidth = false
    return legend
end
function make_scatter_weight_legend_small(f, markersizes, markerstr; strokewidth=1.0)

    group_size = [MarkerElement(marker=:circle, color=:lightgray, markersize=ms, strokewidth=strokewidth) for ms in markersizes]

    legend = Legend(f,
        [group_size],
        [markerstr],
        [L"Weight: $\frac{w_i}{w_\text{max}}$"], tellheight=true)

    legend.orientation = :horizontal
    legend.tellheight = true
    legend.tellwidth = false
    return legend
end
function plot_weighted_locations_2D(samples::SamplesAndWeights, dir_name::String; do_weighted::Bool=true, ms=20, ms2=3, base_unit::String="Hz", factor=10^6, xlabel::String="", ylabel::String="", strokewidth::Real=1.0, size=(600, 525), annotations::Bool=false, kwargs...)
    locations = samples.locations
    x_loc = [l[1] for l in locations]
    y_loc = [l[2] for l in locations]
    range_vec = samples.range_vec
    xlims = range_vec[1]
    ylims = range_vec[2]
    # rescale locations 
    xloc, xlims, x_unit = transform_to_closest_unit(x_loc, xlims, base_unit, factor)
    yloc, ylims, y_unit = transform_to_closest_unit(y_loc, ylims, base_unit, factor)

    if isempty(xlabel)
        xlabel = latexstring(samples.var_strs[1] * "\\text{ [", x_unit, "]}")
    end
    if isempty(ylabel)
        ylabel = latexstring(samples.var_strs[2] * "\\text{ [", y_unit, "]}")
    end

    curr_weights = samples.weights_dict[dir_name]
    max_weight = maximum(abs.(curr_weights))
    all_sizes = [sqrt(abs(w) / max_weight) * ms for w in curr_weights]
    sizes = [sqrt(abs(s)) * ms for s in [0.01, 0.1, 1.0]]
    sizestr = [L"0.01", L"0.1", L"1.0"]
    fig = Figure(size=size)
    if minimum(curr_weights) < 0.0
        fig[1, 1] = make_scatter_weight_legend(fig, sizes, sizestr; strokewidth=strokewidth)
    else
        fig[1, 1] = make_scatter_weight_legend_small(fig, sizes, sizestr; strokewidth=strokewidth)
    end
    ax = Axis(fig[2, 1], xlabel=xlabel, ylabel=ylabel; kwargs...)
    # find indexes that sort curr_weights in ascending order 
    for i in length(xloc):-1:1
        color = curr_weights[i] < 0 ? Cycled(2) : Cycled(1)
        if do_weighted
            scatter!(ax, [xloc[i]], [yloc[i]], markersize=all_sizes[i], color=color, strokewidth=strokewidth)
        else
            scatter!(ax, [xloc[i]], [yloc[i]], markersize=ms2, color=color, strokewidth=strokewidth)
        end
        if annotations
            text!(ax, string(i), position=(xloc[i], yloc[i]), fontsize=14, color=:black, align=(:bottom, :left))
        end
    end
    xlims!(ax, xlims)
    ylims!(ax, ylims)
    return fig
end

###########################################################################################################################
########### Plot Time Evolutions of Operators #############################################################################
###########################################################################################################################

function ribbon!(ax, t, mean_vals, std_vals; label="", color=:black, color_ind::Int=-1, alpha=0.3, do_ribbon::Bool=true, do_clamped::Bool=false, clamps::Tuple{Float64,Float64}=(-1.0, 1.0))
    # color_ind is an alternative to color, if > 0 it will use the corresponding color from the wong palette
    # alpha is only used if color_ind > 0
    if color_ind > 0
        color = Makie.wong_colors()[color_ind]
        color2 = (color, alpha)
    else
        color2 = color
    end
    # Plot mean line
    lines!(ax, t, mean_vals, label=label, color=color)
    # Add ribbon
    if do_ribbon
        min_vals = mean_vals .- std_vals
        max_vals = mean_vals .+ std_vals
        if do_clamped
            min_vals = max.(min_vals, clamps[1])
            max_vals = min.(max_vals, clamps[2])
        end
        band!(ax, t, min_vals, max_vals, color=color2)
    end
    xlims!(ax, t[1], t[end])
end
function clamp_opstr(plot_str::String)
    if !(occursin("+", plot_str) || occursin("-", plot_str))
        return true
    else
        return false
    end
end
function plot_single_expectation_value_over_t(fig::Figure, ax::Axis, i::Int, curr_ind::Int, op_str::String, op_prefix::String, op_suffix::String, solution::DifferentialSolution, op2::Function; do_cumulants::Bool=false, base_unit::String="s", factor=10^6, plot_trajectories::Bool=false, plot_sampled_stat::Bool=false, alpha::Float64=0.5, do_ribbon::Bool=true, do_mean::Bool=true, clamps::Bool=true, my_clamps::Tuple{Float64,Float64}=(-1.0, 1.0), t=-1, threaded::Bool=true)
    samples = solution.samples
    do_clamped = clamp_opstr(op_str)
    op, op_string = get_op_and_string(op_str, opstr)
    curr_str = latexstring(op_prefix * op_string * op_suffix)
    color = Cycled(i)
    curr_values, t_vals, curr_spin_orders = values_at_t(op_str, solution; do_cumulants=do_cumulants, t=t, threaded=threaded)
    t_vals, t_unit = transform_to_closest_unit(t_vals, base_unit, 1 / factor, dolatex=true)
    if plot_trajectories
        normalized_weights = get_weights(samples, curr_spin_orders, normalize=true)
        n_trajectories = size(curr_values, 1)

        # plot <=max_how_many_trajectories transjectories
        max_how_many_trajectories = 500
        max_weight_strength = 0.02
        how_many_trajectories = 0
        for j in 1:n_trajectories  
            if normalized_weights[j] > max_weight_strength
                how_many_trajectories += 1
            end
        end
        plot_every = max(1, Int(round(how_many_trajectories / max_how_many_trajectories)))
        linecolor = :gray
        every = 0
        # sort indexes by weight 
        indexes = sortperm(normalized_weights, rev=true)
        for j in indexes
            if normalized_weights[j] > max_weight_strength
                every += 1
                if every == 1
                    if j == 1
                        lines!(ax, t_vals, op2.(op.(curr_values[j, :])), linewidth=normalized_weights[j], label=nothing, color=linecolor)#, label=label_str, args...; kwargs...)
                    else
                        lines!(ax, t_vals, op2.(op.(curr_values[j, :])), linewidth=normalized_weights[j], label=nothing, color=linecolor)#, label=label_str, args...; kwargs...)
                    end
                end
                if every == plot_every
                    every = 0
                end
            end
        end
    end
    if do_mean
        curr_ind += 1
        # try unsampled
        if haskey(samples.order_to_exp_std_dict, curr_spin_orders)
            has_key = true
        end
        if has_key
            curr_exp, curr_std = get_exp_and_std_to_vals_in_t(curr_values, curr_spin_orders, samples)
        else
            # Do sampled instead
            normalized_weights = get_weights(samples, curr_spin_orders, normalize=true)
            curr_exp, curr_std = weighted_mean_std(curr_values, normalized_weights)
        end
        if length(curr_exp) == 0
            error("No values for ", op_str, " (likely cause: Cumulants or lower order cumulants haven't been saved, activate in solver)")
        end
        ribbon!(ax, t_vals, op2.(op.(curr_exp)), op.(curr_std), label=curr_str, color=color, color_ind=curr_ind, alpha=alpha, do_ribbon=do_ribbon, do_clamped=do_clamped, clamps=my_clamps) # Ribbon plot with error bars
    end
    if plot_sampled_stat
        weights = get_weights(samples, curr_spin_orders, normalize=true)
        mean_vals, std_vals = weighted_mean_std(curr_values, weights)
        curr_ind += 1
        ribbon!(ax, t_vals, op2.(op.(mean_vals)), op.(std_vals), label=L"Sampled %$curr_str", color=color, color_ind=curr_ind, alpha=alpha, do_ribbon=do_ribbon, do_clamped=do_clamped, clamps=my_clamps) # Ribbon plot with error bars
    end
    return fig, ax
end
function plot_expectation_values_over_t(plot_which::Vector{String}, solution::DifferentialSolution; plottype::String="x", opstr::String="real", do_cumulants::Bool=false, plot_trajectories::Bool=false, plot_sampled_stat::Bool=false, base_unit::String="s", factor=10^6, xlabel::String="", ylabel::String="", alpha::Float64=0.5, do_ribbon::Bool=true, do_mean::Bool=true, clamps::Bool=true, t=-1, threaded::Bool=true)
    # find indexes of plot_which in op_index_vec
    # solution is a DifferentialSolution object 
    # plot_which is a vector of strings of the operators to plot, i.e. ["z", "zz"]
    # system is a SystemAndDicts object
    # samples is a SamplesAndWeights object
    # opstr is a string that specifies if the operators are (real, imag, abs, abs2)
    # args and kwargs are arguments to the plot function
    # do_cumulants and do_lower_order_cumulants are booleans to decide if instead of expectation values cumulants or "full cumulants" of lower order should be plotted
    # plottype specifies if <op>, 1-<op>, <op>-1, <op>+1 and 1+<op> should be plotted via "x", "1-x", "x-1", "x+1" or "1+x"
    # alpha is the transparency of the ribbon plot (0.5)
    system = solution.system
    t_vals, do_interpolation = which_times(solution, t=t)
    if !do_interpolation
        t_ind = t_vals[1]:t_vals[end]
        t_vals = solution.t[t_ind]
    elseif isa(t_vals, AbstractRange)
        t_vals = collect(t_vals)
    end
    count_xyz_sum = input_string -> sum(count.(c -> c == 'x' || c == 'y' || c == 'z', input_string))
    spin_orders = maximum([count_xyz_sum(p) for p in plot_which])
    op2strlist::Vector{String} = ["x", "1-x", "x-1", "x+1", "1+x"]
    op2ind::Int = findfirst(x -> x == plottype, op2strlist)
    if op2ind == nothing
        error("plottype must be one of ", op2strlist)
    end

    op2::Function = [x -> x, x -> 1 - x, x -> x - 1, x -> x + 1, x -> 1 + x][op2ind]
    my_clamps = (op2(-1.0), op2(1.0))
    op_prefix::String = ["", "1-", "", "", "1+"][op2ind]
    op_suffix::String = ["", "", "-1", "+1", ""][op2ind]

    new_t_vals, t_unit = transform_to_closest_unit(t_vals, base_unit, 1 / factor, dolatex=true)
    if length(xlabel) == 0
        xlabel = L"$t$ [%$t_unit]"
    end

    if length(ylabel) == 0
        sig_str = get_op_and_string(plot_which[1], opstr)[2]
        if length(plot_which) > 1
            sig_str = get_op_and_string("O", opstr)[2]
        end
        do_lower_order_cumulants::Bool = false
        if haskey(system.operator_strings_dict, op_str) && do_cumulants
            do_lower_order_cumulants = true
        end
        if do_lower_order_cumulants
            ylabel = latexstring(op_prefix * sig_str * "_c" * op_suffix)
        elseif do_cumulants
            signum = "-"
            if occursin("-", op_prefix)
                signum = "+"
            end
            ylabel = latexstring(op_prefix * sig_str * signum * sig_str * "_c" * op_suffix)
        else
            ylabel = latexstring(op_prefix * sig_str * op_suffix)
        end
    end
    n_plot = length(plot_which)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=xlabel, ylabel=ylabel)
    curr_ind = 0
    for i in 1:n_plot
        op_str = plot_which[i]
        fig, ax = plot_single_expectation_value_over_t(fig, ax, i, curr_ind, op_str, op_prefix, op_suffix, solution, op2; base_unit=base_unit, factor=factor, do_cumulants=do_cumulants, plot_trajectories=plot_trajectories, plot_sampled_stat=plot_sampled_stat, alpha=alpha, do_ribbon=do_ribbon, do_mean=do_mean, clamps=clamps, my_clamps=my_clamps, t=t, threaded=threaded)
    end
    # add clamps?
    if clamps
        # only if no +- in at least one operator 
        do_clamps = false
        for plot_str in plot_which
            if !(occursin("+", plot_str) || occursin("-", plot_str))
                do_clamps = true
                break
            end
        end
        # -1 and +1 horizontal dotted lines, but shifted via op2 
        if do_clamps
            #for yval in my_clamps
            #    lines!(ax, t[[1, length(t)]], yval .* ones(length(2)), linestyle=:dot, color=:black, linewidth=1)
            #end
            # change y axis limits to my_clamps 
            ylims!(ax, my_clamps)
        end
    end
    # Add legend 
    if do_mean || plot_sampled_stat
        axislegend(ax, merge=true, unique=true)
    end
    return fig
end
function plot_expectation_values_over_t(plot_which::String, solution::DifferentialSolution; plottype::String="x", opstr::String="real", do_cumulants::Bool=false, plot_trajectories::Bool=false, plot_sampled_stat::Bool=false, base_unit::String="s", factor=10^6, xlabel::String="", ylabel::String="", do_ribbon::Bool=true, do_mean::Bool=true, clamps::Bool=true, t=-1, threaded::Bool=true)
    return plot_expectation_values_over_t([plot_which], solution; plottype=plottype, opstr=opstr, do_cumulants=do_cumulants, plot_trajectories=plot_trajectories, plot_sampled_stat=plot_sampled_stat, base_unit=base_unit, factor=factor, xlabel=xlabel, ylabel=ylabel, do_ribbon=do_ribbon, do_mean=do_mean, clamps=clamps, t=t, threaded=threaded)
end

function plot_cumulants_vs_cumulant_expansion(op_str::String, solution::DifferentialSolution; plottype::String="x", base_unit::String="s", factor=10^6, opstr::String="real", alpha=0.5, difference::Bool=false, trajectories::Bool=false)
    system = solution.system
    samples = solution.samples
    op2strlist::Vector{String} = ["x", "1-x", "x-1", "x+1", "1+x"]
    op2ind::Int = findfirst(x -> x == plottype, op2strlist)
    if op2ind == nothing
        error("plottype must be one of ", op2strlist)
    end
    op2::Function = [x -> x, x -> 1 - x, x -> x - 1, x -> x + 1, x -> 1 + x][op2ind]
    op_prefix::String = ["", "1-", "", "", "1+"][op2ind]
    op_suffix::String = ["", "", "-1", "+1", ""][op2ind]

    # check if the op_Str s in cumulants dictionary
    if haskey(system.cumulant_terms_dict, op_str)
        error("The operator string is already in the cumulants dictionary, we do not have the evolved operator, only it's cumulant expansion.")
    end
    curr_vals, t_pos, curr_spin_orders = values_at_t(op_str, solution, do_cumulants=true, t=t)[1:3]
    cumulant_vals = curr_vals[:, 1] # make vector
    real_vals = get_all_values_to_op_str_in_t(op_str, solution, system)[1][:,1] # make vector

    t, t_unit = transform_to_closest_unit(t_pos, base_unit, 1 / factor, dolatex=true)
    xlabel = L"$t$ [%$t_unit]"
    op, op_string = get_op_and_string(op_str, opstr)
    ylabel = L"%$op_prefix %$op_string %$op_suffix"

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=xlabel, ylabel=ylabel)
    if !difference
        # get the cumulant values and operator statistics
        exp_op_cum, std_op_cum = get_exp_and_std_to_vals_in_t(cumulant_vals, curr_spin_orders, samples)
        exp_op, std_op = get_exp_and_std_to_vals_in_t(real_vals, curr_spin_orders, samples)
        # plot each of them 
        ribbon!(ax, t, op2.(op.(exp_op)), op.(std_op), label=latexstring(op_prefix * op_string * op_suffix), color=Cycled(1), color_ind=1, alpha=alpha)
        ribbon!(ax, t, op2.(op.(exp_op_cum)), op.(std_op_cum), label=latexstring(op_prefix * op_string * "_c" * op_suffix), color=Cycled(2), color_ind=2, alpha=alpha)
    else
        diff_vals = cumulant_vals .- real_vals
        exp_diff, std_diff = get_exp_and_std_to_vals_in_t(diff_vals, curr_spin_orders, samples)
        ribbon!(ax, t, op2.(op.(exp_diff)), op.(std_diff), label=latexstring(op_prefix * op_string * "_c" * op_suffix * " - " * op_prefix * op_string * op_suffix), color=Cycled(1), color_ind=1, alpha=alpha)
        if trajectories
            normalized_weights = get_weights(samples, curr_spin_orders, normalize=true)
            n_trajectories = size(cumulant_vals, 1)
            for j in 1:n_trajectories
                if normalized_weights[j] > 0.02
                    lines!(ax, t, op2.(op.(cumulant_vals[j, :] .- real_vals[j, :])), linewidth=normalized_weights[j], label=nothing, color=:gray)
                end
            end
        end
    end
    axislegend(ax, merge=true, unique=true)
    return fig
end

###########################################################################################################################
########### Plot Basis Functions in 2D ####################################################################################
###########################################################################################################################

function plot_basis_function_coeffs_2D(op_str::String, solution::DifferentialSolution; do_title::Bool=false, do_annotate::Bool=true, plottype::String="coeffs", do_cumulants=false, fontsize::Int=8, t=-2)
    # Plot the coefficients of the basis functions asa function of the order of the basis functions in 2D at the last timestep
    # Required arguments:
    #  op_str is the operator string
    #  solution is a DifferentialSolution object
    # Optional arguments:
    #  do_title is a boolean that specifies if the plot should have a title
    #  do_annotate is a boolean that specifies if the plot should be annotated
    #  plottype is a string that specifies what to plot, "coeffs", "integrals", "coeffs_integrals", 
    #           so that the coefficients, the integrals (of the weighted basis functions) or the product of both are plotted
    #  do_cumulants is a boolean that specifies if the cumulants of the operator should be plotted
    #  do_lower_order_cumulants is a boolean that specifies if the cumulants of lower order should be plotted
    #  fontsize is the fontsize of the annotations (default 8)
    if !isa(t, Int) && !isa(t, Float64)
        error("t must be an integer or float")
    end
    opstr = "abs"
    op, op_string = get_op_and_string(op_str, opstr)
    samples::SamplesAndWeights = solution.samples

    order_comb::Vector{Vector{Int}} = samples.order_comb
    max_orders_x::Int = maximum(order_comb[i][1] for i in 1:length(order_comb))
    max_orders_y::Int = maximum(order_comb[i][2] for i in 1:length(order_comb))

    curr_values = values_at_t(op_str, solution, do_cumulants=do_cumulants, t=t)[1][:,1] # make vector
    vals2coeffs::Matrix{Float64} = samples.vals2coeffs
    coeffs::Vector{ComplexF64} = vals2coeffs * curr_values
    integrals_vector = samples.integrals_vector

    coeffs_by_order::Matrix{Float64} = zeros(max_orders_x + 1, max_orders_y + 1)
    do_c = false
    do_int = false
    for i in 1:length(coeffs)
        ci = order_comb[i][1] + 1
        cj = order_comb[i][2] + 1
        curr_val::ComplexF64 = 1.0
        if contains(plottype, "coeffs")
            curr_val *= coeffs[i]
            do_c = true
        end
        if contains(plottype, "integrals")
            curr_val *= integrals_vector[i]
            do_int = true
        end
        coeffs_by_order[ci, cj] = op.(curr_val)
    end
    text_types = ["Coefficients", "Integrals", "Coefficients and Integrals"]
    eq_types = ["c(p)", "I(p)", "c(p)I(p)"]
    eq_str = eq_types[do_c+2*do_int]
    pre_text = text_types[do_c+2*do_int]
    parameters::Vector{String} = samples.var_strs
    xlabel = L"Harmonic of $%$(parameters[1])$"  # Modes
    ylabel = L"Harmonic of $%$(parameters[2])$"  # Modes
    # plot coeffs in 2d as heatmap
    colorbarstring = latexstring("\\log_{10}(|" * eq_str * "|)")
    fig = Figure()
    if do_title 
        titlestring = L"%$pre_text of $%$op_string$"
        ax = Axis(fig[1, 1], xlabel=xlabel, ylabel=ylabel, title=titlestring)
    else
        ax = Axis(fig[1, 1], xlabel=xlabel, ylabel=ylabel)
    end
    corrected_orders = log10.(coeffs_by_order)'
    # replace values of coeffs_by_order that are 0 with NaN in 
    corrected_orders[abs.(coeffs_by_order').<10^-15] .= NaN
    hm = heatmap!(ax, collect(0:max_orders_x), collect(0:max_orders_y), corrected_orders')
    Colorbar(fig[1, 2], hm, label=colorbarstring)
    if do_annotate
        for i in 1:length(coeffs)
            ci = order_comb[i][1] + 1
            cj = order_comb[i][2] + 1
            if (coeffs_by_order[ci, cj]) > 10^-15
                f = Printf.Format("%0.1f")
                text!(ax, Printf.format(f, log10.(coeffs_by_order[ci, cj])), position=(ci - 1.0, cj - 1.0), color=:white, fontsize=fontsize * 2, align=(:center, :center))
            end
        end
    end
    return fig
end
## Test 
#op_str = "z"
#fig = plot_basis_function_coeffs_2D(op_str, solution, do_annotate=true, plottype="coeffs", do_cumulants=false, fontsize=8)

###########################################################################################################################
########### Plot Function Distribution of Samples #########################################################################
###########################################################################################################################

function rgb2weightedrgb(rgb::Matrix{RGB{Float64}}, weights::Matrix{Float64})
    #weighted_values = Matrix{RGB{N0f8}}(undef, size(rgb))
    weighted_values = Matrix{RGB{Float64}}(undef, size(rgb))
    for i in 1:size(rgb, 1)
        for j in 1:size(rgb, 2)
            r = clamp(rgb[i, j].r * weights[i, j], 0, 1)
            g = clamp(rgb[i, j].g * weights[i, j], 0, 1)
            b = clamp(rgb[i, j].b * weights[i, j], 0, 1)
            weighted_values[i, j] = RGB{Float64}(r, g, b)
            #weighted_values[i, j] = RGB{N0f8}(r, g, b)
        end
    end
    return weighted_values
end
function weighted_heatmap(values::Matrix{Float64}, weights::Matrix{Float64}; minmax::Vector{Float64}=[-1.0, -1.0], color_map=:viridis, weighted::Bool=true, rel_col_range::Float64=1.0)
    # Ensure values and weights matrices are the same size
    if size(values) != size(weights)
        throw(ArgumentError("The dimensions of values and weights matrices must match."))
    end
    if minmax == [-1.0, -1.0]
        minmax = [minimum(values), maximum(values)]
    end
    # Normalize values to [0, 1/rel_col_range] for color mapping
    mean_values = (maximum(values) + minimum(values)) / 2
    delta_values = (maximum(values) - minimum(values)) / 2
    d_values = (values .- mean_values) / rel_col_range
    new_values = d_values .+ mean_values
    norm_values = (new_values .- minmax[1]) ./ ((minmax[2] - minmax[1]))

    # Map normalized values to colors
    norm_values_int = ceil.(Int, norm_values * 256)
    # clamp to [1, 256]
    norm_values_int = max.(1, min.(256, norm_values_int))

    # iterate over the color map and get the color values
    value_array = to_colormap(color_map)
    # from rgba to rgb 
    value_array = [RGB{Float64}(value_array[i].r, value_array[i].g, value_array[i].b) for i in 1:256]
    color_mapped_values = value_array[norm_values_int]
    #color_mapped_values = get_color_palette(color_map, 256)[norm_values_int]

    ## Apply weights to the alpha channel
    min_max_weights = [minimum(weights), maximum(weights)]
    if weighted
        normalized_weights = (weights .- min_max_weights[1]) ./ (min_max_weights[2] - min_max_weights[1])
    else
        normalized_weights = ones(size(weights))
    end
    colors_with_alpha = rgb2weightedrgb(color_mapped_values, normalized_weights)
    return colors_with_alpha, min_max_weights
end
function plot_operator_as_function_2D(op_str::String, solution::DifferentialSolution; t=-2, do_title::Bool=false, weighted::Bool=true, opstr::String="real", res::Int=256, cmap=:viridis, rel_range::Float64=1.0, rel_col_range::Float64=1.0, do_cumulants::Bool=false, base_unit::String="Hz", base_unit_t::String=raw"\mu\text{s}", factor=10^6, xlabel::String="", ylabel::String="", plottype="1+x", figsize=(725, 450), do_locations::Bool=false)
    # Plot the function values of the operator op_str as a function of the integrated parameters using color for amplitude and alpha channel for relative probability of spins at the parameter location
    # Required arguments:
    #  op_str: string, the operator to plot
    #  solution: DifferentialSolution, the solution of the system
    #  system: SystemAndDicts, the system of differential equations
    #  samples: SamplesAndWeights, the samples and weights for the parameters
    # Optional arguments:
    #  t: int, the time index to plot (default -2 = last time index)
    #  do_title: bool, whether to add a title to the plot
    #  weighted::Bool, whether to use the weights in the solution
    #  opstr: string, the operation to perform on the operator, "abs", "real" or "imag"
    #  res: int, the resolution of the heatmap
    #  rel_range: float, the relative range of the heatmap, so that range of plot is range_vec*rel_range using the range_vec from the samples
    #  rel_col_range: float, the relative range of the color map, so that the color map is from min(abs(values)) to rel_col_range*max(abs(values))
    #  do_cumulants: bool, whether to use cumulants in the solution
    #  do_lower_order_cumulants: bool, whether to use lower order cumulants in the solution
    # plottype specifies if <op>, 1-<op>, <op>-1, <op>+1 and 1+<op> should be plotted via "x", "1-x", "x-1", "x+1" or "1+x"
    # figsize: tuple, the size of the figure
    if !isa(t, Int) && !isa(t, Float64)
        error("t must be an integer or float")
    end
    samples = solution.samples
    op2strlist::Vector{String} = ["x", "1-x", "x-1", "x+1", "1+x"]
    op2ind::Int = findfirst(x -> x == plottype, op2strlist)
    if op2ind == nothing
        error("plottype must be one of ", op2strlist)
    end
    op2::Function = [x -> x, x -> 1 - x, x -> x - 1, x -> x + 1, x -> 1 + x][op2ind]
    op_prefix::String = ["", "1-", "", "", "1+"][op2ind]
    op_suffix::String = ["", "", "-1", "+1", ""][op2ind]

    x_samples = [s[1] for s in samples.locations]
    y_samples = [s[2] for s in samples.locations]
    vals2coeffs = samples.vals2coeffs
    order_comb::Vector{Vector{Int}} = samples.order_comb
    n_coeffs::Int = size(vals2coeffs, 1)

    # get the function orders 
    range_vec::Vector{Vector{Float64}} = samples.range_vec
    # renormalize the range_vec to the relative range using the range_vec and delta 
    mean_range::Vector{Float64} = [(range_vec[i][1] + range_vec[i][2]) / 2 for i in 1:length(range_vec)]
    delta_range::Vector{Float64} = [(range_vec[i][2] - range_vec[i][1]) / 2 for i in 1:length(range_vec)]
    range_vec = [[mean_range[i] - delta_range[i] * rel_range, mean_range[i] + delta_range[i] * rel_range] for i in 1:length(range_vec)]

    # for every basis function, create a map of the function values in the range_vec interval 
    dx::Float64 = (range_vec[1][2] - range_vec[1][1]) / res
    dy::Float64 = (range_vec[2][2] - range_vec[2][1]) / res
    x_loc::Vector{Float64} = [range_vec[1][1] + dx / 2 + dx * i for i in 0:(res-1)]
    y_loc::Vector{Float64} = [range_vec[2][1] + dy / 2 + dy * i for i in 0:(res-1)]
    x_pos, x_samples, x_unit = transform_to_closest_unit(x_loc, x_samples, base_unit, factor)
    y_pos, y_samples, y_unit = transform_to_closest_unit(y_loc, y_samples, base_unit, factor)
    range_vec = [[x_pos[1], x_pos[end]], [y_pos[1], y_pos[end]]]
    if length(xlabel) == 0
        xlabel = L"$%$(samples.var_strs[1])$ [%$(x_unit)]"
    end
    if length(ylabel) == 0
        ylabel = L"$%$(samples.var_strs[2])$ [%$(y_unit)]"
    end

    max_order_x::Int = maximum([order_comb[i][1] for i in 1:n_coeffs])
    max_order_y::Int = maximum([order_comb[i][2] for i in 1:n_coeffs])
    x_vals_by_order::Vector{Vector{Float64}} = [fourier_exp2vals(x_loc, exponent, samples.range_vec[1]) for exponent in 0:max_order_x]
    y_vals_by_order::Vector{Vector{Float64}} = [fourier_exp2vals(y_loc, exponent, samples.range_vec[2]) for exponent in 0:max_order_y]


    # Get the operator 
    curr_values_, curr_t = values_at_t(op_str, solution, do_cumulants=do_cumulants, t=t)[[1,2]]
    curr_values = curr_values_[:,1] # make vector
    curr_t = curr_t[1]
    curr_t_rounded = round(curr_t * 10^1) / 10^1
    # determine coefficients for the function values 
    coeffs::Vector{ComplexF64} = vals2coeffs * curr_values
    z_vals::Matrix{ComplexF64} = zeros(res, res)
    for i in 1:n_coeffs
        z_vals .+= coeffs[i] * x_vals_by_order[order_comb[i][1]+1] * y_vals_by_order[order_comb[i][2]+1]'
    end
    # plot the function values
    op, op_string = get_op_and_string(op_str, opstr)
    parameters::Vector{String} = samples.var_strs
    prob_fun = samples.prob_fun
    weights_x = [prob_fun[1](x) for x in x_loc]
    weights_y = [prob_fun[2](y) for y in y_loc]
    weightmap::Matrix{Float64} = weights_x * weights_y'

    op_z_vals = op2.(op.(z_vals))
    mean_z_vals = (maximum(op_z_vals) + minimum(op_z_vals)) / 2
    delta_z = (maximum(op_z_vals) - minimum(op_z_vals)) / 2
    min_max_op_z_vals = [mean_z_vals - delta_z * rel_col_range, mean_z_vals + delta_z * rel_col_range]
    weighted_values, min_max_weights = weighted_heatmap(op_z_vals, weightmap, weighted=weighted, rel_col_range=rel_col_range, color_map=cmap)
    #println(range_vec)
    fig = Figure(size=figsize)
    if do_title
        ax = Axis(fig[1, 1], xlabel=xlabel, ylabel=ylabel, title=L"$t = %$curr_t_rounded$ %$base_unit_t")
    else
        ax = Axis(fig[1, 1], xlabel=xlabel, ylabel=ylabel)
    end
    image!(ax, (range_vec[1][1], range_vec[1][2]), (range_vec[2][1], range_vec[2][2]), weighted_values)
    xlims!(ax, range_vec[1])
    ylims!(ax, range_vec[2])
    # add colorbar of color_map
    ind = 2
    if weighted
        cbar2 = Colorbar(fig[1, 2], limits=(0.0, 1.0), colormap=:grays, label=L"Probability $p/p_{max}$")
        ind = 3
    end
    cbar = Colorbar(fig[1, ind], limits=(min_max_op_z_vals[1], min_max_op_z_vals[2]), colormap=cmap, label=L"%$op_prefix %$op_string %$op_suffix")
    # plot the samples as circles
    if do_locations
        scatter!(ax, x_samples, y_samples, color=:black, markersize=5, strokewidth=0.5, strokecolor=:white)
    end
    return fig, curr_t_rounded
end