using Combinatorics
include("cumulants.jl")

function which_times(solution::DifferentialSolution; t=-1)
    # which times? 
    # if t is Float -> evaluate at that time 
    # if t is Vector -> evaluate at all times in vector
    # if t is Int -> evaluate at as many equidistant time steps 
    # if t is -1.0 -> evaluate at t stored in solution
    sol_t::Vector{Float64} = solution.t
    sol_T::Float64 = solution.T
    if isa(t, Int)
        if t == -1
            return (1, length(sol_t)), false
        elseif t==-2
            return [sol_T], true
        elseif t < 0 || t > length(sol_t)
            error("t needs to be in [0, length(sol_t)], for integer indexing between [1, $length(sol_t)]")
        end
        return range(0, stop=sol_T, length=t), true
    elseif isa(t, Float64)
        if t < 0 || t > sol_T
            error("t needs to be in [0, T]")
        end
        return [t], true
    elseif isa(t, Vector)
        if isa(t[1], Float64)
            return t, true
        elseif isa(t[1], Int) && length(t) == 2
            return (t[1], t[2]), false
        end
    elseif isa(t, Tuple)
        if isa(t[1], Float64)
            # find indexes for the time interval 
            if length(t) < 2
                error("t needs to be an interval")
            end
            if t[1] > t[2]
                t[1], t[2] = t[2], t[1]
            end
            if t[1] < 0 || t[2] > sol_T
                error("t needs to be in [0, T]")
            end
            if length(t) == 2
                # find first and last index in interval 
                ind1 = findfirst(x -> x >= t[1], sol_t)
                ind2 = findfirst(x -> x >= t[2], sol_t)-1
                if ind1 == nothing || ind2 == nothing
                    error("t needs to be in [0, T] and not be an empty set")
                end
                return (ind1, ind2), false
            elseif length(t) == 3 && isa(t[3], Int)
                return range(t[1], stop=t[2], length=t[3]), true
            end 
        elseif isa(t[1], Int)
            return (t[1], t[2]), false
        end
    else
        error("t needs to be Int, Float64, Vector{Float64}, Tuple{Float64, Float64}, Tuple{Int, Int} or Vector{Int} of length 2")
    end
end
## Test 
#which_times(solution, t=-1)
function values_at_t(op_str::String, solution::DifferentialSolution; which_element::Int=-1, do_cumulants::Bool=false, t=-1, threaded::Bool=true, checks::Bool=false)::Tuple{Matrix{ComplexF64}, Vector{Float64}, Vector{Int}, Int}
    # Function to extract correct information from solution using system for a given operator string 
    # op_str is a string of the operator to extract, i.e. "+-z"
    # solution is a DifferentialSolution object
    # do_cumulants and do_lower_order_cumulants are booleans to decide if instead of expectation values cumulants or "full cumulants" of lower order should be plotted
    # check if op_str is key of cumulant_terms_dict or operator_strings_dict
    system = solution.system
    max_spin_order = system.spin_order 
    max_order = system.order 
    vals = solution.vals
    if checks
        check_consistent_operator_type(op_str, system.is_pauli)
        check_operator_ordering(op_str)
        only_know_conjugate(op_str)
        if length(op_str) > max_order+1
            error("Operator string $(op_str) is too long for the current system, neither in operator set nor one order higher (cumulant).")
        end
    end
    t_vals, do_interpolation = which_times(solution, t=t) # interpolation uses sol(t, idxs = indexes) instead of sol.u[t_ind][indexes] notation
    if !do_interpolation
        t_ind = t_vals[1]:t_vals[end]
        t_vals = solution.t[t_ind]
    else
        t_ind_lower = [findlast(x -> x <= t_val, solution.t) for t_val in t_vals]
        t_ind_upper = [findfirst(x -> x >= t_val, solution.t) for t_val in t_vals]
        dt_ind = t_ind_upper .- t_ind_lower
        t_lower = solution.t[t_ind_lower]
        t_upper = solution.t[t_ind_upper]
        t_coeff::Vector{Float64} = Vector{Float64}(undef, length(t_vals))
        for (i, t_l, t_u) in zip(1:length(t_vals), t_lower, t_upper)
            if t_u-t_l < 1e-10
                t_coeff[i] = 0.0
            else
                t_coeff[i] = (t_vals[i] - t_l) / (tu - t_l)
            end
        end
    end
    function sol(i::Int; idxs::Union{Vector{Int}, UnitRange{Int}})
        # linear interpolation of the solution
        return vals[t_ind_lower[i]][idxs] .+ t_coeff[i] * (vals[t_ind_upper[i]][idxs] .- vals[t_ind_lower[i]][idxs])
    end
    function sols(idxs::Vector{Int})
        # linear interpolation of the solution
        myval::Matrix{ComplexF64} = zeros(ComplexF64, length(idxs), length(t_lower))
        for i in 1:length(t_vals)
            myval[:,i] = vals[t_ind_lower[i]][idxs] .+ t_coeff[i] * (vals[t_ind_upper[i]][idxs] .- vals[t_ind_lower[i]][idxs])
        end
        return myval
    end
    function index_conj_prod_interp(indexes::Vector{Int})
        # calculates the product of the values that are indexed, with the extra condition, that negative indexes are conjugated
        res::Vector{ComplexF64} = zeros(ComplexF64, length(t_vals))
        for ind in indexes
            val = sols(idxs=ind)
            if ind < 0
                res .*= conj.(val)
            else
                res .*= val
            end
        end
        return res
    end
    function index_conj_prod_stored(indexes::Vector{Int})
        # calculates the product of the values that are indexed, with the extra condition, that negative indexes are conjugated
        res::Vector{ComplexF64} = zeros(ComplexF64, length(t_ind))
        for ind in indexes
            val = [vals[i][ind] for i in t_ind]
            if ind < 0
                res .*= conj.(val)
            else
                res .*= val
            end
        end
        return res
    end
    function calc_cumulant_value_interp(cumulant_indexed::Cumulant_indexed)
        val::Vector{ComplexF64} = zeros(ComplexF64, length(t_vals))
        for i in 1:length(cumulant_indexed.weights)
            val .+= cumulant_indexed.weights[i] .* index_conj_prod_interp(cumulant_indexed.partitions[i])
        end
        return val * cumulant_indexed.factor
    end
    function calc_full_cumulant_value_interp(cumulant_indexed::Full_Cumulant_indexed)
        val::Vector{ComplexF64} = zeros(ComplexF64, length(t_vals))
        for i in 1:length(cumulant_indexed.weights)
            val .+= cumulant_indexed.weights[i] .* index_conj_prod_interp(cumulant_indexed.partitions[i])
        end
        return val
    end
    function calc_cumulant_value_stored(cumulant_indexed::Cumulant_indexed)
        val::Vector{ComplexF64} = zeros(ComplexF64, length(t_ind))
        for i in 1:length(cumulant_indexed.weights)
            val .+= cumulant_indexed.weights[i] .* index_conj_prod_stored(cumulant_indexed.partitions[i])
        end
        return val * cumulant_indexed.factor
    end
    function calc_full_cumulant_value_stored(cumulant_indexed::Full_Cumulant_indexed)
        val::Vector{ComplexF64} = zeros(ComplexF64, length(t_ind))
        for i in 1:length(cumulant_indexed.weights)
            val .+= cumulant_indexed.weights[i] .* index_conj_prod_stored(cumulant_indexed.partitions[i])
        end
        return val
    end

    in_operators::Bool = false
    if haskey(system.operator_strings_dict, op_str)
        in_operators = true
        curr_op_group = system.operator_strings_dict[op_str]
    elseif haskey(system.cumulant_terms_dict, op_str)
        curr_op_group = system.cumulant_terms_dict[op_str]
    else
        error("Operator string $(op_str) not found in system.operator_strings_dict or system.cumulant_terms_dict")
    end
    first_index::Int = curr_op_group.first_index
    how_many::Int = curr_op_group.how_many
    if which_element == -1
        indexes = first_index:first_index+how_many-1
    else
        indexes = [first_index + which_element - 1]
    end
    curr_spin_orders::Vector{Int} = curr_op_group.op_spin_orders
    curr_solutions::Matrix{ComplexF64} = Matrix{ComplexF64}(undef, how_many, length(t_vals))
    is_cumulant::Int = 0  # 0 = not, 1 = cumulant, 2 = full cumulant
    if in_operators
        if !do_cumulants
            # just extract the values 
            if !do_interpolation
                @usethreads threaded for i in 1:length(t_ind) #(i, curr_t_ind) in enumerate(t_ind)
                    curr_t_ind = t_ind[i]
                    curr_solutions[:,i] = vals[curr_t_ind][indexes]
                end
            else
                @usethreads threaded for i in 1:length(t_vals) # (i, curr_t) in enumerate(t_vals)
                    curr_solutions[:,i] = sol(i, idxs=indexes)
                end
            end
        else # do lower order cumulant
            is_cumulant = 2
            if !do_interpolation
                @usethreads threaded for i in 1:length(indexes) 
                    index = indexes[i]
                    curr_cumulant::Full_Cumulant_indexed = system.lower_order_cumulant_terms[index]
                    curr_solutions[i,:] = calc_full_cumulant_value_stored(curr_cumulant)
                end
            else
                @usethreads threaded for i in 1:length(indexes)
                    index = indexes[i]
                    curr_cumulant::Full_Cumulant_indexed = system.lower_order_cumulant_terms[index]
                    curr_solutions[i,:] = calc_full_cumulant_value_interp(curr_cumulant)
                end
            end
        end
    else # do cumulant
        is_cumulant = 1
        if !do_interpolation
            @usethreads threaded for i in 1:length(indexes) 
                index = indexes[i]
                curr_cumulant::Cumulant_indexed = system.cumulant_terms[index]
                curr_solutions[i,:] = calc_cumulant_value_stored(curr_cumulant)
            end
        else
            @usethreads threaded for i in 1:length(indexes) 
                index = indexes[i]
                curr_cumulant::Cumulant_indexed = system.cumulant_terms[index]
                curr_solutions[i,:] = calc_cumulant_value_interp(curr_cumulant)
            end
        end
    end
    return curr_solutions, t_vals, curr_spin_orders, is_cumulant
end
## Test
#values_at_t("---zz", solution, do_cumulants=false, t=10, threaded=true)

function get_exp_and_std_to_vals_in_t(curr_solutions::Matrix{ComplexF64}, curr_spin_orders::Vector{Int}, samples::SamplesAndWeights)::Tuple{Vector{ComplexF64},Vector{ComplexF64}}
    # Function to calculate expectation values and standard deviations for a given operator string
    # curr_solutions is a matrix of the solutions for the operator string
    # curr_spin_orders is a vector of the spin orders for the operator string
    # samples is a SamplesAndWeights object
    n_t::Int = length(curr_solutions[1, :])
    curr_exp::Vector{ComplexF64} = Vector{ComplexF64}(undef, n_t)
    curr_std::Vector{ComplexF64} = Vector{ComplexF64}(undef, n_t)
    for j in 1:n_t
        curr_exp[j], curr_std[j] = Ex_and_Stdx(curr_solutions[:, j], samples, curr_spin_orders)
    end
    return curr_exp, curr_std
end


#### Samples Statistics for comparison
# function that calculates the weighted standard deviation
function weighted_mean_std(x::Vector{Float64}, w::Vector{Float64})
    mean_x = sum(x .* w) / sum(w)
    std_x = sqrt(sum((x .- mean_x) .^ 2 .* w) / sum(w))
    return mean_x, std_x
end
#  now for matrix row by row 
function weighted_mean_std(x::Matrix{Float64}, w::Vector{Float64})
    n_t::Int = size(x, 2)
    mean_vals::Vector{Float64} = zeros(n_t)
    std_vals::Vector{Float64} = zeros(n_t)
    for i in 1:n_t
        mean_vals[i], std_vals[i] = weighted_mean_std(x[:, i], w)
    end
    return mean_vals, std_vals
end
# Now the same for complex numbers
function weighted_mean_std(x::Vector{ComplexF64}, w::Vector{Float64})
    mean_x = sum(x .* w) / sum(w)
    std_x = sqrt(sum(abs2.(x .- mean_x) .* w) / sum(w))
    return mean_x, std_x
end
#  now for matrix row by row
function weighted_mean_std(x::Matrix{ComplexF64}, w::Vector{Float64})
    n_t::Int = size(x, 2)
    mean_vals::Vector{ComplexF64} = zeros(ComplexF64, n_t)
    std_vals::Vector{Float64} = zeros(n_t)
    for i in 1:n_t
        mean_vals[i], std_vals[i] = weighted_mean_std(x[:, i], w)
    end
    return mean_vals, std_vals
end

function get_op_and_string(op_str::String, opstr::String, dolatex::Bool=true)::Tuple{Function,String}
    op_str_options::Vector{String} = ["real", "imag", "abs", "abs2"]
    ops::Vector{Function} = [real, imag, abs, abs2]
    op_str_start_list::Vector{String} = ["\\mathrm{Re}\\,\\langle ", "\\mathrm{Im}\\,\\langle ", "|\\langle ", "|\\langle "]
    op_str_end_list::Vector{String} = ["\\rangle", "\\rangle", "\\rangle|" , "\\rangle|^2"]
    # get op_ind via op_str_options 
    op_ind::Int = -1
    for (i, op_str) in enumerate(op_str_options)
        if op_str == opstr
            op_ind = i
            break
        end
    end
    if op_ind == -1
        error("Operator string not found, options are " * string(op_str_options))
    end

    if dolatex
        substitutions = Dict('+' => raw"\hat{a}^\dagger", '-' => raw"\hat{a}", 'n' => raw"\hat{n}", 'x' => raw"\hat{\sigma}_{x}", 'y' => raw"\hat{\sigma}_{y}", 'z' => raw"\hat{\sigma}_{z}", 'p' => raw"\hat{\sigma}_{+}", 'm' => raw"\hat{\sigma}_{-}", 'O' => raw"\hat{\mathcal{O}}")
    else
        boson_substitutions = Dict('+' => raw"a†", '-' => raw"a")
        subscript_indexes = Dict('i' => "ᵢ", 'j' => "ⱼ", 'k' => "ₖ", 'l' => "ₗ", 'm' => "ₘ")
    end

    new_op_str::String = ""
    op_str2 = replace(op_str, "+-" => "n")  # replace "+-" with "n" in op_str
    ind_list = ["i", "j", "k", "l"]
    subscript_ind_list = ["ᵢ", "ⱼ", "ₖ", "ₗ"]
    last_ind = 1
    for i in 1:length(op_str2)
        if !dolatex
            new_op_str *= op_str2[i]
        else
            new_op_str *= substitutions[op_str2[i]]
        end
        if op_str2[i] in ['x', 'y', 'z', 'p', 'm']
            if !dolatex
                new_op_str *= subscript_ind_list[last_ind]
            else
                new_op_str *= "^{(" * ind_list[last_ind] * ")}"
                #new_op_str *= "_" * ind_list[last_ind]
            end
        end
    end

    op_string = op_str_start_list[op_ind] * new_op_str * op_str_end_list[op_ind]
    return ops[op_ind], op_string
end

# a function that takes a vector, and units, e.g. "Hz" a factor, e.g. 10^6, and a String and transforms the vector to the closest exponent 10^(3n) and the corresponding unit
function transform_to_closest_unit(v::Vector{Float64}, base_unit::String="Hz", factor=10^6; dolatex::Bool=false, return_scale::Bool=false)
    new_v::Vector{Float64} = v .* factor # rescale to base unit 
    # find unit of maximum exponent of abs(new_v)
    max_exp::Int = floor(Int, log10(maximum(abs.(new_v))))
    # closest multiple of 3 (floor)
    unit_strs::Vector{String} = ["n", "μ", "m", "", "k", "M", "G", "T"]
    if dolatex
        unit_strs = ["n", L"\mu", "m", "", "k", "M", "G", "T"]
    end
    exp_ind::Int = 0
    if max_exp > 0
        exp_ind = max(1, min(8, max_exp ÷ 3 + 4))  # 1 is n, 4 is "", 8 is T
    else
        exp_ind = floor(Int, max_exp / 3) + 4
    end
    new_factor::Float64 = 10.0^(3 * (exp_ind - 4))
    new_unit::String = unit_strs[exp_ind] * base_unit
    new_v = new_v ./ new_factor
    if !return_scale
        return new_v, new_unit
    else
        return new_v, new_unit, factor / new_factor
    end
end
# Version with complex v -> use norm 
function transform_to_closest_unit(v::Vector{ComplexF64}, base_unit::String="Hz", factor=10^6; dolatex::Bool=false, return_scale::Bool=false)
    return transform_to_closest_unit(norm.(v), base_unit, factor; dolatex=dolatex, return_scale=return_scale)
end

function transform_to_closest_unit(v::Vector{Float64}, w::Vector{Float64}, base_unit::String="Hz", factor=10^6; dolatex::Bool=false)
    # same but for two vectors v and w
    if base_unit == "Hz"
        factor /= 2 * pi
    end
    new_v::Vector{Float64} = v .* factor # rescale to base unit 
    new_w::Vector{Float64} = w .* factor # rescale to base unit
    # find unit of maximum exponent of abs(new_v)
    max_val = max(maximum(abs.(new_v)), maximum(abs.(new_w)))
    max_exp::Int = floor(Int, log10(max_val))
    # closest multiple of 3 (floor)
    unit_strs::Vector{String} = ["n", "μ", "m", "", "k", "M", "G", "T"]
    if dolatex
        unit_strs = ["n", L"\mu", "m", "", "k", "M", "G", "T"]
    end
    exp_ind::Int = 0
    if max_exp > 0
        exp_ind = max(1, min(8, max_exp ÷ 3 + 4))  # 1 is n, 4 is "", 8 is T
    else
        exp_ind = floor(Int, max_exp / 3) + 4
    end
    new_factor::Float64 = 10.0^(3 * (exp_ind - 4))
    new_unit::String = unit_strs[exp_ind] * base_unit
    new_v = new_v ./ new_factor
    new_w = new_w ./ new_factor
    return new_v, new_w, new_unit
end
# Version with complex v -> use norm, and w -> use norm
function transform_to_closest_unit(v::Vector{ComplexF64}, w::Vector{ComplexF64}, base_unit::String="Hz", factor=10^6; dolatex::Bool=false)
    return transform_to_closest_unit(norm.(v), norm.(w), base_unit, factor; dolatex=dolatex)
end
