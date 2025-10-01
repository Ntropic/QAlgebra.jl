include("operator_terms.jl")
include("print_terms.jl")
include("pulses.jl")
include("cumulants.jl")
#using Distributed
using Base.Threads
using ForwardDiff
using DifferentialEquations
using SciMLSensitivity

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

mutable struct DifferentialProblem
    param_generator::Function
    all_eqs_indexed::Vector{DE_Term_indexed} # Add Vector{Multi_DE_Term_indexed}
    cumulant_vector::Vector{Cumulant_indexed}
    lower_order_cumulants_vector::Vector{Full_Cumulant_indexed}
    samples::SamplesAndWeights
    function DifferentialProblem(param_generator::Function, all_eqs_indexed::Vector{DE_Term_indexed}, cumulant_vector::Vector{Cumulant_indexed}, lower_order_cumulants_vector::Vector{Full_Cumulant_indexed}, samples::SamplesAndWeights)
        new(param_generator, all_eqs_indexed, cumulant_vector, lower_order_cumulants_vector, samples)
    end
    function DifferentialProblem(param_generator::Function, system::SystemAndDicts, samples::SamplesAndWeights)
        new(param_generator, system.all_eqs_indexed, system.cumulant_terms, system.lower_order_cumulant_terms, samples)
    end
end

@inline function print_time(print_dt, starttime, t, T, iterations, last_print_time; any_printing::Bool=false)
    do_print = false
    if print_dt > 0.0
        if (t - last_print_time > print_dt && (t - last_print_time) < (3 * print_dt)) 
            do_print = true
        end
    elseif print_dt==0.0
        do_print = true
    end
    curr_elapsed = -1.0
    if do_print
        # separate into minutes and seconds
        curr_time = time()
        elapsed_str, remaining_str = elapsed_remaining_time_str(starttime, t, T)
        if any_printing
            T_str = string("$(round(T, digits=1))")
            tstr = string("$(round(t, digits=1))")
            t_str = " "^max(0, length(T_str) - length(tstr)) * tstr
            #println("  ", t_str, "/", T_str, " ( Elapsed: ", elapsed_str, " - Remaining: ", remaining_str, " )")
            # Add iterations 
            println("  ", t_str, "/", T_str, " ( Elapsed: ", elapsed_str, " - Remaining: ", remaining_str, " ) | Its.: ", iterations)
        end
        last_print_time = copy(t)
    end
    return last_print_time, do_print
end
@inline function exceeded_limits_warning(op_exp_values, diff_problem::DifferentialProblem, t)
    print_warning = false
    max_val = 0.0
    for i in 1:length(op_exp_values)
        if diff_problem.all_eqs_indexed[i].clamped
            curr_abs = abs(op_exp_values[i])
            if curr_abs > 1.0
                print_warning = true
                if curr_abs > max_val
                    max_val = curr_abs
                end
            end
        end
    end
    if print_warning
        println("Warning: Some values are out of bounds at time ", t, " (max. abs. val. : ", max_val, ").")
    end
end

# To calculate cumulants for a set of initial states
@inline function calc_cumulants(initial_conditions, diff_problem::DifferentialProblem; threaded::Bool=true)
    n_cumulant_values::Int = length(diff_problem.cumulant_vector)
    cumulant_values::Vector{ComplexF64} = zeros(ComplexF64, n_cumulant_values)
    @inline function ext_calc_cumulant_value(cumulant_indexed::Cumulant_indexed)
        val::ComplexF64 = 0.0
        for i in 1:length(cumulant_indexed.weights)
            val += cumulant_indexed.weights[i] * index_conj_prod(initial_conditions, cumulant_indexed.partitions[i])
        end
        return val
    end
    @usethreads threaded for j in 1:n_cumulant_values
        cumulant_values[j] = ext_calc_cumulant_value(diff_problem.cumulant_vector[j])
    end
    return DifferentialSolutionAtT(initial_conditions, cumulant_values, 1)
end
# To calculate cumulants for a set of initial states
@inline function calc_cumulants(initial_conditions, system::SystemAndDicts; threaded::Bool=true)
    n_cumulant_values::Int = length(system.cumulant_terms)
    cumulant_values::Vector{ComplexF64} = zeros(ComplexF64, n_cumulant_values)
    @inline function ext_calc_cumulant_value(cumulant_indexed::Cumulant_indexed)
        val::ComplexF64 = 0.0
        for i in 1:length(cumulant_indexed.weights)
            val += cumulant_indexed.weights[i] * index_conj_prod(initial_conditions, cumulant_indexed.partitions[i])
        end
        return val
    end
    @usethreads threaded for j in 1:n_cumulant_values
        cumulant_values[j] = ext_calc_cumulant_value(system.cumulant_terms[j])
    end
    return DifferentialSolutionAtT(initial_conditions, cumulant_values, 1)
end

mutable struct DifferentialSolution
    T::Float64
    t::Vector{Float64}
    how_many_evaluations::Int
    vals::Vector{Vector{ComplexF64}}
    total_calc_time::Float64
    calc_times::Vector{Float64}
    calc_t_vals::Vector{Float64}
    calc_evals::Vector{Int}
    system::SystemAndDicts
    samples::SamplesAndWeights
    is_reduced::Bool
    function DifferentialSolution(T::Float64, t::Vector{Float64}, how_many_evaluations::Int, total_calc_time::Float64, calc_times::Vector{Float64}, calc_t_vals::Vector{Float64}, calc_evals::Vector{Int}, sol::ODESolution, system::SystemAndDicts, samples::SamplesAndWeights, is_reduced::Bool=false)
        vals = sol.u
        new(T, t, how_many_evaluations, vals, total_calc_time, calc_times, calc_t_vals, calc_evals, system, samples, is_reduced)
    end
end

#############################################################################################################
#### Solver for a set of initial coniditions ################################################################
function evolve_system(initial_conditions, diff_problem::DifferentialProblem, T::Float64; solvetype::Symbol=:Vern9, bounds::Bool=true, dt::Float64=0.05, n_t::Int=-100, abstol=1e-6, reltol=1.0, autodiff::Bool=true, threaded::Bool=true, print_dt::Float64=-1.0, any_printing::Bool=true, reduced_save::Bool=true, kwargs...)::DifferentialSolution
    # initial_conditions is a vector of initial conditions
    # diff_problem is a DifferentialProblem
    # T is a vector of times
    # solvetype is a symbol for the solver, currently only supports :nonstiff for Tsit5 and :stiff for Rodas5P, 
    #                   Beyond this: :Vern6, :Vern7, :Vern8, :Vern9 (default), (:Feagin10, :Feagin12, :Feagin14)
    # abstol and reltol are tolerances for the solver
    # autodiff is a boolean, that determines whether to use autodiff or not, currently uses ForwardDiffSensitivity()
    # save_cummulants = 0, 1, -1 determines whether to save cummulants or not, 0 = no, 1 = yes, -1 = only save the last one (default is 0 - don't save)
    # save_lower_cummulants = 0, 1, -1 determines whether to save lower cummulants or not, 0 = no, 1 = yes, -1 = only save the last one (default is 0 - don't save)
    # supports other solver arguments via kwargs (named tuples)
    # returns a vector of solutions
    if n_t > 0
        dt = T / (n_t - 1)
    else
        n_t::Int = ceil(Int, T / dt)+1
    end
    save_at::Vector{Float64} = range(0.0, stop=T, length=n_t)

    n_cumulant_values::Int = length(diff_problem.cumulant_vector)
    cumulant_values::Vector{ComplexF64} = Vector{ComplexF64}(undef, n_cumulant_values)

    last_print_time::Float64 = -print_dt * 1.5
    calc_times::Vector{Float64} = []
    calc_t_vals::Vector{Float64} = []
    calc_evals::Vector{Int} = []
    starttime = time()
    # print elapsed time via starttime in mm:ss
    how_many_evaluations::Int = 0
    # if save_cumulants and save_lower_cumulants are both = 1, add them to the initial_conditions vector, to also store them in the solution
    function timestep_generator!(d_op_exp_values, op_exp_values, diff_problem::DifferentialProblem, t)
        how_many_evaluations += 1
        if bounds # check if the values are within bounds for values with corresponding clamped=true in the operator_terms
            exceeded_limits_warning(op_exp_values, diff_problem, t)
            #for i in 1:length(op_exp_values)   # deactivated the clamping correction
            #    if diff_problem.all_eqs_indexed[i].clamped
            #        # mirror the values at the borders of the bounds
            #        if real(op_exp_values[i]) > 1.0
            #            op_exp_values[i] = 2.0 - op_exp_values[i]
            #        elseif real(op_exp_values[i]) < -1.0
            #            op_exp_values[i] = -2.0 - op_exp_values[i]
            #        end
            #        #op_exp_values[i] = max(-1.0, min(1.0, abs(op_exp_values[i]))) + 0.0im
            #    end
            #end
        end
        # First generated the parameters

        last_print_time, do_print = print_time(print_dt, starttime, t, T, how_many_evaluations, last_print_time, any_printing=any_printing)

        var_values = param_generator(t)
        @inline function calc_cumulant_value(cumulant_indexed::Cumulant_indexed)
            val::ComplexF64 = 0.0
            for i in 1:length(cumulant_indexed.weights)
                val += cumulant_indexed.weights[i] * index_conj_prod(op_exp_values, cumulant_indexed.partitions[i])
            end
            return val * cumulant_indexed.factor
        end
        @inline function get_term_val(term::DE_Term_indexed) #::ComplexF64
            val::ComplexF64 = 0.0
            curr_term::ComplexF64 = 0.0
            # first get linear_terms
            for lin_term in term.linear_terms
                # lin_term = (operator_index, conjugate, variables_index, coefficient)
                if !lin_term.conjugate
                    curr_term = op_exp_values[lin_term.operator_index]
                else
                    curr_term = conj(op_exp_values[lin_term.operator_index])
                end
                val += lin_term.coefficient * var_values[lin_term.variables_index] * curr_term
            end
            # then get cumulant_terms
            for cum_term in term.cumulant_terms
                # cum_term = (operator_index, conjugate, variables_index, coefficient)
                if !cum_term.conjugate
                    curr_term = cumulant_values[cum_term.operator_index]
                else
                    curr_term = conj(cumulant_values[cum_term.operator_index])
                end
                val += cum_term.coefficient * var_values[cum_term.variables_index] * curr_term
            end
            # now finally add constant_terms
            for const_term in term.constant_terms
                # const_term = (variables_index, coefficient)
                val += const_term.coefficient * var_values[const_term.variables_index]
            end
            return val
        end
        # Then update the cumulant_vector
        @usethreads threaded for i in 1:length(diff_problem.cumulant_vector)
            cumulant_values[i] = calc_cumulant_value(diff_problem.cumulant_vector[i])
        end
        # Then update the op_exp_values
        @usethreads threaded for i in 1:length(op_exp_values)
            d_op_exp_values[i] = get_term_val(diff_problem.all_eqs_indexed[i])
        end
        if do_print > 0.0
            push!(calc_times, time() - starttime)
            push!(calc_t_vals, t)
            push!(calc_evals, how_many_evaluations)
        end
    end

    problem = ODEProblem(timestep_generator!, initial_conditions, (0.0, T), diff_problem)
    #kwargs::NamedTuple = NamedTuple()
    if abstol > 0.0
        kwargs = (; kwargs..., abstol=abstol)
    end
    if reltol > 0.0
        kwargs = (; kwargs..., reltol=reltol)
    end
    kwargs = (; kwargs..., dt=dt/10, save_everystep=false, saveat=save_at) # removed dt, as it is the initial step size. and 

    if solvetype == :nonstiff 
        solution = solve(problem, Tsit5(); kwargs...)
    elseif solvetype == :Tsit5
        solution = solve(problem, Tsit5(); kwargs...)
    elseif solvetype == :Vern6
        solution = solve(problem, Vern6(); kwargs...)
    elseif solvetype == :Vern7
        solution = solve(problem, Vern7(); kwargs...)
    elseif solvetype == :Vern8
        solution = solve(problem, Vern8(); kwargs...)  # 
    elseif solvetype == :Vern9
        solution = solve(problem, Vern9(); kwargs...)
    elseif solvetype == :Feagin10
        solution = solve(problem, Feagin10(); kwargs...)
    elseif solvetype == :Feagin12
        solution = solve(problem, Feagin12(); kwargs...)
    elseif solvetype == :Feagin14
        solution = solve(problem, Feagin14(); kwargs...)
    else
        error("solvetype must be either :stiff (Rodas5P), :nonstiff (Tsit5), :Tsit5, :Vern6, :Vern7, :Vern8, :Vern9, :Feagin10, :Feagin12 or :Feagin14")
    end
    total_calc_time::Float64 = time() - starttime
    last_print_time, _ = print_time(0.0, starttime, T, T, how_many_evaluations, last_print_time, any_printing=any_printing)
    t::Vector{Float64} = solution.t
    if !reduced_save
        return DifferentialSolution(T, t, how_many_evaluations, total_calc_time, calc_times, calc_t_vals, calc_evals, solution, system, diff_problem.samples, reduced_save)
    else
        sys = reduced_system(system)    
        return DifferentialSolution(T, t, how_many_evaluations, total_calc_time, calc_times, calc_t_vals, calc_evals, solution, sys, diff_problem.samples, reduced_save)
    end
end

### Needs to be updated
struct DifferentialSolutionAtT{T}   # only the final element
    vals::Vector{T} # [op_exp_values, times]
    how_many_evaluations::Int
    samples::SamplesAndWeights
    function DifferentialSolutionAtT(vals::Vector{T}, cumulant_vals::Vector{T}, how_many_evaluations::Int, samples::SamplesAndWeights) where {T}
        new{T}(vals, how_many_evaluations, samples)
    end
end
function evolve_system_final(initial_conditions, diff_problem::DifferentialProblem, T::Float64; solvetype::Symbol=:Vern9, bounds::Bool=true, print_evaluations::Int=-1, abstol=1e-13, reltol=1.0, autodiff::Bool=true, threaded::Bool=true, print_dt::Float64=-1.0, kwargs...)::DifferentialSolutionAtT
    # initial_conditions is a vector of initial conditions
    # diff_problem is a DifferentialProblem
    # T is a vector of times
    # solvetype is a symbol for the solver, currently only supports :nonstiff for Tsit5 and :stiff for Rodas5P, 
    #                   Beyond this: :Vern6, :Vern7, :Vern8, :Vern9 (default), (:Feagin10, :Feagin12, :Feagin14)
    # abstol and reltol are tolerances for the solver
    # autodiff is a boolean, that determines whether to use autodiff or not, currently uses ForwardDiffSensitivity()
    # save_cummulants = 0, 1, -1 determines whether to save cummulants or not, 0 = no, 1 = yes, -1 = only save the last one (default is 0 - don't save)
    # save_lower_cummulants = 0, 1, -1 determines whether to save lower cummulants or not, 0 = no, 1 = yes, -1 = only save the last one (default is 0 - don't save)
    # supports other solver arguments via kwargs (named tuples)
    # returns a vector of solutions

    n_cumulant_values::Int = length(diff_problem.cumulant_vector)
    cumulant_values::Vector{ComplexF64} = Vector{ComplexF64}(undef, n_cumulant_values)
    n_lower_order_cumulant_values::Int = length(diff_problem.lower_order_cumulants_vector)

    last_print_time::Float64 = -print_dt * 1.5
    starttime = time()
    # print elapsed time via starttime in mm:ss
    how_many_evaluations::Int = 0
    # if save_cumulants and save_lower_cumulants are both = 1, add them to the initial_conditions vector, to also store them in the solution
    function timestep_generator!(d_op_exp_values, op_exp_values, diff_problem::DifferentialProblem, t)
        how_many_evaluations += 1
        if print_evaluations > 0
            if how_many_evaluations % print_evaluations == 0
                println("Evaluations: ", how_many_evaluations)
            end
        end
        if bounds # check if the values are within bounds for values with corresponding clamped=true in the operator_terms
            exceeded_limits_warning(op_exp_values, diff_problem, t)
            for i in 1:length(op_exp_values)
                if diff_problem.all_eqs_indexed[i].clamped
                    # mirror the values at the borders of the bounds
                    if real(op_exp_values[i]) > 1.0
                        op_exp_values[i] = 2.0 - op_exp_values[i]
                    elseif real(op_exp_values[i]) < -1.0
                        op_exp_values[i] = -2.0 - op_exp_values[i]
                    end
                    #op_exp_values[i] = max(-1.0, min(1.0, abs(op_exp_values[i]))) + 0.0im
                end
            end
        end
        # First generated the parameters
        last_print_time = print_time(print_dt, starttime, t, T, how_many_evaluations, last_print_time)
        var_values = param_generator(t)
        @inline function calc_cumulant_value(cumulant_indexed::Cumulant_indexed)
            val::ComplexF64 = 0.0
            for i in 1:length(cumulant_indexed.weights)
                val += cumulant_indexed.weights[i] * index_conj_prod(op_exp_values, cumulant_indexed.partitions[i])
            end
            return val
        end
        @inline function get_term_val(term::DE_Term_indexed) #::ComplexF64
            val::ComplexF64 = 0.0
            curr_term::ComplexF64 = 0.0
            # first get linear_terms
            for lin_term in term.linear_terms
                # lin_term = (operator_index, conjugate, variables_index, coefficient)
                if !lin_term.conjugate
                    curr_term = op_exp_values[lin_term.operator_index]
                else
                    curr_term = conj(op_exp_values[lin_term.operator_index])
                end
                val += lin_term.coefficient * var_values[lin_term.variables_index] * curr_term
            end
            # then get cumulant_terms
            for cum_term in term.cumulant_terms
                # cum_term = (operator_index, conjugate, variables_index, coefficient)
                if !cum_term.conjugate
                    curr_term = cumulant_values[cum_term.operator_index]
                else
                    curr_term = conj(cumulant_values[cum_term.operator_index])
                end
                val += cum_term.coefficient * var_values[cum_term.variables_index] * curr_term
            end
            # now finally add constant_terms
            for const_term in term.constant_terms
                # const_term = (variables_index, coefficient)
                val += const_term.coefficient * var_values[const_term.variables_index]
            end
            return val
        end
        # Then update the cumulant_vector
        @usethreads threaded for i in 1:length(diff_problem.cumulant_vector)
            cumulant_values[i] = calc_cumulant_value(diff_problem.cumulant_vector[i])
        end
        # Then update the op_exp_values
        @usethreads threaded for i in 1:length(op_exp_values)
            d_op_exp_values[i] = get_term_val(diff_problem.all_eqs_indexed[i])
        end
    end

    problem = ODEProblem(timestep_generator!, initial_conditions, (0.0, T), diff_problem)
    #kwargs::NamedTuple = NamedTuple()
    if abstol > 0.0
        kwargs = (; kwargs..., abstol=abstol)
    end
    if reltol > 0.0
        kwargs = (; kwargs..., reltol=reltol)
    end
    if dt > 0.0
        kwargs = (; kwargs..., dt=dt)
    end
    kwargs = (; kwargs..., save_everystep=false, save_end=true)
    if solvetype == :nonstiff
        solution = solve(problem, Tsit5(); kwargs...)
    elseif solvetype == :Vern6
        solution = solve(problem, Vern6(); kwargs...)
    elseif solvetype == :Vern7
        solution = solve(problem, Vern7(); kwargs...)
    elseif solvetype == :Vern8
        solution = solve(problem, Vern8(); kwargs...)  # 
    elseif solvetype == :Vern9
        solution = solve(problem, Vern9(); kwargs...)
    elseif solvetype == :Feagin10
        solution = solve(problem, Feagin10(); kwargs...)
    elseif solvetype == :Feagin12
        solution = solve(problem, Feagin12(); kwargs...)
    elseif solvetype == :Feagin14
        solution = solve(problem, Feagin14(); kwargs...)
    else
        error("solvetype must be either :stiff (Rodas5P), :nonstiff (Tsit5) or :tenth (Feagin10)")
    end
    last_print_time = print_time(1e-15, starttime, T, T, how_many_evaluations, last_print_time, force=true)
    # Calculate the cumulant terms
    curr_values::Vector{ComplexF64} = solution[end]
    return DifferentialSolutionAtT(curr_values, how_many_evaluations, diff_problem.samples)
end