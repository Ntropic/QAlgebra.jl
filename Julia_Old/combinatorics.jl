using Combinatorics
using Base.Threads
include("operator_terms.jl")
include("diff_Eq.jl")


function mrange(lengths::Vector{Int}; startpoint::Int=1)
    # computes the following type of nested loops in a single loop
    # for i in startpoint:lengths[1]
    #     for j in startpoint:lengths[2]
    #         # n times
    #         ...
    #     end
    # end
    curr_index::Vector{Int} = ones(Int, length(lengths)) * startpoint
    max_iterations = prod(lengths .+ (1 - startpoint))
    counter = 0
    chnl = Channel() do channel
        while counter < max_iterations
            put!(channel, copy(curr_index))
            counter += 1
            for i in length(lengths):-1:startpoint
                curr_index[i] += 1
                if curr_index[i] > lengths[i]
                    curr_index[i] = startpoint
                    if i == 1
                        return
                    end
                else
                    break
                end
            end
        end
    end
    return chnl
end

function mgreater_range(max_val::Int, n::Int; startpoint::Int=1, offdiag::Int=0)
    # computes the following type of nested loops in a single loop
    # for i in startpoint:max_val
    #     for j in (i+offdiag):max_val
    #         # n times
    #         ...
    #     end
    # end
    if n == 0
        return [[]]
    end
    curr_index::Vector{Int} = Vector{Int}(undef, n)
    curr_index[1] = startpoint
    for i in 2:n
        curr_index[i] = curr_index[i-1] + offdiag
    end
    #curr_index = curr_index[end:-1:1]
    max_vals::Vector{Int} = Vector{Int}(undef, n)
    max_vals[1] = max_val
    for i in 2:n
        max_vals[i] = max_vals[i-1] - offdiag
    end
    max_vals = max_vals[end:-1:1]
    changed::Bool = false
    chnl = Channel() do channel
        put!(channel, copy(curr_index))
        # change last element 
        while true
            if curr_index[n] < max_vals[n]
                curr_index[n] += 1
                put!(channel, copy(curr_index))
            else
                changed = false
                for i in (n-1):-1:1
                    if curr_index[i] + 1 <= max_vals[i]
                        curr_index[i] += 1
                        for j in (i+1):n
                            curr_index[j] = curr_index[j-1] + offdiag
                        end
                        put!(channel, copy(curr_index))
                        changed = true
                        break
                    end
                end
                if !changed
                    # last element is max_val but should be changed
                    return
                end
            end
        end
    end
    return chnl
end

function multi_mgreater_range(max_val::Int, n_s::Vector{Int}; startpoint::Int=1, offdiag::Int=0)
    # for every n in n_s construct the nested mgreater_range loop from right to left 
    # and return the channel of all combinations
    if length(n_s) == 0
        error("n_s must have length > 0")
    end
    chnl = Channel() do channel
        if length(n_s) == 1
            for curr_i in mgreater_range(max_val, n_s[1], startpoint=startpoint, offdiag=offdiag)
                put!(channel, copy([curr_i]))
            end
        else
            curr_k::Vector{Vector{Int}} = Vector{Vector{Int}}(undef, length(n_s))
            for curr_i in mgreater_range(max_val, n_s[1], startpoint=startpoint, offdiag=offdiag)
                curr_k[1] = curr_i
                for curr_j in multi_mgreater_range(max_val, n_s[2:end], startpoint=startpoint, offdiag=offdiag)
                    curr_k[2:end] = curr_j
                    put!(channel, copy(curr_k))
                end
            end
        end
    end
    return chnl
end

function multi_mgreater_range_flat(max_val::Int, n_s::Vector{Int}; startpoint::Int=1, offdiag::Int=0)
    # for every n in n_s construct the nested mgreater_range loop from right to left 
    # and return the channel of all combinations
    if length(n_s) == 0
        error("n_s must have length > 0")
    end
    chnl = Channel() do channel
        if length(n_s) == 1
            for curr_i in mgreater_range(max_val, n_s[1], startpoint=startpoint, offdiag=offdiag)
                put!(channel, copy(curr_i))
            end
        else
            last_entry::Int = 0
            for (i, curr_i) in enumerate(mgreater_range(max_val, n_s[1], startpoint=startpoint, offdiag=offdiag))
                for curr_j in multi_mgreater_range_flat(max_val, n_s[2:end], startpoint=startpoint, offdiag=offdiag)
                    put!(channel, append!(copy(curr_i), curr_j))
                end
            end
        end
    end
    return chnl
end

# function to generate all combinations of numbers from 1 to max_num with length n
function all_int_combinations(max_num::Int, n::Int)
    # uses mrange
    max_nums::Vector{Int} = max_num * ones(Int, n)
    all_combinations::Vector{Vector{Int}} = collect(mrange(max_nums))
    return all_combinations
end

# A function that constructs all combinations of increasing or equal numbers from 1 to max_num with length n
function geq_int_combinations(max_num::Int, n::Int)
    # uses mgreater_range
    all_combinations::Vector{Vector{Int}} = collect(mgreater_range(max_num, n, offdiag=0))
    return all_combinations
end
## Test 
#all_int_combinations(3, 2)
#geq_int_combinations(3, 2)

function all_vec_sum_n(max_sum::Int, len::Int)
    # Construct all integer vectors of integers > 0 and length len, so that the sum of the vector is <= max_sum
    all_comb::Vector{Vector{Int}} = []
    for i in 0:max_sum
        if len == 1
            push!(all_comb, [i])
        elseif len == 0
            push!(all_comb, [])
        else
            for comb in all_vec_sum_n(max_sum - i, len - 1)
                push!(all_comb, [i; comb])
            end
        end
    end
    return all_comb
end
## Test 
#max_sum = 3
#len = 2
#all_comb = all_vec_sum_n(max_sum, len)

#### Spin Combinations Construction ################################################################################################
function geq_spin_combinations(number_of_spins::Int; is_pauli::Bool=true)::Vector{String}
    # all combinations of 'x,y,z' for a given number of spins with (greater equal ordering)
    spin_types::Vector{String} = []
    if is_pauli
        spin_types = ["x", "y", "z"]
    else
        spin_types = ["p", "m", "z"]
    end
    ind_strings::Vector{String} = ["i", "j", "k", "l", "m", "n", "o", "p", "q", "r"]

    if number_of_spins > 10
        error("Too many spins. Maximum 10 possible.")
    end

    comb::Vector{String} = []
    if number_of_spins == 0
        # Add empty string
        push!(comb, "")
        return comb
    end

    ind_substring = ind_strings[1:number_of_spins]
    spin_ind_comb = geq_int_combinations(length(spin_types), number_of_spins)
    for i in 1:length(spin_ind_comb)
        curr_comb = spin_ind_comb[i]
        curr_string = ""
        for j in 1:number_of_spins
            curr_string *= spin_types[curr_comb[j]] * ind_substring[j]
            if j < number_of_spins
                curr_string *= "*"
            end
        end
        push!(comb, curr_string)
    end
    return comb
end
## Test
#geq_spin_combinations(2)

function bosonic_combinations_of_order(n::Int)::Tuple{Vector{String},Vector{Tuple{Int,Int}}}
    # Only consider bosonic operators of type $(a^\dagger)^p a^q$ with $q >= p$ (fun fact: p=q is hermitian)
    all_comb::Vector{String} = []
    pm_nums::Vector{Tuple{Int,Int}} = []
    for i in 0:n÷2
        push!(all_comb, "+"^i * "-"^(n - i))
        push!(pm_nums, (i, n - i))
    end
    return all_comb, pm_nums
end
## test
#bosonic_combinations_of_order(4)

function generate_all_combinations_of_combinations(orders::Vector{Int}, all_combinations_list::Vector{Vector{Vector{Int}}})
    # A function to construc t index combinations from a list of index combinations, for combinations of orders 
    len_orders::Int = length(orders)
    how_many_each::Vector{Int} = Vector{Int}(undef, len_orders)
    for i in 1:len_orders
        if orders[i] >= 1
            how_many_each[i] = length(all_combinations_list[orders[i]])
        else
            how_many_each[i] = 1
        end
    end
    combined_order::Int = sum(orders)
    curr_combination::Vector{Int} = Vector{Int}(undef, combined_order)
    curr_index::Int = 1
    chnl = Channel() do channel
        for i_s in mrange(how_many_each)
            curr_index = 1
            for i in 1:len_orders
                next_index = curr_index + orders[i]
                if orders[i] >= 1
                    curr_combination[curr_index:(next_index-1)] = all_combinations_list[orders[i]][i_s[i]]
                end
                curr_index = next_index
            end
            put!(channel, copy(curr_combination))
        end
    end
    return chnl
end

#### Generate the set of all operator products considered in an expansion #########################################################

###################################################################################################################################
#### Multi Spin Systems ###########################################################################################################
###################################################################################################################################

function multi_spin_gen_operator_strings(order::Int, max_spins::Int=0; is_pauli::Bool=true)::Vector{String}
    # spin_order is the order of spin operators,
    # max_cavity_order is the maximum number of bosonic operators in the system 
    if order < 0
        error("order must be positive.")
    end
    operators::Vector{String} = []
    min_c = max(0, order - max_spins) # minimum number of bosonic operators, so that the total order is order and there are at most max_spins bosonic operators
    for c in min_c:order
        cavity_combinations, pm_nums = bosonic_combinations_of_order(c)
        spin_order::Int = order - c
        spin_combinations = geq_spin_combinations(spin_order, is_pauli=is_pauli)
        curr_spin_pm_geq::Bool = false
        for spin_comb in spin_combinations
            num_p, num_m = 0, 0
            for s in spin_comb
                if s == 'p'
                    num_p += 1
                elseif s == 'm'
                    num_m += 1
                end
            end
            curr_spin_m_geq_p = num_p < num_m
            for (cavity_comb, pm_num) in zip(cavity_combinations, pm_nums)
                #println(spin_comb * cavity_comb, " - ", pm_num[1] == pm_num[2], " - ", curr_spin_pm_geq)
                pm_hermitian = pm_num[1] == pm_num[2]
                if !(!is_pauli && pm_hermitian && curr_spin_m_geq_p)
                    new_term = spin_comb * cavity_comb
                    if new_term != ""
                        push!(operators, new_term)
                    end
                end
            end
        end
    end
    return operators
end
function multi_spin_gen_operator_strings_up_to_order(max_order::Int, less_spin::Int=0; is_pauli::Bool=true)::Vector{String}
    # spin_order is the order of spin operators,
    # max_cavity_order is the maximum number of bosonic operators in the system 
    max_spins::Int = max_order - less_spin
    operators::Vector{String} = []
    for order in 1:max_order
        append!(operators, multi_spin_gen_operator_strings(order, max_spins, is_pauli=is_pauli))
    end
    return operators
end

function multi_spin_operators(max_order::Int, max_spins::Int=0; is_pauli::Bool=true)::Vector{Term}
    # create operator strings
    operator_strings::Vector{String} = multi_spin_gen_operator_strings(max_order, max_spins, is_pauli=is_pauli)
    # create operators
    operators::Vector{Term} = [make_term(op, is_pauli=is_pauli) for op in operator_strings]
    return operators
end
function multi_spin_operators_up_to_order(max_order::Int, less_spins::Int=0; is_pauli::Bool=true)::Vector{Term}
    # create operators
    operator_strings::Vector{String} = multi_spin_gen_operator_strings_up_to_order(max_order, less_spins, is_pauli=is_pauli)
    #if (max_order - less_spins > 4) && is_pauli == false  # no longer an issue, as indexing was changed, and operators are now separated by * or indexes by _
    #    println("Indexes reaching m, which is an operator")
    #end
    operators::Vector{Term} = [make_term(op, is_pauli=is_pauli) for op in operator_strings]
    return operators
end

function multi_spin_DE_operators(max_order::Int, max_spins::Int=0; is_pauli::Bool=true)::Vector{DE_Term_Multi}
    # create operator strings
    operators::Vector{Term} = multi_spin_operators(max_order, max_spins, is_pauli=is_pauli)
    # create DE operators
    DE_operators::Vector{DE_Term_Multi} = [time_derivative_expectation_of_operator_multispin(op) for op in operators]
    return DE_operators
end
function multi_spin_DE_operators_up_to_order(max_order::Int, less_spins::Int=0; is_pauli::Bool=true)::Vector{DE_Term_Multi}
    # create DE operators
    operators::Vector{Term} = multi_spin_operators_up_to_order(max_order, less_spins, is_pauli=is_pauli)
    DE_operators::Vector{DE_Term_Multi} = [time_derivative_expectation_of_operator_multispin(op) for op in operators]
    return DE_operators
end


# a function that takes a list and returns a dictionary with its elements as keys and their indexes as values
function list2dict(list::Vector)::Dict
    dict = Dict()
    for i in 1:length(list)
        dict[list[i]] = i
    end
    return dict
end

function specified_index_combinations_up_to_order(indexes::Vector{Int}, max_order::Int)::Tuple{Vector{Dict{Vector{Int},Int}},Vector{Vector{Vector{Int}}}}
    all_combinations_list::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}[]
    max_num::Int = length(indexes)
    for i in 1:max_order
        all_combinations = geq_int_combinations(max_num, i)
        # replace numbers by indexes 
        for comb in all_combinations
            for (j, c) in enumerate(comb)
                comb[j] = indexes[c]
            end
        end
        push!(all_combinations_list, all_combinations)
    end
    return all_combinations_list
end
## Test 
#dicter, lister = all_index_combinations_up_to_order(3, 2)
#dicter, lister = specified_index_combinations_up_to_order([2,3,4], 2)



###################################################################################################################################
#### Single Spin Systems ##########################################################################################################
###################################################################################################################################

function single_spin_gen_operator_strings_of_order(n::Int; is_pauli::Bool=true)::Vector{String}
    operators::Vector{String} = []
    pauli::Vector{String} = []
    if is_pauli
        pauli = ["x", "y", "z"]
    else
        pauli = ["p", "z"]
    end

    # All combinations of bosonic operators of order n
    all_comb_n = bosonic_combinations_of_order(n)[]
    append!(operators, all_comb_n)

    # All combinations of bosonic operators of order n-1, appended with each Pauli operator
    all_comb_n_1 = bosonic_combinations_of_order(n - 1)[1]
    for p in pauli
        append!(operators, [comb * p for comb in all_comb_n_1])
    end
    return operators
end
## test
#single_spin_gen_operator_strings_of_order(3)

function single_spin_gen_operator_strings_up_to_order(n::Int; separate_orders::Bool=true, is_pauli::Bool=true)::Union{Vector{String},Vector{Vector{String}}}
    operators::Any = separate_orders ? Vector{Vector{String}}() : Vector{String}() # = separate_orders ? Vector{String}[] : String[]
    if separate_orders
        for i in 1:n
            push!(operators, single_spin_gen_operator_strings_of_order(i, is_pauli=is_pauli))
        end
        return operators
    else
        for i in 1:n
            append!(operators, single_spin_gen_operator_strings_of_order(i, is_pauli=is_pauli))
        end
        return operators
    end
end
## test
#display(single_spin_gen_operator_strings_up_to_order(3, separate_orders=true))
#display(single_spin_gen_operator_strings_up_to_order(3, separate_orders=false))

function single_spin_gen_operators_of_order(n::Int; is_pauli::Bool=true)::Vector{Term}
    # create operator strings
    operator_strings::Vector{String} = single_spin_gen_operator_strings_of_order(n, is_pauli=is_pauli)
    # create operators
    operators::Vector{Term} = []
    for op in operator_strings
        push!(operators, make_term(op))
    end
    return operators
end
## test
#display(single_spin_gen_operators_of_order(3))

function single_spin_gen_operators_up_to_order(n::Int; is_pauli::Bool=true)::Vector{Vector{Term}}
    # create operators
    operators::Vector{Vector{Term}} = [single_spin_gen_operators_of_order(i, is_pauli=is_pauli) for i in 1:n]
    return operators
end
## test
#display(single_spin_gen_operators_up_to_order(3))

function single_spin_gen_DE_of_order(n::Int; is_pauli::Bool=true)::Vector{DE_Term}
    # create operator strings
    operators::Vector{Term} = single_spin_gen_operators_of_order(n, is_pauli=is_pauli)
    # create DE operators
    DE_operators::Vector{DE_Term} = []
    for op in operators
        push!(DE_operators, time_derivative_expectation_of_operator(op))
    end
    return DE_operators
end
## test
#display(single_spin_gen_DE_of_order(1))

function single_spin_gen_DE_up_to_order(n::Int; is_pauli::Bool=true)::Vector{Vector{DE_Term}}
    # create DE operators
    DE_operators::Vector{Vector{DE_Term}} = [single_spin_gen_DE_of_order(i, is_pauli=is_pauli) for i in 1:n]
    return DE_operators
end
## test
#eq_of_eq = single_spin_gen_DE_up_to_order(3)
#for eqs in eq_of_eq
#    for eq in eqs
#        display(eq)
#    end
#end