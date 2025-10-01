 # Construct the commulant expansion of the term
# assuming that all correlations are already captured by the lower order terms, the cumulant expansion of an expectation values <X1 X2 ... Xn>_c = 0, so that we can expand <X1 X2 ... Xn> into
# <X1 X2 ... Xn> = \sum_{p in P(I)} (|p|-1)! (-1)^{|p|} \prod_{B in p} <\prod_{i \in B} X_i>
# where |p| is the number of terms in the partition and P(I) is the set of all partitions of the set I = {1, 2, ..., n},
# and B is a block in the partition p (i.e. a subset of I), with the i in B the indices of the operators in the block B
# (notice, that the operators don't get reordered)
# Example: <X1 X2 X3> = <X1 X2><X3> + <X1 X3><X2> + <X2 X3><X1> - 2<X1><X2><X3>
include("operator_terms.jl")
include("preprocessing.jl")
include("indexing.jl")
using Combinatorics
using StatsBase

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

# a function that constructs all combinations of numbers (1 to n-1) that sum up to n
# only combinations from large to small (decreasing but not strictly decreasing)
function combinations_summing_to_n(n::Int, max_i::Int=-1)::Vector{Vector{Int}}
    # which combinations of numbers sum up to n
    # for example, for n=3, we have [[1, 1, 1], [2, 1], [3]] # from large to small, to avoid double counting
    # n is the number we want to sum up to
    # max_i is the smallest number to the right (as we want the numbers to decrease it sets a limit to following numbers)
    # returns a Vector of Vectors, where each Vector is a combination of numbers that sum up to n (Vector{Vector{Int}})
    all_combinations::Vector{Vector{Int}} = Vector{Vector{Int}}()
    if max_i == -1
        max_i = n
    elseif max_i > n
        max_i = n
    end
    remaining::Int = n
    curr_combination::Array{Int,1} = Array{Int,1}()
    for i in 1:max_i
        remaining = n - i
        if remaining == 0
            push!(all_combinations, [i])
        else
            # get all combinations of remaining numbers
            remaining_combinations = combinations_summing_to_n(remaining, i)
            for j in 1:length(remaining_combinations)
                curr_combination = [i]
                append!(curr_combination, remaining_combinations[j])
                push!(all_combinations, curr_combination)
            end
        end
    end
    return all_combinations
end
## Test
#combinations_summing_to_n(3)

# a function that returns modified combinations_summing_to_n,
# if multiple elements have the same value, we group them in a sublist
# for example, for n=3, we turn [[1, 1, 1], [2, 1], [3]] into [[3, 1]], [[1,2], [1,1]], [[1,3]]]   ([counts, value]])
# or [2,1,1] into [[1,2], [2, 1]]
# as each expectation value is a scalar, their order does not matter, and we don't wish to double count,
# so we group them, and then generate their value combinations for the complete groups, to ensure single ocunting
function grouped_combinations_summing_to_n(n::Int, max_i::Int=-1; including_n::Bool=false)::Vector{Vector{Tuple{Int,Int}}}
    # get all combinations of integers summing to n
    if including_n
        all_combinations = combinations_summing_to_n(n, max_i)
    else
        all_combinations = combinations_summing_to_n(n, max_i)[1:end-1]   # remove last one - all in one
    end
    # group the combinations
    grouped_combinations::Vector{Vector{Tuple{Int,Int}}} = Vector{Vector{Tuple{Int,Int}}}()
    curr_comb::Vector{Tuple{Int,Int}} = Tuple{Int,Int}[]
    for comb in all_combinations
        # check if any of the unique values in comb appear more than once
        # get unique values and their counts
        unique_values = sort(unique(comb), rev=true)
        counts = [count(x -> x == i, comb) for i in unique_values]
        # if all counts are 1, then we don't need to group
        curr_comb = Tuple{Int,Int}[]
        for (val, counts) in zip(unique_values, counts)
            push!(curr_comb, (counts, val))
        end
        push!(grouped_combinations, curr_comb)
    end
    return grouped_combinations
end
## Test
#grouped_combinations_summing_to_n(3)

# A function that splits n elements into k equally sized groups, 
# but in all combinations of splits, where different orders of the groups are not counted as different 
# and different orders within the groups also do not count as different
function all_equal_array_splits(arr::Vector{Int}, size_of_boxes::Int)::Vector{Vector{Vector{Int}}}
    n::Int = length(arr)
    if n % size_of_boxes != 0
        error("Array length not divisible by size_of_boxes")
    end
    number_of_boxes::Int = Int(n / size_of_boxes)  # number of boxes
    if number_of_boxes == 1
        return [[arr]]
    end
    # get all combinations of numbers that sum up to number_of_boxes
    # first box starts with first elements, second box starts with second remaining element, etc.
    curr_boxing::Vector{Vector{Int}} = Vector{Vector{Int}}()
    all_boxings::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}()
    curr_box::Vector{Int} = Vector{Int}()
    anti_array::Vector{Int} = Vector{Int}()
    for comb in combinations(arr[2:end], size_of_boxes - 1)
        curr_box = [arr[1]]
        append!(curr_box, comb)
        anti_array = setdiff(arr, curr_box)
        subsequent_boxes::Vector{Vector{Vector{Int}}} = all_equal_array_splits(anti_array, size_of_boxes)
        for sub_box in subsequent_boxes
            curr_boxing = Vector{Vector{Int}}()
            push!(curr_boxing, curr_box)
            for sub in sub_box
                push!(curr_boxing, sub)
            end
            push!(all_boxings, curr_boxing)
        end
    end
    return all_boxings
end
## Test
#all_equal_array_splits([1,2,3,4], 2)

# A function that partitions the numbers from 1 to n into all possible partitions (<X1 X2><X3> -> function output is [[1, 2], [3]], for <X1 X2 X3> -> function output is [[1, 2, 3]] and <X1><X2><X3> -> function output is [[1], [2], [3]])
function partitions(n_array::Vector{Int}, comb::Union{Vector{Vector{Tuple{Int,Int}}},Nothing}=nothing; including_n::Bool=false)::Vector{Vector{Vector{Int}}}
    n = length(n_array)
    if isa(comb, Nothing)
        comb = grouped_combinations_summing_to_n(n, n, including_n=including_n)
    end
    all_partitions::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}()
    new_partition::Vector{Vector{Int}} = Vector{Vector{Int}}()
    distributed_partitions::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}()
    for c in comb
        c0 = c[1]
        for partition in combinations(n_array, c0[1] * c0[2])
            # get anti partition for remaining elements
            distributed_partitions = all_equal_array_splits(partition, c0[2])
            if length(c) > 1
                # get anti partition for remaining elements
                remaining_elements = setdiff(n_array, partition)
                remaining_partition = partitions(remaining_elements, [c[2:end]])
                for p in remaining_partition
                    # concatenate distributed_partition and p into new_partition
                    for i in 1:length(distributed_partitions)
                        new_partition = deepcopy(distributed_partitions[i])
                        for i in 1:length(p)
                            push!(new_partition, p[i])
                        end
                        push!(all_partitions, new_partition)
                    end
                end
            else
                for partition in distributed_partitions
                    push!(all_partitions, partition)
                end
            end
        end
    end
    return all_partitions
end
## Test
#partitions([1,2,3,4])

mutable struct Cumulant
    partitions::Vector{Vector{Vector{Int}}}   # Cumulant terms 
    weights::Vector{Int}                      # Weights of the cumulant terms 
    operator::Vector{Int}                     # The operator for which we expanded the cumulant
end

mutable struct Full_Cumulant
    partitions::Vector{Vector{Vector{Int}}}   # Cumulant terms 
    weights::Vector{Int}                      # Weights of the cumulant terms 
    operator::Vector{Int}                     # The operator for which we expanded the cumulant
end

function subscript_int(num::Int)::String
    subscript_dict::Dict{Char,String} = Dict('1' => "₁", '2' => "₂", '3' => "₃", '4' => "₄", '5' => "₅", '6' => "₆", '7' => "₇", '8' => "₈", '9' => "₉", '0' => "₀", '-' => "₋", '+' => "₊")
    num_str::String = string(num)
    sub_str::String = ""
    for c in num_str
        sub_str *= subscript_dict[c]
    end
    return sub_str
end

function cumulant_braket_to_str(t::Vector{Int}; do_latex::Bool=true, operator_names::Union{String,Vector{String}}="\\hat{X}")::String
    term_str::String = ""
    name_is_string::Bool = isa(operator_names, String)
    for i in 1:length(t)
        if name_is_string
            if do_latex
                term_str *= operator_names * "_" * string(t[i])
            else
                term_str *= operator_names * subscript_int(t[i])
            end
        else
            term_str *= operator_names[t[i]]
        end
    end
    if do_latex
        term_str = "\\braket{" * term_str * "}"
    else
        #term_str = "〈" * term_str * "〉"
        term_str = "⟨" * term_str * "⟩"
    end
    return term_str
end

function cumulant_term_to_str(term::Vector{Vector{Int}}, weight::Int; do_latex::Bool=true, operator_names::Union{String,Vector{String}}="\\hat{X}")::String
    str::String = ""
    name_is_string::Bool = isa(operator_names, String)
    if !name_is_string
        cum_length = 0
        for t in term
            cum_length += length(t)
        end
        if length(operator_names) < cum_length
            error("Length of operator_names does not match number of operators in term")
        end
    end
    if weight >= 0
        str *= "+"
    else
        str *= "-"
    end
    if abs(weight) != 1
        str *= string(abs(weight))
    end
    for t in term
        str *= cumulant_braket_to_str(t, do_latex=do_latex, operator_names=operator_names)
    end
    return str
end

function cumulant_to_str(cumulant::Union{Cumulant,Full_Cumulant}; do_latex::Bool=true, operator_names::Union{String,Vector{String}}="\\hat{X}")::String
    str::String = ""
    curr_str::String = ""
    is_cumulant = isa(cumulant, Cumulant)
    str = cumulant_braket_to_str(cumulant.operator, do_latex=do_latex, operator_names=operator_names)
    if is_cumulant
        str *= " = "
    else
        if do_latex
            str *= "_c = "
        else
            str *= "ᶜ = "
        end
    end
    for i in 1:length(cumulant.partitions)
        curr_str = cumulant_term_to_str(cumulant.partitions[i], cumulant.weights[i], do_latex=do_latex, operator_names=operator_names)
        if i == 1 && curr_str[1] == '+'
            str *= " " * curr_str[2:end]
        else
            if curr_str[1] == '+'
                str *= " + " * curr_str[2:end]
            else
                str *= " - " * curr_str[2:end]
            end
        end
    end
    return str
end

function Base.show(io::IO, ::MIME"text/plain", cumulant::Cumulant)
    str = cumulant_to_str(cumulant, do_latex=false)
    print(io, str)
end
function Base.show(io::IO, ::MIME"text/latex", cumulant::Cumulant)
    str = latexstring(cumulant_to_str(cumulant, do_latex=true))
    print(io, str)
end
function Base.show(io::IO, ::MIME"text/plain", cumulant::Full_Cumulant)
    str = cumulant_to_str(cumulant, do_latex=false)
    print(io, str)
end
function Base.show(io::IO, ::MIME"text/latex", cumulant::Full_Cumulant)
    str = latexstring(cumulant_to_str(cumulant, do_latex=true))
    print(io, str)
end

function cumulant_expansion(n_array::Vector{Int}; including_n::Bool=false)::Union{Cumulant,Full_Cumulant}
    # Returns partitions and a vector of their weights, to construct the cumulant expansion
    # n_array is an array of numbers from 1 to n, where n is the number of operators in the term
    # returns a tuple of partitions and their weights
    # partitions is a vector of vectors of vectors, where each vector is a partition of the numbers from 1 to n
    # weights is a vector of integers, where each integer is the weight of the corresponding partition
    # the weight of a partition is (-1)^|p| (|p|-1)! where |p| is the number of terms in the partition
    # for example, for <X1 X2 X3> = <X1 X2><X3> + <X1 X3><X2> + <X2 X3><X1> - 2<X1><X2><X3>
    # we have partitions = [[[1, 2], [3]], [[1, 3], [2]], [[2, 3], [1]], [[1], [2], [3]]]
    # and weights = [1, 1, 1, -2]
    all_partitions::Vector{Vector{Vector{Int}}} = partitions(n_array, including_n=including_n)
    lengths::Vector{Int} = [length(p) for p in all_partitions]
    if including_n
        weights = Int[(-1)^(l + 1) * factorial(l - 1) for l in lengths]
        cumulant = Full_Cumulant(all_partitions, weights, n_array)
    else
        weights = Int[(-1)^l * factorial(l - 1) for l in lengths]
        cumulant = Cumulant(all_partitions, weights, n_array)
    end
    return cumulant
end
## Test
#cumulant = cumulant_expansion([1,2,3])

function n_th_order_cumulant(n::Int; include_n::Bool=false)::Union{Cumulant,Full_Cumulant}
    arr::Vector{Int} = [i for i in 1:n]
    return cumulant_expansion(arr, including_n=include_n)
end
## Test
#n_th_order_cumulant(3)

##### Scaling
#for n in 2:6
#    println(n, " & ", length(partitions(collect(1:n))), "  \\\\")
#end

##### Print the cumulants 
#for n in 2:6
#    cumulant = cumulant_expansion(collect(1:n))
#    display(cumulant)
#    #println(cumulant_to_str(cumulant))
#end


# Similar as Cumulant, but indexed with respect to operator indexes -> only Vector of Vector{Int}
struct Cumulant_indexed
    partitions::Vector{Vector{Int}}           # Cumulant terms   # negative signs for complex conjugate
    weights::Vector{Int}                      # Weights of the cumulant terms 
    #operator::String                            # The operator for which we expanded the cumulant
    factor::Float64
    clamped::Bool
end
struct Full_Cumulant_indexed
    partitions::Vector{Vector{Int}}           # Cumulant terms  # negative sign for complex conjugate
    weights::Vector{Int}                      # Weights of the cumulant terms 
    #operator::String                            # The operator for which we expanded the cumulant
    clamped::Bool
end

#### Now for reduced cumulants ##################################################################
function index_conj_prod(values::Vector{ComplexF64}, indexes::Vector{Int})
    # calculates the product of the values that are indexed, with the extra condition, that negative indexes are conjugated
    res::ComplexF64 = 1.0
    for ind in indexes
        val = values[abs(ind)]
        if ind < 0
            res *= conj(val)
        else
            res *= val
        end
    end
    return res
end
## Test
#values = rand(rng, 10) + im*rand(rng, 10)
#indexes = [1, -2]
#display(index_conj_prod(values, indexes))

#### Calculate the reduced cumulant terms for a given set of cumulant terms
function get_cumulant_value(cumulant_indexed::Cumulant_indexed)::ComplexF64
    #if !isdefined(Main, :op_exp_values)
    #    error("op_exp_values not defined. (can be provided as an input argument with the get_cumulant_value_provided function. ")
    #end
    val::ComplexF64 = 0.0
    for i in 1:length(cumulant_indexed.weights)
        val += cumulant_indexed.weights[i] * index_conj_prod(op_exp_values, cumulant_indexed.partitions[i])
    end
    return val*cumulant_indexed.factor
end
function get_cumulant_value(cumulant_indexed::Full_Cumulant_indexed)::ComplexF64
    #if !isdefined(Main, :op_exp_values)
    #    error("op_exp_values not defined. (can be provided as an input argument with the get_cumulant_value_provided function. ")
    #end
    val::ComplexF64 = 0.0
    for i in 1:length(cumulant_indexed.weights)
        val += cumulant_indexed.weights[i] * index_conj_prod(op_exp_values, cumulant_indexed.partitions[i])
    end
    return val
end
function get_cumulant_value_provided(cumulant_indexed::Cumulant_indexed, op_exp_values)::ComplexF64
    #if !isdefined(Main, :op_exp_values)
    #    error("op_exp_values not defined. (can be provided as an input argument with the get_cumulant_value_provided function. ")
    #end
    val::ComplexF64 = 0.0
    for i in 1:length(cumulant_indexed.weights)
        val += cumulant_indexed.weights[i] * index_conj_prod(op_exp_values, cumulant_indexed.partitions[i])
    end
    return val*cumulant_indexed.factor
end
function get_cumulant_value(cumulant_indexed::Full_Cumulant_indexed, op_exp_values)::ComplexF64
    #if !isdefined(Main, :op_exp_values)
    #    error("op_exp_values not defined. (can be provided as an input argument with the get_cumulant_value_provided function. ")
    #end
    val::ComplexF64 = 0.0
    for i in 1:length(cumulant_indexed.weights)
        val += cumulant_indexed.weights[i] * index_conj_prod(op_exp_values, cumulant_indexed.partitions[i])
    end
    return val
end

function get_cumulant_values(cumulant_indexed_vector::Vector{Cumulant_indexed})::Vector{ComplexF64}
    len = length(cumulant_indexed_vector)
    values::Vector{ComplexF64} = Vector{ComplexF64}(undef, len)
    for i in 1:len
        values[i] = get_cumulant_value(cumulant_indexed_vector[i])
    end
    return values
end
function get_cumulant_values(cumulant_indexed_vector::Vector{Full_Cumulant_indexed})::Vector{ComplexF64}
    len = length(cumulant_indexed_vector)
    values::Vector{ComplexF64} = Vector{ComplexF64}(undef, len)
    for i in 1:len
        values[i] = get_cumulant_value(cumulant_indexed_vector[i])
    end
    return values
end
function get_cumulant_values_provided(cumulant_indexed_vector::Vector{Cumulant_indexed}, op_exp_values)::Vector{ComplexF64}
    len = length(cumulant_indexed_vector)
    values::Vector{ComplexF64} = Vector{ComplexF64}(undef, len)
    for i in 1:len
        values[i] = get_cumulant_value_provided(cumulant_indexed_vector[i], op_exp_values)
    end
    return values
end

function cumulant_order_from_term_str(term_str::String)::Int
    # count the sum of the number of occurances of x,y,z,p,m
    # if there is at least one + or - in the term, then the order is that plus 1 (cavity operators are not separated)
    # if there is no + or - in the term, then the order is the sum of the number of occurances of x,y,z,p,m

    # count the number of x,y,z,p,m
    counter::Int = 0
    for s in term_str
        if s in ['x', 'y', 'z', 'p', 'm']
            counter += 1
        end
    end
    if term_str[1] in ['+', '-']
        counter += 1
    end
    return counter
end

###################################################################################################
##### Higher Order Reduced Cumulants ##############################################################
###################################################################################################

function prepare_cumulant_to_cumulant_indexed(term_str::String, curr_op_spin_orders::Vector{Int}, op_index_dicts::Dict{String,Op_Group_Type}, cumulant::Cumulant, lower_cumulant::Cumulant, max_spin_order::Int)::Tuple{String,String,Vector{Int}, Float64,Vector{Vector{Int}},Vector{Vector{Vector{Vector{Int}}}},Vector{Vector{Bool}}}
    how_many_operators = cumulant_order_from_term_str(term_str)
    how_many_spin_operators::Int = sum(curr_op_spin_orders)
    how_many_cavity_operators::Int = length(term_str) - how_many_spin_operators
    any_cavity_operators::Int = 0
    do_averaged_cumulants::Bool = false
    if how_many_spin_operators > max_spin_order && how_many_cavity_operators > 0
        do_averaged_cumulants = true
    end
    if how_many_cavity_operators > 0
        any_cavity_operators = 1
    end
    #how_many_cavity_operators::Int = how_many_operators - how_many_spin_operators
    partitions::Vector{Vector{Vector{Int}}} = cumulant.partitions
    lower_partitions::Vector{Vector{Vector{Int}}} = lower_cumulant.partitions
    weights::Vector{Int} = cumulant.weights
    lower_weights::Vector{Int} = lower_cumulant.weights 
    # create X_term_str which replaces the cavity substring (all + and -) with a single X
    X_term_strs::Vector{String} = []
    cavity_str::String = term_str[1:how_many_cavity_operators]
    spin_str::String = term_str[how_many_cavity_operators+1:end]
    X_strs::Vector{String} = []
    if !do_averaged_cumulants
        curr_X_term_str = ""
        if how_many_cavity_operators > 0
            curr_X_term_str *= "X"
        end
        curr_X_term_str *= spin_str
        push!(X_term_strs, curr_X_term_str)
        push!(X_strs, "")
        curr_partitions = partitions
        curr_weights = weights
        curr_factor = 1.0
    else
        for i in 1:how_many_spin_operators
            curr_X_term_str = spin_str[1:i-1]*"X"*spin_str[i+1:end]
            push!(X_term_strs, curr_X_term_str)
            push!(X_strs, string(spin_str[i]))
        end
        curr_partitions = lower_partitions
        # repeat lower_weights how_many_spin_operators times
        curr_weights = repeat(lower_weights, how_many_spin_operators)
        curr_factor = 1/how_many_spin_operators
    end

    first_indexes::Vector{Vector{Int}} = []
    index_subsets_by_op_type::Vector{Vector{Vector{Vector{Int}}}} = []
    conjugate_terms::Vector{Vector{Bool}} = []
    for (X_term_str, X_str) in zip(X_term_strs, X_strs)
        X_ind::Int = 0
        # if X in X_term_str, then X_ind is the index of it's position in X_term_str 
        if length(X_term_str) > 0 
            for i in 1:length(X_term_str)
                if X_term_str[i] == 'X'
                    X_ind = i
                    break
                end
            end
        end
        for partition in curr_partitions
            curr_part_indexes::Vector{Int} = Vector{Int}(undef, length(partition))
            curr_do_conjugate::Vector{Bool} = Vector{Bool}(undef, length(partition))
            curr_index_subsets_by_op_type::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}(undef, length(partition))
            for (i, part) in enumerate(partition)
                # if 1 in part, then X needs to replaced with X_str 
                if maximum(part) > length(X_term_str)
                    error("Error: index larger than length of X_term_str (", X_term_str, ") in part ", part, " of partition ", partition, " of cumulant ", cumulant)
                end
                curr_str::String = ""
                for p in part
                    if p == X_ind 
                        curr_str = cavity_str
                        break
                    end
                end
                for p in part
                    if p != X_ind
                        curr_str *= X_term_str[p]
                    else
                        curr_str *= X_str
                    end
                end
                curr_do_conjugate[i], new_str = conjugate_to_basis(curr_str)
                # check if we need to conjugate the string
                #println(curr_str, " -> ", new_str)
                curr_op::Op_Group_Type = op_index_dicts[new_str]
                curr_part_indexes[i] = curr_op.first_index
                curr_how_many_op_each = curr_op.op_spin_orders
                cum_how_many_each = cumsum(curr_how_many_op_each)
                # only spin operators need indexing subpart is the indexes of part that are larger than how_many_cavity_operators
                subpart::Vector{Int} = []
                if !do_averaged_cumulants
                    subpart = part[part.>any_cavity_operators] .- any_cavity_operators
                else
                    subpart = part
                end
                curr_index_subsets_by_op_type[i] = [subpart[1:cum_how_many_each[1]], subpart[cum_how_many_each[1]+1:cum_how_many_each[2]], subpart[cum_how_many_each[2]+1:end]]

            end
            push!(first_indexes, curr_part_indexes)
            push!(index_subsets_by_op_type, curr_index_subsets_by_op_type)
            push!(conjugate_terms, curr_do_conjugate)
        end
    end
    return spin_str, cavity_str, curr_weights, curr_factor, first_indexes, index_subsets_by_op_type, conjugate_terms
end
function prepare_cumulant_to_cumulant_indexed(term_str::String, curr_op_spin_orders::Vector{Int}, op_index_dicts::Dict{String,Op_Group_Type}, cumulant::Full_Cumulant)::Tuple{String,String,Vector{Int},Vector{Vector{Int}},Vector{Vector{Vector{Vector{Int}}}},Vector{Vector{Bool}}}
    how_many_operators = cumulant_order_from_term_str(term_str)
    how_many_spin_operators::Int = sum(curr_op_spin_orders)
    how_many_cavity_operators::Int = length(term_str) - how_many_spin_operators
    any_cavity_operators::Int = 0
    if how_many_cavity_operators > 0
        any_cavity_operators = 1
    end
    #how_many_cavity_operators::Int = how_many_operators - how_many_spin_operators
    partitions::Vector{Vector{Vector{Int}}} = cumulant.partitions
    weights::Vector{Int} = cumulant.weights
    # create X_term_str which replaces the cavity substring (all + and -) with a single X
    X_term_strs::Vector{String} = []
    cavity_str::String = term_str[1:how_many_cavity_operators]
    spin_str::String = term_str[how_many_cavity_operators+1:end]

    X_strs::Vector{String} = []
    curr_X_term_str = ""
    if how_many_cavity_operators > 0
        curr_X_term_str *= "X"
    end
    curr_X_term_str *= spin_str
    push!(X_term_strs, curr_X_term_str)
    push!(X_strs, "")
    curr_partitions = partitions
    curr_weights = weights

    first_indexes::Vector{Vector{Int}} = []
    index_subsets_by_op_type::Vector{Vector{Vector{Vector{Int}}}} = []
    conjugate_terms::Vector{Vector{Bool}} = []
    if length(spin_str) > 0
        for (X_term_str, X_str) in zip(X_term_strs, X_strs)
            X_ind::Int = 0
            # if X in X_term_str, then X_ind is the index of it's position in X_term_str 
            if length(X_term_str) > 0 
                for i in 1:length(X_term_str)
                    if X_term_str[i] == 'X'
                        X_ind = i
                        break
                    end
                end
            end
            for partition in curr_partitions
                curr_part_indexes::Vector{Int} = Vector{Int}(undef, length(partition))
                curr_do_conjugate::Vector{Bool} = Vector{Bool}(undef, length(partition))
                curr_index_subsets_by_op_type::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}(undef, length(partition))
                for (i, part) in enumerate(partition)
                    # if 1 in part, then X needs to replaced with X_str 
                    if maximum(part) > length(X_term_str)
                        error("Error: index larger than length of X_term_str (", X_term_str, ") in part ", part, " of partition ", partition, " of cumulant ", cumulant)
                    end
                    curr_str::String = ""
                    for p in part
                        if p == X_ind 
                            curr_str = cavity_str
                            break
                        end
                    end
                    for p in part
                        if p != X_ind
                            curr_str *= X_term_str[p]
                        else
                            curr_str *= X_str
                        end
                    end
                    curr_do_conjugate[i], new_str = conjugate_to_basis(curr_str)
                    # check if we need to conjugate the string
                    #println(curr_str, " -> ", new_str)
                    curr_op::Op_Group_Type = op_index_dicts[new_str]
                    curr_part_indexes[i] = curr_op.first_index
                    curr_how_many_op_each = curr_op.op_spin_orders
                    cum_how_many_each = cumsum(curr_how_many_op_each)
                    subpart = part[part.>any_cavity_operators] .- any_cavity_operators
                    # only spin operators need indexing subpart is the indexes of part that are larger than how_many_cavity_operators
                    curr_index_subsets_by_op_type[i] = [subpart[1:cum_how_many_each[1]], subpart[cum_how_many_each[1]+1:cum_how_many_each[2]], subpart[cum_how_many_each[2]+1:cum_how_many_each[3]]]

                end
                push!(first_indexes, curr_part_indexes)
                push!(index_subsets_by_op_type, curr_index_subsets_by_op_type)
                push!(conjugate_terms, curr_do_conjugate)
            end
        end
    else
        curr_op = op_index_dicts[term_str]
        curr_first_index = [[curr_op.first_index], [curr_op.first_index]]
        index_subsets_by_op_type = [[[], [], []], [[], [],[]]]
        conjugate_terms = [[false], [false]]
    end
    return spin_str, cavity_str, curr_weights, first_indexes, index_subsets_by_op_type, conjugate_terms
end

function execute_cumulant_to_cumulant_indexed(curr_how_many::Int, spin_str::String, weights::Vector{Int}, factor::Float64, first_indexes::Vector{Vector{Int}}, index_subsets_by_op_type::Vector{Vector{Vector{Vector{Int}}}}, curr_do_conjugate::Vector{Vector{Bool}}, curr_op_spin_orders::Vector{Int}, all_combinations_list::Vector{Vector{Vector{Int}}}, multi_determine_indexes::Function; threaded::Bool=true)::Vector{Cumulant_indexed}
    # Generate all combinations of the first indexes
    cumulants::Vector{Cumulant_indexed} = Vector{Cumulant_indexed}(undef, curr_how_many)
    curr_order, spin_order = string2order_spinorder(spin_str)[1:2]
    do_clamp::Bool = false
    if curr_order == spin_order
        do_clamp = true
    end
    all_comb::Vector{Vector{Int}} = collect(generate_all_combinations_of_combinations(curr_op_spin_orders, all_combinations_list))
    if length(all_comb) != curr_how_many
        error("Error: number of combinations does not match curr_how_many, expected ", curr_how_many, " got ", length(all_comb))
    end
    @usethreads threaded for i in 1:curr_how_many
        comb::Vector{Int} = all_comb[i]
        # Turn comb and curr_cumulant into a Cumulant_indexed using the indexes for the terms in the cumulant
        # use indexes_to_index_vector to turn subsets of comb into index vectors
        curr_partitions::Vector{Vector{Int}} = Vector{Vector{Int}}(undef, length(first_indexes))
        outer_counter::Int = 1
        for (first, index_subset_by_op, do_conjugates) in zip(first_indexes, index_subsets_by_op_type, curr_do_conjugate)
            subpart::Vector{Int} = Vector{Int}(undef, length(index_subset_by_op))
            counter::Int = 1
            for (fi, ind_by_op, conj) in zip(first, index_subset_by_op, do_conjugates)
                index_vector::Vector{Vector{Int}} = indexes_to_index_vector(comb, ind_by_op)
                part::Int = fi + multi_determine_indexes(index_vector) - 1
                if conj
                    subpart[counter] = -part # negative index for conjugate
                else
                    subpart[counter] = part
                end
                counter += 1
            end
            curr_partitions[outer_counter] = subpart
            outer_counter += 1
        end
        # construct cumulant_indexed
        cumulants[i] = Cumulant_indexed(curr_partitions, weights, factor, do_clamp)# new_term_str)
    end
    return cumulants
end

function cumulant_to_cumulant_indexed(curr_how_many::Int, term_str::String, curr_op_spin_orders::Vector{Int}, op_index_dicts::Dict{String,Op_Group_Type}, cumulant::Cumulant, lower_cumulant::Cumulant, all_combinations_list::Vector{Vector{Vector{Int}}}, multi_determine_indexes::Function, max_spins::Int; threaded::Bool=true)::Vector{Cumulant_indexed}
    # Turn comb and curr_cumulant into a Cumulant_indexed using the indexes for the terms in the cumulant
    # First, find the first indexes of the suboperators in term_str 
    spin_str, cavity_str, weights, factor, first_indexes, index_subsets_by_op_type, curr_do_conjugate = prepare_cumulant_to_cumulant_indexed(term_str, curr_op_spin_orders, op_index_dicts, cumulant, lower_cumulant, max_spins)
    return execute_cumulant_to_cumulant_indexed(curr_how_many, spin_str, weights, factor, first_indexes, index_subsets_by_op_type, curr_do_conjugate, curr_op_spin_orders, all_combinations_list, multi_determine_indexes, threaded=threaded)
end

function cumulant_terms_from_dict(cumulant_terms_dict::Dict{String,Op_Group_Type}, op_index_dicts::Dict{String,Op_Group_Type}, multi_determine_indexes::Function, all_combinations_list::Vector{Vector{Vector{Int}}}, max_spins::Int; printing::Bool=false, threaded::Bool=true)::Vector{Cumulant_indexed}
    # construct the cumulant expansions to calculate cumulants from other operators
    # Using Op_Group_Type to calculate the precise indexing 
    # Extract cumulant terms from cumulant_terms_dict
    # if loadsave, generate hash code from inputs (from samples take sample_locations). If hash code exists, load from file. If not, save to file
    first_index::Vector{Int} = []
    how_many::Vector{Int} = []
    term_str::Vector{String} = []
    op_spin_orders::Vector{Vector{Int}} = []
    for key in keys(cumulant_terms_dict)
        cumulant_term = cumulant_terms_dict[key]
        push!(first_index, cumulant_term.first_index)
        push!(how_many, cumulant_term.how_many)
        push!(op_spin_orders, cumulant_term.op_spin_orders)
        push!(term_str, key)
    end
    # sort by first_index and sort how_many accordingly (store the order of the sorting)
    sort_order = sortperm(first_index)
    first_index = first_index[sort_order]
    how_many = how_many[sort_order]
    term_str = term_str[sort_order]
    op_spin_orders = op_spin_orders[sort_order]
    # construct cumulant terms
    # lengths of cumulant terms 
    max_cumulant_order::Int = maximum(length.(term_str))
    max_spin_order::Int = maximum([sum(op_spin_order) for op_spin_order in op_spin_orders])
    #println("cumulant_orders: ", cumulant_orders)

    cumulants_by_order::Vector{Cumulant} = [n_th_order_cumulant(new_cumulant_order, include_n=false) for new_cumulant_order in 1:max_cumulant_order]

    # construct cumulant terms
    how_many_cumulants = sum(how_many)
    how_many_printing = floor(Int, how_many_cumulants / 10^6)
    cumulant_terms::Vector{Cumulant_indexed} = Vector{Cumulant_indexed}(undef, how_many_cumulants)
    for i in progress(1:length(term_str), how_many_printing, text="   -> ", printing=printing, alternative_counter=how_many, min_for_counting=10^6)
        curr_term_str::String = term_str[i]
        curr_first_index::Int = first_index[i]
        curr_how_many::Int = how_many[i]
        curr_op_spin_orders::Vector{Int} = op_spin_orders[i]
        curr_cumulant_order::Int = cumulant_order_from_term_str(curr_term_str)

        curr_cumulant::Cumulant = cumulants_by_order[curr_cumulant_order]
        curr_lower_cumulant::Cumulant = cumulants_by_order[max(1, curr_cumulant_order-1)]
        # Generate index combinations for cumulant
        new_cumulants = cumulant_to_cumulant_indexed(curr_how_many, curr_term_str, curr_op_spin_orders, op_index_dicts, curr_cumulant, curr_lower_cumulant, all_combinations_list, multi_determine_indexes, max_spins, threaded=threaded)
        if !(length(new_cumulants) == curr_how_many)
            error("Error: number of cumulants does not match how_many, expected ", curr_how_many, " for ", curr_term_str, " got ", length(new_cumulants))
        end
        cumulant_terms[curr_first_index:curr_first_index+curr_how_many-1] = new_cumulants
    end
    return cumulant_terms
end
## Test 
#number_of_samples = length(locations)
#cumulant_terms = cumulant_terms_from_dict(number_of_samples, cumulant_terms_dict, operator_strings_dict)

###### This Script is for validation of indexed cumulants, to undo the indexing ################################################################
# Generalize to cumulant_indexed -> turn into Cumulant
function cumulant_indexed_to_string(cumulant_indexed::Cumulant_indexed, operator_strings_dict::Dict{String,Op_Group_Type}, sample_num::Int)
    # Transform first into an indexed op then into a string and then into a term
    partitions_indexed = cumulant_indexed.partitions
    weights = cumulant_indexed.weights
    operator_str = cumulant_indexed.operator
    str = term2str(make_term(operator_str; optim=false); do_braket=true)
    str *= " = "
    for (i, (part, weight)) in enumerate(zip(partitions_indexed, weights))
        if i > 1 && weight > 0
            str *= " + "
        end
        if weight != 1
            str *= string(weight) * " "
        elseif weight == -1
            str *= "-"
        end
        for (j, ind) in enumerate(part)
            if j > 1
                str *= " "
            end
            new_term = index_to_term(ind, operator_strings_dict, sample_num)
            braket_new_term = term2str(new_term; do_braket=true)
            str *= braket_new_term
        end
    end
    return latexstring(str)
end

###################################################################################################
##### Lower Order Full Cumulants ##################################################################
###################################################################################################

# Possible generalisation: use random components instead of the first curr_how_many
function full_execute_cumulant_to_cumulant_indexed(curr_how_many::Int, spin_str::String, weights::Vector{Int}, first_indexes::Vector{Vector{Int}}, index_subsets_by_op_type::Vector{Vector{Vector{Vector{Int}}}}, curr_do_conjugate::Vector{Vector{Bool}}, curr_op_spin_orders::Vector{Int}, all_combinations_list::Vector{Vector{Vector{Int}}}, multi_determine_indexes::Function; threaded::Bool=true)::Vector{Full_Cumulant_indexed}
    # Generate curr_how_many combinations of indexes, create all combinations for how_many=-1
    cumulants::Vector{Full_Cumulant_indexed} = Vector{Full_Cumulant_indexed}(undef, curr_how_many)
    curr_order, spin_order = string2order_spinorder(spin_str)[1:2]
    do_clamp = false
    if curr_order == spin_order
        do_clamp = true
    end
    all_comb = collect(generate_all_combinations_of_combinations(curr_op_spin_orders, all_combinations_list))
    if length(all_comb) != curr_how_many
        error("Error: number of combinations does not match curr_how_many, expected ", curr_how_many, " got ", length(all_comb))
    end
    @usethreads threaded for i in 1:curr_how_many
        comb::Vector{Int} = all_comb[i]
        # Turn comb and curr_cumulant into a Cumulant_indexed using the indexes for the terms in the cumulant
        # use indexes_to_index_vector to turn subsets of comb into index vectors        
        curr_partitions::Vector{Vector{Int}} = Vector{Vector{Int}}(undef, length(first_indexes))
        outer_counter::Int = 1
        for (first, index_subset_by_op, do_conjugates) in zip(first_indexes, index_subsets_by_op_type, curr_do_conjugate)
            subpart::Vector{Int} = Vector{Int}(undef, length(index_subset_by_op))
            counter::Int = 1
            for (fi, ind_by_op, conj) in zip(first, index_subset_by_op, do_conjugates)
                index_vector::Vector{Vector{Int}} = indexes_to_index_vector(comb, ind_by_op)
                part::Int = fi + multi_determine_indexes(index_vector) - 1
                if conj
                    subpart[counter] = -part # negative index for conjugate
                else
                    subpart[counter] = part
                end
                counter += 1
            end
            curr_partitions[outer_counter] = subpart
            outer_counter += 1
        end
        # construct cumulant_indexed
        cumulants[i] = Full_Cumulant_indexed(curr_partitions, weights, do_clamp)
    end
    return cumulants
end

# A modified cumulant_to_cumulant_indexed that calculates Full_Cumulants instead of Cumulants, and doesn't necessarily track all index combinations
function full_cumulant_to_cumulant_indexed(curr_how_many::Int, term_str::String, curr_op_spin_orders::Vector{Int}, op_index_dicts::Dict{String,Op_Group_Type}, cumulant::Full_Cumulant, all_combinations_list::Vector{Vector{Vector{Int}}}, multi_determine_indexes::Function; threaded::Bool=true)::Vector{Full_Cumulant_indexed}
    # Turn comb and curr_cumulant into a Cumulant_indexed using the indexes for the terms in the cumulant
    # First, find the first indexes of the suboperators in term_str 
    spin_str, cavity_str, weights, first_indexes, index_subsets_by_op_type, curr_do_conjugate = prepare_cumulant_to_cumulant_indexed(term_str, curr_op_spin_orders, op_index_dicts, cumulant)
    return full_execute_cumulant_to_cumulant_indexed(curr_how_many, spin_str, weights, first_indexes, index_subsets_by_op_type, curr_do_conjugate, curr_op_spin_orders, all_combinations_list, multi_determine_indexes, threaded=threaded)
end
function lower_order_full_cumulant_terms(operator_index_dict::Dict{String,Op_Group_Type}, max_order::Int, multi_determine_indexes::Function, all_combinations_list::Vector{Vector{Vector{Int}}}, max_spins::Int; threaded::Bool=true)::Vector{Full_Cumulant_indexed}
    #for each operator, how_many, specifies the (maximum) number of index combinations that are evaluated,
    # the ordering procedure from lowest to increasing indexes should esure that we can construct child cumulants and compare orders better than via random choice 
    # if how_many is -1, generate for all possibilities
    # generate cumulant-vector 

    first_index::Vector{Int} = []
    how_many::Vector{Int} = []
    term_str::Vector{String} = []
    op_spin_orders::Vector{Vector{Int}} = []
    cumulant_dict::Dict{String,Op_Group_Type} = Dict{String,Op_Group_Type}()
    for key in keys(operator_index_dict)
        cumulant_term = operator_index_dict[key]
        #curr_order = length(key)
        # removed condition of it being larger than order 1
        push!(first_index, cumulant_term.first_index)
        push!(how_many, cumulant_term.how_many)
        push!(op_spin_orders, cumulant_term.op_spin_orders)
        push!(term_str, key)
    end
    sort_order = sortperm(first_index)
    how_many = how_many[sort_order]
    term_str = term_str[sort_order]
    op_spin_orders = op_spin_orders[sort_order]
    how_many_cumulants = sum(how_many)
    # recalculate first_index from sorted how_many starting from 1 
    first_index = [1]
    for i in 2:length(how_many)
        push!(first_index, first_index[i-1] + how_many[i-1])
    end

    cumulants_by_order::Vector{Full_Cumulant} = [n_th_order_cumulant(new_cumulant_order, include_n=true) for new_cumulant_order in 1:max_order]
    cumulant_vector::Vector{Full_Cumulant_indexed} = Vector{Full_Cumulant_indexed}(undef, how_many_cumulants)

    for i in 1:length(term_str)
        # figure out the order, spin order and cavity order of the operator
        curr_term_str::String = term_str[i]
        curr_first_index::Int = first_index[i]
        curr_how_many::Int = how_many[i]
        curr_op_spin_orders::Vector{Int} = op_spin_orders[i]
        curr_cumulant_order::Int = cumulant_order_from_term_str(curr_term_str)
        #println("curr_term_str: ", curr_term_str, " -> ", curr_cumulant_order)
        # curr_cumulant_order::Int = length(curr_term_str)  # wrong, only  for comparison with literature mistakes
        
        #if curr_cumulant_order > 1
        curr_cumulant::Full_Cumulant = cumulants_by_order[curr_cumulant_order]
        # use cumulant_to_cumulant_indexed to generate cumulants for different index combinations
        cumulants = full_cumulant_to_cumulant_indexed(curr_how_many, curr_term_str, curr_op_spin_orders, operator_index_dict, curr_cumulant, all_combinations_list, multi_determine_indexes, threaded=threaded)
        cumulant_vector[curr_first_index:curr_first_index+curr_how_many-1] = cumulants
        #end
    end
    return cumulant_vector
end