using Combinatorics
using Base.Threads
using Random
include("operator_terms.jl")
include("diff_Eq.jl")
include("combinatorics.jl")
include("indexing_binomial.jl")

# create a list of all operators in set of equations and index them. then create a list of all higher order operators that need to be cumulant expanded
# the operators and indexes need to be matched via dictionaries. each dictionary matches the operator to it's index for operators of a particular order and spin order
function gen_operator_index_dictionaries(all_eqs::Union{Vector{DE_Term},Vector{DE_Term_Multi},Vector{Vector{DE_Term}},Vector{Vector{DE_Term_Multi}}}, max_order::Int=-1, max_spin_order::Int=-1)::Vector{Vector{Dict{String,Int}}}
    # if all_eqs isa Vector{DE_Term} or  isa Vector{DE_Term_Multi}, do not need to sort order of operators
    ### assume last term is highest order (normal ordering)
    if isa(all_eqs, Vector{DE_Term}) || isa(all_eqs, Vector{DE_Term_Multi})
        println("Transforming input into matrix for further processing")
        all_eqs = make_terms_matrix(all_eqs)
    end
    if max_order < 0 || max_spin_order < 0
        _, max_order, max_spin_order = term2str_identifier(all_eqs[end][end].exp_op)
    end
    # Walk through vectors or vectors of vectors and make a list of all operators
    nth_order_term_indexes::Vector{Dict{String,Int}} = Vector{Dict{String,Int}}()
    op_index_dicts::Vector{Vector{Dict{String,Int}}} = Vector{Vector{Dict{String,Int}}}(undef, max_order)  # orders tstart at 1
    curr_max_spin_order::Int = 0
    i::Int = 0
    for curr_order in 1:max_order   # if slow use @threads here
        curr_max_spin_order = min(curr_order, max_spin_order)
        nth_order_term_indexes = [Dict{String,Int}() for i in 1:curr_max_spin_order+1]  # order of operators starts at 1
        # Walk through all equations of nth order
        curr_eqs = all_eqs[curr_order]
        for curr_eq in curr_eqs
            curr_str, order, spin_order = term2str_identifier(curr_eq.exp_op)
            if order != curr_order
                error("Order of operator does not match expected order")
            end
            # add term to dictionary
            i += 1
            nth_order_term_indexes[spin_order+1][curr_str] = i
        end
        op_index_dicts[curr_order] = nth_order_term_indexes
    end
    return op_index_dicts
end
## test
#all_eqs = single_spin_gen_DE_up_to_order(2)
#op_index_dicts = gen_operator_index_dictionaries(all_eqs)

# A set of functions that generalizes this to the multi qubit case: 
function extract_all_exp_op_vector(all_eqs::Vector{DE_Term_Multi})::Tuple{Vector{String},Vector{Vector{Int}},Vector{Vector{Int}}}
    # Create a vector of all operator term types (i.e. +-xxyz) and the number of times each spin type operator appears in the term
    operator_strings::Vector{String} = []
    operator_spin_orders::Vector{Vector{Int}} = []
    operator_cavity_orders::Vector{Vector{Int}} = []
    curr_occurances::Vector{Int} = zeros(Int, 3)
    curr_cavity_occurances::Vector{Int} = zeros(Int, 2)
    for curr_eq in all_eqs
        curr_term = curr_eq.exp_op
        curr_op_str, curr_occurances, curr_cavity_occurances = term2reducedstr(curr_term)
        push!(operator_strings, curr_op_str)
        push!(operator_spin_orders, curr_occurances)
        push!(operator_cavity_orders, curr_cavity_occurances)
    end
    return operator_strings, operator_spin_orders, operator_cavity_orders
end

function get_max_orders(eqs::Vector{DE_Term_Multi})::Tuple{Int,Int,Int}
    # get the maximum order of operators and maximum spin order of operators in the equations
    max_order::Int = 0
    max_spin_order::Int = 0
    max_cavity_order::Int = 0
    for eq in eqs
        _, order, spin_order = term2str_identifier(eq.exp_op)
        max_order = max(max_order, order)
        max_spin_order = max(max_spin_order, spin_order)
        max_cavity_order = max(max_cavity_order, order - spin_order)
    end
    return max_order, max_spin_order, max_cavity_order
end
mutable struct Indexed_Op
    op_vec::Vector{Int} # how many +, - , x, y, z 
    op_indexes::Vector{Vector{Int}} # indexes of operators x, y, z
    function Indexed_Op(op_vec::Vector{Int}, op_indexes::Vector{Vector{Int}})
        new(op_vec, op_indexes)
    end
end
mutable struct Op_Group_Type
    first_index::Int
    op_spin_orders::Vector{Int}
    how_many::Int
    function Op_Group_Type(first_index::Int, op_spin_orders::Vector{Int}, how_many::Int)
        new(first_index, op_spin_orders, how_many)
    end
end
function prepare_all_exp_op_vector(all_eqs::Vector{DE_Term_Multi}, number_of_samples::Int, how_many_orders_more::Int=1)
    # genreate the expectation value vector and lists about the content from the all_eqs
    operator_strings, operator_spin_orders, _ = extract_all_exp_op_vector(all_eqs)
    # find maximum cavity order number of plus and minus operators
    # first create vectors of index combinations
    max_order, max_spin_order, max_cavity_order = get_max_orders(all_eqs)
    multi_determine_indexes, determine_combined_indexes_with_zero = determine_multi_indexes_gen(number_of_samples, max_spin_order + how_many_orders_more)
    all_combinations_list::Vector{Vector{Vector{Int}}} = [geq_int_combinations(number_of_samples, i) for i in 1:max_spin_order+how_many_orders_more]
    all_combinations_lengths::Vector{Int} = ones(Int, length(all_combinations_list) + 1)
    for i in 1:length(all_combinations_list)
        all_combinations_lengths[i+1] = length(all_combinations_list[i])
    end
    operator_strings_numbers::Vector{Int} = zeros(Int, length(operator_strings))
    operator_how_many::Vector{Int} = zeros(Int, length(operator_strings))
    curr_index::Int = 1
    how_many_each::Vector{Int} = zeros(Int, 3)
    how_many_total::Int = 0
    for (i, op_order) in enumerate(operator_spin_orders)
        # save the first index for the operator string
        how_many_each = [all_combinations_lengths[occurance+1] for occurance in op_order]
        how_many_total = prod(how_many_each)
        operator_strings_numbers[i] = curr_index
        operator_how_many[i] = how_many_total
        curr_index += how_many_total
    end
    how_many_total = curr_index - 1
    # now create the expectation value vector
    # Create a dictionary with the operator strings as keys and , (first_index, op_spin_orders) as values
    operator_strings_dict::Dict{String,Op_Group_Type} = Dict{String,Op_Group_Type}()
    for (i, (op_inds, op_str, op_order)) in enumerate(zip(operator_strings_numbers, operator_strings, operator_spin_orders))
        operator_strings_dict[op_str] = Op_Group_Type(op_inds, op_order, operator_how_many[i])
    end
    return how_many_total, operator_strings_dict, multi_determine_indexes, determine_combined_indexes_with_zero, all_combinations_list, all_combinations_lengths
end

function op_group_type_dict_to_vector(dict::Dict{String,Op_Group_Type}, n_samples::Int)::Vector{Indexed_Op}
    # get first indexes
    first_indexes::Vector{Int} = Vector{Int}()
    strings::Vector{String} = Vector{String}()
    op_spin_orders::Vector{Vector{Int}} = Vector{Vector{Int}}()
    how_many::Vector{Int} = Vector{Int}()
    for key in keys(dict)
        vals = dict[key]
        push!(first_indexes, vals.first_index)
        push!(strings, key)
        push!(op_spin_orders, vals.op_spin_orders)
        push!(how_many, vals.how_many)
    end
    # sort by first index
    order = sortperm(first_indexes)
    first_indexes = first_indexes[order]
    strings = strings[order]
    op_spin_orders = op_spin_orders[order]
    how_many = how_many[order]
    how_many_total = first_indexes[end] + how_many[end] - 1
    # Generate a vector with Strings and indexes for the operators
    new_vec::Vector{Indexed_Op} = Vector{String}(undef, how_many_total)
    for (curr_str, first_ind, spin_ord, n) in zip(strings, first_indexes, op_spin_orders, how_many)
        # iterate through spin_ord and generate all combinations of operators
        op_vec::Vector{Int} = zeros(Int, 5)
        op_vec[1] = count(==('+'), curr_str)
        op_vec[2] = count(==('-'), curr_str)
        op_vec[3:5] = spin_ord
        i::Int = 0
        for curr_inds in multi_mgreater_range(n_samples, spin_ord)
            # generate the operator vector
            i += 1
            new_vec[first_ind+i-1] = Indexed_Op(op_vec, curr_inds)
        end
        if !(i == n)
            println("Error in op_group_type_dict_to_vector. Expected to construct $n operators, but constructed $i operators")
        end
    end
    return new_vec
end
# In order to find indexes corresponding to a operator
function preprocess_index_inversion(vars_dict::Dict{String,Int})::Vector{Vector{String}}

    inverse_vars_vec_pre::Vector{String} = index_dict_to_vector(invert_dict(vars_dict))
    inverse_vars_vec::Vector{Vector{String}} = Vector{Vector{String}}()
    for elem in inverse_vars_vec_pre
        #println(elem)
        # split the strings at ? to create vector of strings
        push!(inverse_vars_vec, split(elem, "?"))
    end
    return inverse_vars_vec
end
## Test
#all_eqs = single_spin_gen_DE_up_to_order(2)
#conjugate_to_basis!(all_eqs)
#all_eqs_indexed, op_index_dicts, vars_dict, cumulant_terms_dict = parse_DE_to_indexes(all_eqs)
#inverse_op_index_vec, inverse_vars_vec, inverse_cumulant_terms_vec = preprocess_index_inversion(op_index_dicts, vars_dict, cumulant_terms_dict)
#display(inverse_op_index_vec)
#display(inverse_vars_vec)
#display(inverse_cumulant_terms_vec)

function Random_Integers(min::Int=1, max::Int=10, how_many::Int=1; rng=MersenneTwister(1234))::Vector{Int}
    # if rng exists as a variable in the global scope, use it, otherwise create a new one
    return rand(rng, min:max, how_many)
end
## Test
#Random_Integers(1,10,10, printing=true)



#### Indexing Operators #########################################################################
@inline function indexes_to_index_vector(indexes::Vector{Int}, index_order::Vector{Int})::Vector{Vector{Int}}
    # use the index_order to construct a grouped indexing vector
    # construct a vector of same content vectors as index_order 
    inds::Vector{Int} = cumsum([1, index_order...])
    index_vector::Vector{Vector{Int}} = [indexes[inds[i]:inds[i+1]-1] for i in 1:length(index_order)]
    return index_vector
end
@inline function indexes_to_index_vector(indexes::Vector{Int}, index_order::Vector{Vector{Int}})::Vector{Vector{Int}}
    # use the index_order to construct a grouped indexing vector
    # construct a vector of same content vectors as index_order 
    index_vector::Vector{Vector{Int}} = Vector{Vector{Int}}(undef, length(index_order))
    curr_indexes::Vector{Int} = Vector{Int}(undef, length(index_order))
    for (i, ind) in enumerate(index_order)
        curr_indexes = sort!(indexes[ind])
        index_vector[i] = curr_indexes
    end
    return index_vector
end
@inline function indexes_to_index_vectors_with_zero(indexes::Vector{Int}, index_order::Vector{Vector{Int}})::Tuple{Vector{Vector{Int}},Int}
    # use the index_order to construct a grouped indexing vector
    # construct a vector of same content vectors as index_order 
    index_vector::Vector{Vector{Int}} = Vector{Vector{Int}}(undef, length(index_order))
    where_zero::Int = -1
    curr_indexes::Vector{Int} = Vector{Int}[]
    for (i, ind) in enumerate(index_order)
        if (0 in ind)
            where_zero = i
            # remove zero from ind 
            sort!(ind)
            if length(ind) > 1
                curr_indexes = sort!(indexes[ind[2:end]])
            else
                curr_indexes = []
            end
        else
            curr_indexes = sort!(indexes[ind])
        end
        index_vector[i] = curr_indexes
    end
    return index_vector, where_zero
end
## Test 
#indexes = [2,3,5]
#index_order = Vector{Int}[[1,2], [], [0,2]]
#indexes, where_zero = indexes_to_index_vectors_with_zero(indexes, index_order)


#### Get the indexes from an indexed equations ###################################################
# This function is meant to debug the code of this package, to see if the equations are transformed correctly
# a function that computes the operator and index from just the index 
function index_to_indexed_op(index::Int, operator_strings_dict::Dict{String,Op_Group_Type}, sample_num::Int)
    # not fast but easy to use
    # go through the keys and find the right one [interval from .first_index to .first_index+.how_many]
    op_found::String = ""
    for op in keys(operator_strings_dict)
        start = operator_strings_dict[op].first_index
        stop = operator_strings_dict[op].first_index + operator_strings_dict[op].how_many - 1
        if index >= start && index <= stop
            op_found = op
            break
        end
    end
    if length(op_found) == 0
        error("Operator not found")
    end
    # now we know the operator, find the index
    op = operator_strings_dict[op_found]
    op_spin_order = op.op_spin_orders # vector with orders of spins
    # iterate through indexes until we find the index
    curr_index = op.first_index - 1

    for i_s in multi_mgreater_range(sample_num, op_spin_order)
        curr_index += 1
        if curr_index == index
            return op_found, i_s
        end
    end
    error("Index not found")
end

function index_to_term(index::Int, operator_strings_dict::Dict{String,Op_Group_Type}, sample_num::Int)
    # Transform first into an indexed op then into a string and then into a term
    op_str, op_index = index_to_indexed_op(index, operator_strings_dict, sample_num)
    # make string "+--x1y2z3" and so on
    str::String = ""
    index_count::Vector{Int} = ones(Int, 3)
    for elem in op_str
        if elem == '+'
            str *= "+"
        elseif elem == '-'
            str *= "-"
        elseif elem == 'x'
            str *= "x" * string(op_index[1][index_count[1]])
            index_count[1] += 1
        elseif elem == 'y'
            str *= "y" * string(op_index[2][index_count[2]])
            index_count[2] += 1
        elseif elem == 'z'
            str *= "z" * string(op_index[3][index_count[3]])
            index_count[3] += 1
        end
    end
    return make_term(str, optim=false)
end
## Check if indexing is correct
#display(all_eqs[1])
## find the operators for all terms in the term
#eq_ind = all_eqs_indexed[1]
#index_to_term(eq_ind.exp_index, operator_strings_dict, sample_num)

# Generalize index_to_Term to Indexed_Term and DE_Term_indexed
function index_to_term(indexed_term::Indexed_Term, operator_strings_dict::Dict{String,Op_Group_Type}, vars_vec::Vector{Vector{String}}, sample_num::Int)::Term
    # Transform first into an indexed op then into a string and then into a term
    operator_index = indexed_term.operator_index
    conjugate = indexed_term.conjugate
    variables_index = indexed_term.variables_index
    coefficients = indexed_term.coefficient
    op_str, op_index = index_to_indexed_op(operator_index, operator_strings_dict, sample_num)
    spin_str::String = "x"^length(op_index[1]) * "y"^length(op_index[2]) * "z"^length(op_index[3])
    spin_indices::String = ""
    for op in op_index[1:3]
        for i in op
            spin_indices *= string(i)
        end
    end
    bosons_str::String = ""
    exponents::Vector{Int} = []
    # how many + and - are there
    num_p::Int = count(letter -> letter == '+', op_str)
    num_m::Int = count(letter -> letter == '-', op_str)
    if num_p > 0
        bosons_str *= "+"
        push!(exponents, num_p)
    end
    if num_m > 0
        bosons_str *= "-"
        push!(exponents, num_m)
    end
    vars = vars_vec[variables_index]
    return Term(spin_types=spin_str, spin_indices=spin_indices, bosons=bosons_str, exponents=exponents, coeff=coefficients, vars=vars, conjugate=conjugate)
end
# Generalize to Constant_Term
function index_to_term(constant_term::Constant_Term, vars_vec::Vector{Vector{String}})::Term
    return Term(spin_types="", spin_indices="", bosons="", exponents=Vector{Int}(), coeff=constant_term.coefficient, vars=vars_vec[constant_term.variables_index], conjugate=false)
end
# Generalize to DE_Term_indexed
function index_to_term(de_term_indexed::DE_Term_indexed, operator_strings_dict::Dict{String,Op_Group_Type}, cumulant_terms_dict::Dict{String,Op_Group_Type}, vars_vec::Vector{Vector{String}}, sample_num::Int)::DE_Term
    exp_index::Term = index_to_term(de_term_indexed.exp_index, operator_strings_dict, sample_num)
    which_ind::Vector{String} = []
    terms::Vector{Term} = []
    lin_terms::Vector{Indexed_Term} = de_term_indexed.linear_terms
    cum_terms::Vector{Indexed_Term} = de_term_indexed.cumulant_terms
    const_terms::Vector{Constant_Term} = de_term_indexed.constant_terms
    for ind_term in lin_terms
        push!(terms, index_to_term(ind_term, operator_strings_dict, vars_vec, sample_num))
    end
    for ind_term in const_terms
        push!(terms, index_to_term(ind_term, vars_vec))
    end
    for ind_term in cum_terms
        push!(terms, index_to_term(ind_term, cumulant_terms_dict, vars_vec, sample_num))
    end
    return DE_Term(exp_index, which_ind, terms=terms)
end
# Test 
#display(all_eqs[1])
#display(index_to_term(all_eqs_indexed[1], operator_strings_dict, cumulant_terms_dict, vars_vec, sample_num))
