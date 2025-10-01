using Combinatorics
using Base.Threads
include("operator_terms.jl")
include("diff_Eq.jl")
include("indexing.jl")

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
#### Conjugation Scripts ############################################################################################

function is_normally_ordered(term::Term)
    exponents::Vector{Int} = term.exponents
    bosons::String = term.bosons
    if length(exponents) > 2
        error("Error: Not normally ordered (bosons " * bosons * ", exponents " * string(exponents) * ")")
    elseif length(exponents) == 2
        if !(bosons == "+-")
            error("Error: Not normally ordered (bosons " * bosons * ", exponents " * string(exponents) * ")")
        end
    end
end
function is_normally_ordered_return_counts(term::Term)::Tuple{Int,Int}
    exponents::Vector{Int} = term.exponents
    bosons::String = term.bosons
    num_create::Int = 0
    num_annihilate::Int = 0
    if length(exponents) == 1
        if bosons[1] == '+'
            num_create = exponents[1]
        else
            num_annihilate = exponents[1]
        end
    elseif length(exponents) == 2
        if !(bosons == "+-")
            error("Error: Not normally ordered (bosons " * bosons * ", exponents " * string(exponents) * ")")
        end
        num_create = exponents[1]
        num_annihilate = exponents[2]
    elseif length(exponents) > 2
        error("Error: Not normally ordered (bosons " * bosons * ", exponents " * string(exponents) * ")")
    end
    return num_create, num_annihilate
end
function pm_nums(term::Term)
    # determines the number of p and m operators in a term 
    spin_types = term.spin_types
    num_p::Int = 0
    num_m::Int = 0
    for s in spin_types
        if s == 'p'
            num_p += 1
        elseif s == 'm'
            num_m += 1
        end
    end
    return num_p, num_m
end

# No2w the same for strings
function is_normally_ordered(term::String)
    # check if - after a plus 
    had_minus::Bool = false
    for s in term
        if s == '-'
            had_minus = true
        elseif s == '+' && had_minus
            error("Not normally ordered String term: $term")
        end
    end
end
function is_normally_ordered_return_counts(term::String)
    num_create::Int = 0
    num_annihilate::Int = 0
    had_minus::Bool = false
    for s in term
        if s == '-'
            had_minus = true
            num_annihilate += 1
        elseif s == '+'
            num_create += 1
            if had_minus
                error("Not normally ordered String term: $term")
            end
        end
    end
    return num_create, num_annihilate
end
function pm_nums(term::String)
    # determines the number of p and m operators in a term 
    num_p::Int = 0
    num_m::Int = 0
    for s in term
        if s == 'p'
            num_p += 1
        elseif s == 'm'
            num_m += 1
        end
    end
    return num_p, num_m
end
function pmxyz_nums(term::String)::Vector{Int}
    # determines the number of p and m operators in a term 
    nums::Vector{Int} = [0, 0, 0, 0, 0, 0]
    chars::Vector{Char} = ['p', 'm', 'x', 'y', 'z']
    for s in term
        for (i, c) in enumerate(chars)
            if s == c
                nums[i] += 1
                break
            end
        end
    end
    return nums
end
## Test 
#is_normally_ordered_return_counts("++---")    
#pm_nums("ppmm")  

# if an operator has more creation (a^\dagger) than annihilation (a) operators, 
# it is conjugated, and the conjugate flag is set to true
function is_hermitian(term::Term)::Bool
    num_create, num_annihilate = is_normally_ordered_return_counts(term)
    is_hermitian::Bool = num_create == num_annihilate
    is_pauli::Bool = term.is_pauli
    if !is_pauli
        # check hermiticity of spin operators + and - operators present?
        num_p, num_m = pm_nums(term)
        if !(num_p == num_m)
            is_hermitian = false
        end
    end
    return is_hermitian
end
function is_hermitian_output(term::Term)::Tuple{Bool,Int,Int}
    num_create, num_annihilate = is_normally_ordered_return_counts(term)
    is_hermitian::Bool = num_create == num_annihilate
    is_pauli::Bool = term.is_pauli
    num_p::Int, num_m::Int = 0, 0
    if !is_pauli
        # check hermiticity of spin operators + and - operators present?
        num_p, num_m = pm_nums(term)
        if !(num_p == num_m)
            is_hermitian = false
        end
    end
    return is_hermitian, num_create, num_annihilate
end

function is_hermitian(term::String)::Bool
    num_create, num_annihilate = is_normally_ordered_return_counts(term)
    is_hermitian::Bool = num_create == num_annihilate
    # check hermiticity of spin operators + and - operators present?
    num_p, num_m = pm_nums(term)
    if !(num_p == num_m)
        is_hermitian = false
    end
    return is_hermitian
end
function is_hermitian_output(term::String)::Tuple{Bool,Int,Int}
    num_create, num_annihilate = is_normally_ordered_return_counts(term)
    is_hermitian::Bool = num_create == num_annihilate
    # check hermiticity of spin operators + and - operators present?
    num_p, num_m = pm_nums(term)
    if !(num_p == num_m)
        is_hermitian = false
    end
    return is_hermitian, num_create, num_annihilate
end
function is_hermitian_output_all(term::String)::Tuple{Bool, Int, Int, Vector{Int}}
    num_create, num_annihilate = is_normally_ordered_return_counts(term)
    is_hermitian::Bool = num_create == num_annihilate
    # check hermiticity of spin operators + and - operators present?
    nums = pmxyz_nums(term)
    if !(nums[1] == nums[2])
        is_hermitian = false
    end
    return is_hermitian, num_create, num_annihilate, nums
end

function need_to_conjugate(term::Term)::Bool
    num_create, num_annihilate = is_normally_ordered_return_counts(term)
    if term.is_pauli
        if num_create <= num_annihilate
            return false # in the set 
        else
            return true
        end
    else
        # two cases:
        num_p, num_m = pm_nums(term)
        if num_create == num_annihilate # check if  num_p >= num_m
            if num_p >= num_m
                return false
            else
                return true # not in the set 
            end
        else
            if num_create > num_annihilate
                return true
            else
                return false
            end
        end
    end
end
function need_to_conjugate_output(term::Term)::Tuple{Bool,Int,Int}
    num_create, num_annihilate = is_normally_ordered_return_counts(term)
    conjugate::Bool = false
    if term.is_pauli
        if num_create <= num_annihilate
            conjugate = false # in the set 
        else
            conjugate = true
        end
    else
        # two cases:
        num_p, num_m = pm_nums(term)
        if num_create == num_annihilate # check if  num_p >= num_m
            if num_p >= num_m
                conjugate = false
            else
                conjugate = true # not in the set 
            end
        else
            if num_create > num_annihilate
                conjugate = true
            else
                conjugate = false
            end
        end
    end
    return conjugate, num_create, num_annihilate
end

# Generalize need_to_conjugate scripts to strings instead of Terms
function need_to_conjugate(term::String)::Bool
    num_create, num_annihilate = is_normally_ordered_return_counts(term)
    num_p, num_m = pm_nums(term)
    if num_create == num_annihilate # check if  num_p >= num_m
        if num_p >= num_m
            return false
        else
            return true # not in the set 
        end
    else
        if num_create > num_annihilate
            return true
        else
            return false
        end
    end
end
function need_to_conjugate_output(term::String)::Tuple{Bool,Int,Int}
    num_create, num_annihilate = is_normally_ordered_return_counts(term)
    num_p, num_m = pm_nums(term)
    conjugate::Bool = false
    if num_create == num_annihilate # check if  num_p >= num_m
        if num_p >= num_m
            conjugate = false
        else
            conjugate = true # not in the set 
        end
    else
        if num_create > num_annihilate
            conjugate = true
        else
            conjugate = false
        end
    end
    return conjugate, num_create, num_annihilate
end
function need_to_conjugate_output_all(term::String)::Tuple{Bool, Int, Int, Vector{Int}}
    num_create, num_annihilate = is_normally_ordered_return_counts(term)
    nums = pmxyz_nums(term)
    conjugate::Bool = false
    if num_create == num_annihilate # check if  num_p >= num_m
        if nums[1] >= nums[2]
            conjugate = false
        else
            conjugate = true # not in the set 
        end
    else
        if num_create > num_annihilate
            conjugate = true
        else
            conjugate = false
        end
    end
    return conjugate, num_create, num_annihilate, nums
end

function conjugate_to_basis(term::String)::Tuple{Bool, String}
    # use only on normal ordered terms
    # check if no creation operators present after annihilation operators
    # conjugate::Bool = false
    need, num_creation, num_annihilation, nums = need_to_conjugate_output_all(term)
    if need
        # operators::Array{String} = ["p", "m", "x", "y", "z"]
        operators::Array{String} = ["m", "p", "x", "y", "z"]  # switched m and p -> for conjugation
        new_string::String = "+"^num_annihilation * "-"^num_creation 
        for (num, op) in zip(nums, operators)
            new_string *= op^num
        end
        return true, new_string
    end
    return false, term
end
## Test 
#conjugate_to_basis("m")

function conjugate_to_basis!(term::Term)
    # use only on normal ordered terms
    # check if no creation operators present after annihilation operators
    # conjugate::Bool = false
    need, num_creation, num_annihilation = need_to_conjugate_output(term)
    if need
        # conjugate = true
        # make new term
        new_bosons::String = ""
        new_exponents::Vector{Int} = Int[]
        if num_annihilation > 0
            new_bosons *= "+"
            push!(new_exponents, num_annihilation)
        end
        if num_creation > 0
            new_bosons *= "-"
            push!(new_exponents, num_creation)
        end
        term.bosons = new_bosons
        term.exponents = new_exponents
        if !term.is_pauli
            term.spin_types = flip_m_p(term.spin_types)
        end
        # flip conjugate flag
        term.conjugate = !term.conjugate
    end
end
function conjugate_to_basis!(term::Vector{Term})
    for t in term
        conjugate_to_basis!(t)
    end
end
function conjugate_to_basis!(eq::DE_Term)
    # do for all terms in DE_Term
    if need_to_conjugate(eq.exp_op)
        println("Warning: Operator (", term2reducedstr(eq.exp_op), ") should not be in the set of equations")
        conjugate_to_basis!(eq.exp_op)
    end
    # conjugate all operators in the eq
    conjugate_to_basis!(eq.terms)
end
function conjugate_to_basis!(eqs::Vector{DE_Term})
    for eq in eqs
        conjugate_to_basis!(eq)
    end
end
function conjugate_to_basis!(eqs::Vector{Vector{DE_Term}})
    for eq in eqs
        conjugate_to_basis!(eq)
    end
end
function conjugate_to_basis!(eq::DE_Term_Multi)
    if need_to_conjugate(eq.exp_op)
        println("Warning: Operator (", term2reducedstr(eq.exp_op), ") should not be in the set of equations")
        conjugate_to_basis!(eq.exp_op)
    end
    conjugate_to_basis!(eq.sum_i_terms)
    conjugate_to_basis!(eq.sum_i_neq_j_terms)
    conjugate_to_basis!(eq.non_sum_terms)
end
function conjugate_to_basis!(eqs::Vector{DE_Term_Multi})
    for eq in eqs
        conjugate_to_basis!(eq)
    end
end
function conjugate_to_basis!(eqs::Vector{Vector{DE_Term_Multi}})
    for eq in eqs
        conjugate_to_basis!(eq)
    end
end

# Test
#termA = make_term("++xy")
#termB = make_term("+-yx")
#conjugate_to_basis!(termA) # true
#conjugate_to_basis!(termB) # false
#conjugate_to_basis!([termA, termB])
#display(termA)
#display(termB)
#all_eqs = single_spin_gen_DE_up_to_order(2)
#conjugate_to_basis!(all_eqs)
#for all_eq in all_eqs
#    for eq in all_eq
#        display(eq)
#    end
#end

# Add str_conjugate script -- conjugate a str_identifier for a term
function str_conjugate(term_str::String)::String
    # conjugate a term string
    # Example: " "+-xiyj" -> " "+-yixj"
    # for conjugation, only bosons are changed, assums normal ordering of bosons
    # (a^\dagger)^p a^q -> (a^\dagger)^q a^p
    # conjugation only if p > q
    # count +
    num_creation::Int = 0
    num_annihilation::Int = 0
    i::Int = 1
    while i <= length(term_str)
        if term_str[i] == '+'
            num_creation += 1
        else
            break
        end
        i += 1
    end
    while i <= length(term_str)
        if term_str[i] == '-'
            num_annihilation += 1
        else
            break
        end
        i += 1
    end
    if num_creation < num_annihilation
        println("Warning: Term should not be conjugated. Conjugating anyways. ")
    end
    # conjugate
    new_str::String = "+"^num_annihilation * "-"^num_creation
    if i <= length(term_str)
        new_str *= term_str[i:end]
    end
    return new_str
end
## Test
#println(str_conjugate("++-xiyj"))

function separate_index_from_str(str::String, which_ind::Vector{String})::Tuple{String,Int}
    # separate str into pre and post _
    if occursin("_", str)
        pre, post = split(str, "_")
        # find ind via position of post in which_ind
        if post in which_ind
            ind = findfirst(x -> x == post, which_ind)
        else
            ind = 0 # sum ind index (Assumed)
        end
    else
        pre = str
        ind = -1   # not indexed
    end
    return pre, ind
end
## Test

function term_order!(term::Term, which_ind::Vector{String})
    # check if term.spin_indices end on the same elements as which_ind
    spin_indices::String = term.spin_indices
    spin_index_order::Vector{Int} = Vector{Int}(undef, length(spin_indices))
    for (i, ind) in enumerate(spin_indices)
        if string(ind) in which_ind
            spin_index_order[i] = findfirst(x -> x == string(ind), which_ind)
        else
            spin_index_order[i] = 0
        end
    end
    # if spin-index_order is not sorted, sort it, and get indexes of the sorted vector
    if spin_index_order != sort(spin_index_order)
        spin_index_order = sortperm(spin_index_order)
        # reorder spin_indices and spin_types
        term.spin_indices = spin_indices[spin_index_order]
        term.spin_types = term.spin_types[spin_index_order]
    end
end
function term_order!(term::Vector{Term}, which_ind::Vector{String})
    for (i, t) in enumerate(term)
        term_order!(t, which_ind)
        term[i] = t
    end
end
function order_terms_in_DE(eq::DE_Term_Multi)
    # order terms in DE_Term_Multi
    # order terms in sum_i_terms, sum_i_neq_j_terms, non_sum_terms
    # order terms in sum_i_terms
    sum_i_terms = eq.sum_i_terms
    term_order!(sum_i_terms, eq.which_ind)
    eq.sum_i_terms = sum_i_terms
    # order terms in sum_i_neq_j_terms
    sum_i_neq_j_terms = eq.sum_i_neq_j_terms
    term_order!(sum_i_neq_j_terms, eq.which_ind)
    eq.sum_i_neq_j_terms = sum_i_neq_j_terms
    # order terms in non_sum_terms
    non_sum_terms = eq.non_sum_terms
    term_order!(non_sum_terms, eq.which_ind)
    eq.non_sum_terms = non_sum_terms
    return eq
end

function multi_parse_terms_to_indexes!(summed::Int, min_ind::Char, max_order::Int, max_spin_order::Int, how_many_combinations::Vector{Int}, exp_op_str::String, terms::Vector{Term}, constant_terms::Vector{Abstract_Constant_Term}, linear_terms::Vector{Abstract_Term}, cumulant_terms::Vector{Abstract_Term}, vars_dict::Dict{String,Int}, operator_strings_dict::Dict{String,Op_Group_Type}, how_many_cumulant::Int, cumulant_terms_dict::Dict{String,Op_Group_Type})
    how_many_each::Vector{Int} = Vector{Int}()
    how_many_total::Int = 0
    for term in terms
        # figure out if it's a constant term, a cumulant term or a linear term
        term_str, occurances, cavity_occurances = term2reducedstr(term)
        term_spin_order::Int = sum(occurances)
        term_cavity_order = sum(cavity_occurances)
        term_order::Int = term_spin_order + term_cavity_order
        at_max_order::Bool = term_order == max_order
        adder = 0
        if !at_max_order
            adder = 1
        end
        term_order = term_spin_order + term_cavity_order
        # get index of term_str
        if term_order == 0
            do_coeff = true'
            # if constant_terms is a boolean return error
            if isa(constant_terms, Bool)
                error("Found undesired constant term in differential equation of operator (most likely in a summed term, which should not contain constant terms)")
            end
        else
            do_coeff = false
            if term_order <= max_order && term_spin_order <= max_spin_order
                linear = true
                if !haskey(operator_strings_dict, term_str)   # alternative to try, catch
                    error("The operator ( " * term_str * ") was not found in the operator set (term is in differential equation of operator " * exp_op_str * ")")
                end
            else  # cumulant term
                linear = false
                if !haskey(cumulant_terms_dict, term_str)   # alternative to try, catch
                    # add to dict
                    curr_op_str, curr_occurances, _ = term2reducedstr(term)
                    how_many_each = [how_many_combinations[occurance+1] for occurance in curr_occurances]
                    how_many_total = prod(how_many_each)
                    cumulant_terms_dict[term_str] = Op_Group_Type(how_many_cumulant + 1, curr_occurances, how_many_total)
                    how_many_cumulant += how_many_total
                end
            end
        end
        # get coefficient
        if !do_coeff
            if linear
                # add to curr_linear_terms
                push!(linear_terms, term_to_abstract_term(term, vars_dict, summed, min_ind))
            else
                # add to curr_cumulant_terms
                push!(cumulant_terms, term_to_abstract_term(term, vars_dict, summed, min_ind))
            end
        else
            # add to curr_coefficient_terms
            push!(constant_terms, term_to_abstract_constant_term(term, vars_dict, summed, min_ind))
        end
    end
    return how_many_cumulant
end

# for every linear term (cumulant_term, constant term), create the index combination and the term 
function abstract_to_index(curr_terms::Vector{Abstract_Term}, comb::Vector{Int}, how_many_samples::Int, first_indexes::Vector{Int}, multi_determine_indexes::Function, determine_combined_indexes_with_zero::Function, sample_vars::Dict{String,Vector{ComplexF64}}, sample_weights::Dict{String,Vector{Float64}})::Vector{Indexed_Term}
    how_many_summations::Int = 0
    for term in curr_terms
        summed = term.summed
        if !(summed == 0)
            how_many_summations += 1
        end
    end
    total_length::Int = length(curr_terms) + how_many_summations * (how_many_samples - 1)
    curr_terms_indexed::Vector{Indexed_Term} = Vector{Indexed_Term}(undef, total_length)
    curr_pos::Int = 1
    mult_weight::Float64 = 1.0
    for (j, term) in enumerate(curr_terms)
        curr_first_index = first_indexes[j]
        curr_spin_order = term.spin_indices
        do_summation = term.summed # 0 = no sum, 1 = sum_i, 2 = sum_i_neq_j
        curr_sample_vars = term.sample_vars
        curr_coeff = term.coeff
        if do_summation == 0
            inserted_indexes = indexes_to_index_vector(comb, curr_spin_order)
            curr_index = multi_determine_indexes(inserted_indexes) + curr_first_index - 1
            # check if sample vars are present 
            if length(curr_sample_vars) > 0
                # replace sample vars with indexed version of sample_vars dict
                for (var_str, index) in curr_sample_vars
                    curr_coeff *= sample_vars[var_str][comb[index]]
                end
            end
            # create the term
            new_term = Indexed_Term(curr_index, term.conjugate, term.vars, curr_coeff)
            curr_terms_indexed[curr_pos] = new_term
            curr_pos += 1
        else # do summation!
            inserted_indexes, where_zero = indexes_to_index_vectors_with_zero(comb, curr_spin_order)
            curr_indexes::Vector{Int} = determine_combined_indexes_with_zero(inserted_indexes, where_zero)
            for (zero_index, curr_index) in enumerate(curr_indexes)  # the zero index is the index that is summed over
                # check if sample vars are present 
                curr_coeff = term.coeff
                if length(curr_sample_vars) > 1
                    error("Not implemented: more than one sample var in a summed term. Weighting not derived for this case.")
                elseif length(curr_sample_vars) == 1
                    # replace sample vars with indexed version of sample_vars dict 
                    var_str, index = curr_sample_vars[1]
                    if index > 0
                        error("Not implemented: sample variable in a summed term that does not depend on the first term")
                    else
                        mult_weight = sample_weights[var_str][zero_index]
                    end
                    if do_summation == 2  # sum_i_neq_j
                        # check if zero index is present in inserted_indexes[where_zero] -> if yes, subtract from the weight
                        collisions = count(x -> x == zero_index, inserted_indexes[where_zero])
                        if collisions > 0
                            mult_weight -= sample_vars[var_str][zero_index] * collisions
                        end
                    end
                    curr_coeff *= mult_weight
                end
                # create the term
                new_term = Indexed_Term(curr_index + curr_first_index - 1, term.conjugate, term.vars, curr_coeff)
                curr_terms_indexed[curr_pos] = new_term
                curr_pos += 1
            end
        end
    end
    return curr_terms_indexed
end

function abstract_to_index(curr_terms::Vector{Abstract_Constant_Term}, comb::Vector{Int}, how_many_samples::Int, sample_vars::Dict{String,Vector{ComplexF64}})::Vector{Constant_Term}
    how_many_summations::Int = 0
    for term in curr_terms
        summed = term.summed
        if !(summed == 0)
            how_many_summations += 1
        end
    end
    total_length::Int = length(curr_terms) + how_many_summations * (how_many_samples - 1)
    curr_terms_indexed::Vector{Constant_Term} = Vector{Constant_Term}(undef, total_length)
    curr_pos::Int = 1
    for (j, term) in enumerate(curr_terms)
        do_summation = term.summed # 0 = no sum, 1 = sum_i, 2 = sum_i_neq_j
        curr_sample_vars = term.sample_vars
        curr_coeff = term.coefficient
        if do_summation == 0
            # check if sample vars are present 
            if length(curr_sample_vars) > 0
                # replace sample vars with indexed version of sample_vars dict 
                for (var_str, index) in curr_sample_vars
                    curr_coeff *= sample_vars[var_str][comb[index]]
                end
            end
            # create the term
            new_term = Constant_Term(term.variables_index, curr_coeff)
            curr_terms_indexed[curr_pos] = new_term
            curr_pos += 1
        else # do summation! -> no distinction between sum_i and sum_i_neq_j (turned out to be theoretically unnecessary (in the cases we tackled here))
            error("Constant term cannot be summed in current version!")
        end
    end
    return curr_terms_indexed
end

# a function that transforms a DE term into indexes and weights
# also creates a dictionary for the coefficients, and the higher order terms, that need to be cumulant expanded
function prepare_parse_DE_to_indexes(all_eqs::Vector{DE_Term_Multi}, sample_weights::Dict{String,Vector{Float64}}, max_order::Int, max_spins::Int)
    # Outputs for every equation a DE_Term_indexed
    # Preprocessing
    # get first key of sample_weights 
    first_key = first(keys(sample_weights))
    how_many_samples = length(sample_weights[first_key])
    #max_order, max_spin_order, max_cavity_order = get_max_orders(all_eqs)
    how_many_total, operator_strings_dict, multi_determine_indexes, determine_combined_indexes_with_zero, all_combinations_list, how_many_combinations = prepare_all_exp_op_vector(all_eqs, how_many_samples)
    how_many_cumulant::Int = 0
    # Create a dictionary of the vars and of the higher order terms (we will call them cumulant terms)
    vars_dict::Dict{String,Int} = Dict{String,Symbol}()
    cumulant_terms_dict::Dict{String,Op_Group_Type} = Dict{String,Op_Group_Type}()
    # For every DE term, we create a a vector of the indexes of terms, the corresponding vars index and coefficients, a similar vector is created for the cumulant terms
    # still only need these terms -> the sum terms and neqsum terms are expanded into Indexed_Term here to simplify everything to a form that is already covered  by the indexed single spin case

    min_ind::Char = 'h'
    linear_terms::Vector{Vector{Abstract_Term}} = Vector{Vector{Abstract_Term}}()
    cumulant_terms::Vector{Vector{Abstract_Term}} = Vector{Vector{Abstract_Term}}()
    constant_terms::Vector{Vector{Abstract_Constant_Term}} = Vector{Vector{Abstract_Constant_Term}}()
    for eq in all_eqs
        min_ind = eq.sum_ind[1]
        curr_linear_terms = Vector{Abstract_Term}()
        curr_cumulant_terms = Vector{Abstract_Term}()
        curr_constant_terms = Vector{Abstract_Constant_Term}()

        exp_op_str, spin_occurances, _ = term2reducedstr(eq.exp_op)
        order::Int = sum(spin_occurances)
        exp_op_types = operator_strings_dict[exp_op_str] # returns Op_Group_Type(starting_index, operator_orders, how_many)

        for (curr_terms, summed) in zip([eq.non_sum_terms, eq.sum_i_terms, eq.sum_i_neq_j_terms], [0, 1, 2])
            how_many_cumulant = multi_parse_terms_to_indexes!(summed, min_ind, max_order, max_spins, how_many_combinations, exp_op_str, curr_terms, curr_constant_terms, curr_linear_terms, curr_cumulant_terms, vars_dict, operator_strings_dict, how_many_cumulant, cumulant_terms_dict)
        end
        push!(linear_terms, curr_linear_terms)
        push!(cumulant_terms, curr_cumulant_terms)
        push!(constant_terms, curr_constant_terms)
    end
    return vars_dict, how_many_cumulant, (how_many_total, operator_strings_dict, multi_determine_indexes, determine_combined_indexes_with_zero, cumulant_terms_dict, linear_terms, cumulant_terms, constant_terms, how_many_combinations, how_many_samples, all_combinations_list)
end
function parse_DE_to_indexes(all_eqs::Vector{DE_Term_Multi}, sample_weights::Dict{String,Vector{Float64}}, sample_vars::Dict{String,Vector{ComplexF64}}, big_tuple; threaded::Bool=true)
    how_many_total::Int = big_tuple[1]
    operator_strings_dict::Dict{String,Op_Group_Type} = big_tuple[2]
    multi_determine_indexes::Function = big_tuple[3]
    determine_combined_indexes_with_zero::Function = big_tuple[4]
    cumulant_terms_dict::Dict{String,Op_Group_Type} = big_tuple[5]
    linear_terms::Vector{Vector{Abstract_Term}} = big_tuple[6]
    cumulant_terms::Vector{Vector{Abstract_Term}} = big_tuple[7]
    constant_terms::Vector{Vector{Abstract_Constant_Term}} = big_tuple[8]
    how_many_combinations::Vector{Int} = big_tuple[9]
    how_many_samples::Int = big_tuple[10]
    all_combinations_list::Vector{Vector{Vector{Int}}} = big_tuple[11]

    all_eqs_indexed::Vector{DE_Term_indexed} = Vector{DE_Term_indexed}(undef, how_many_total)
    for (eq, curr_linear_terms, curr_cumulant_terms, curr_constant_terms) in zip(all_eqs, linear_terms, cumulant_terms, constant_terms)
        exp_op_str, spin_occurances, _ = term2reducedstr(eq.exp_op)
        order::Int = sum(spin_occurances)
        exp_op_types = operator_strings_dict[exp_op_str] # returns Op_Group_Type(starting_index, operator_orders, how_many)

        do_clamped = false # if only spin operators are present -> do_clamped = true
        if order == length(exp_op_str)
            do_clamped = true
        end
        first_index::Int = exp_op_types.first_index
        how_many::Int = exp_op_types.how_many

        op_spin_orders::Vector{Int} = exp_op_types.op_spin_orders
        # vector of the combinations 
        first_indexes_linear_terms::Vector{Int} = Vector{Int}(undef, length(curr_linear_terms))
        first_indexes_cumulant_terms::Vector{Int} = Vector{Int}(undef, length(curr_cumulant_terms))
        for (i, term) in enumerate(curr_linear_terms)
            first_indexes_linear_terms[i] = operator_strings_dict[term.term_str].first_index
        end
        for (i, term) in enumerate(curr_cumulant_terms)
            first_indexes_cumulant_terms[i] = cumulant_terms_dict[term.term_str].first_index
        end
        # create all
        all_comb::Vector{Vector{Int}} = collect(generate_all_combinations_of_combinations(op_spin_orders, all_combinations_list))
        @usethreads threaded for i in 1:how_many
            comb::Vector{Int} = all_comb[i]
            # for every linear_term (cumulant_term, constant_term), create the index combination and the term 
            curr_linear_terms_indexed = abstract_to_index(curr_linear_terms, comb, how_many_samples, first_indexes_linear_terms, multi_determine_indexes, determine_combined_indexes_with_zero, sample_vars, sample_weights)
            curr_cumulant_terms_indexes = abstract_to_index(curr_cumulant_terms, comb, how_many_samples, first_indexes_cumulant_terms, multi_determine_indexes, determine_combined_indexes_with_zero, sample_vars, sample_weights)
            curr_constant_terms_indexed = abstract_to_index(curr_constant_terms, comb, how_many_samples, sample_vars)
            # create the DE_Term_indexed
            new_de_term = DE_Term_indexed(first_index + i - 1, curr_linear_terms_indexed, curr_cumulant_terms_indexes, curr_constant_terms_indexed, do_clamped)
            all_eqs_indexed[first_index+i-1] = new_de_term
        end
    end
    #println(counter)
    return all_eqs_indexed, operator_strings_dict, cumulant_terms_dict
end

#### Invert dictionaries and flatten them ####################################################################

# flatten vector of dictionaries into a single dictionary for any type of dict
function flatten_dict(dict_vect::Dict, top_level::Bool=true)
    # return keys and vals
    key_type = keytype(dict_vect)
    val_type = valtype(dict_vect)
    keys = Vector{key_type}()
    vals = Vector{val_type}()
    for (k, v) in dict_vect
        push!(keys, k)
        push!(vals, v)
    end
    return keys, vals
end

function flatten_dict(dict_vect::Vector, top_level::Bool=true)
    # check type
    # get keys from subtypes
    keys, vals = flatten_dict(dict_vect[1], false)
    for i in 2:length(dict_vect)
        curr_keys, curr_vals = flatten_dict(dict_vect[i], false)
        append!(keys, curr_keys)
        append!(vals, curr_vals)
    end
    if !top_level
        return keys, vals
    else
        # create new dict
        key_type = typeof(keys[1])
        val_type = typeof(vals[1])
        new_dict = Dict{key_type,val_type}()
        for i in 1:length(keys)
            new_dict[keys[i]] = vals[i]
        end
        return new_dict
    end
end
## Test
#dict_vec = [Dict{String, Int}("a"=>1, "b"=>2), Dict{String, Int}("c"=>3, "d"=>4)]
#println(flatten_dict(dict_vec))

function invert_dict(dict)
    # Inverts dictionaries and vectors of (vectors of...) dictionaries val => key
    # requires bijective dictionaries (no duplicate values)
    if isa(dict, Dict)
        key_type = keytype(dict)
        value_type = valtype(dict)
        inverted_dict::Dict{value_type,key_type} = Dict{value_type,key_type}()
        for (k, v) in dict
            if !haskey(inverted_dict, v)
                inverted_dict[v] = k
            else
                error("Dictionary is not bijective ($v is already a key in inverted dict). ")
            end
        end
        return Dict(v => k for (k, v) in dict)
    elseif isa(dict, Vector)
        return [invert_dict(d) for d in dict]
    else
        throw(ArgumentError("Unsupported type"))
    end
end
## Test with a single Dict
#d1 = Dict("a" => 1, "b" => 2)
#println(invert_dict(d1))
## Test with a Vector of Vectors of Dicts
#d2 = [[Dict("a" => 1, "b" => 2), Dict("c" => 3, "d" => 4)], [Dict("e" => 5, "f" => 6)]]
#println(invert_dict(d2))

function index_dict_to_vector(Dict)::Vector
    key_type = keytype(Dict)
    val_type = valtype(Dict)
    # key_type needs to be an integer
    if !(key_type <: Integer) # <: is a subtype operator
        error("Key type is not an integer")
    end
    vec::Vector{val_type} = Vector{val_type}(undef, length(Dict))
    for (k, v) in Dict
        if k > length(Dict)
            error("Key is larger than length of vector")
        end
        vec[k] = v
    end
    return vec
end
##Test
#d1 = Dict(1 => "a", 2 => "b", 3 => "c")
#println(index_dict_to_vector(d1))

function reduce_partitions2weighted_unique(partitions::Vector{Vector{Int}}, weights::Vector{Int})::Tuple{Vector{Vector{Int}},Vector{Int}}#, Vector{Bool}}
    added_indexes = Vector{Int}()
    reduced_partitions = Vector{Vector{Int}}()
    reduced_weights = Vector{Int}()
    curr_weight::Int = 0
    #curr_weight_diff::Int = 0
    #reduced_real::Vector{Bool} = Vector{Bool}()
    for i in 1:length(partitions)
        if !(i in added_indexes)
            # add partition to reduced_partitions
            push!(reduced_partitions, partitions[i])
            # add weight to reduced_weights
            curr_weight = weights[i]
            curr_weight_diff = weights[i]
            # add all partitions that are equal to the current partition to added_indexes
            for j in i+1:length(partitions)
                if !(j in added_indexes)
                    if partitions[i] == partitions[j]
                        push!(added_indexes, j)
                        curr_weight += weights[j]
                        #curr_weight_diff += weights[j]
                        #elseif partitions[i] == -partitions[j]
                        #    push!(added_indexes, j)
                        #    curr_weight += weights[j]
                        #    curr_weight_diff -= weights[j]
                    end
                end
            end
            #if curr_weight_diff != 0
            #    push!(reduced_real, false)
            #else
            #    push!(reduced_real, true)
            #end
            push!(reduced_weights, curr_weight)
        end
    end
    return reduced_partitions, reduced_weights#, reduced_real
end
## Test
#part = [[1, 1, 4], [5, 4], [9, 1], [9, 1]]
#weights = [-2, 1, 1, 1]
#println(reduce_partitions2weighted_unique(part, weights))
#part2 = [[-1, 1, 4], [6, 4], [-9, 1], [9, -1]]
#weights2 = [-2, 1, 1, 1]
#println(reduce_partitions2weighted_unique(part2, weights2))


##### Remove Zero terms from Differential Equation Terms ##########################################################################
#### Remove Empty Terms => Terms with parameters that are zero (for simplified evaluation of DE)
function make_keep_element_vec(terms::Vector{Indexed_Term}, remove_indexes::Vector{Int})
    keep_elements::Vector{Int} = Vector{Int}()
    for (i, term) in enumerate(terms)
        if !(term.variables_index in remove_indexes)
            push!(keep_elements, i)
        end
    end
    new_terms = terms[keep_elements]
    return new_terms
end
function make_keep_element_vec(terms::Vector{Constant_Term}, remove_indexes::Vector{Int})
    keep_elements::Vector{Int} = Vector{Int}()
    for (i, term) in enumerate(terms)
        if !(term.variables_index in remove_indexes)
            push!(keep_elements, i)
        end
    end
    new_terms = terms[keep_elements]
    return new_terms
end


function remove_terms_with_zero_parameters(terms::DE_Term_indexed, remove_indexes::Vector{Int})::DE_Term_indexed
    # go through terms, remove terms that have indexes in remove_indexes
    new_linear_terms = make_keep_element_vec(terms.linear_terms, remove_indexes)
    new_cumulant_terms = make_keep_element_vec(terms.cumulant_terms, remove_indexes)
    new_constant_terms = make_keep_element_vec(terms.constant_terms, remove_indexes)
    return DE_Term_indexed(terms.exp_index, new_linear_terms, new_cumulant_terms, new_constant_terms, terms.clamped)
end

function find_zero_parameters(param_generator::Function, constant_value_indexes::Vector{Int}; how_many_eps::Int=5)::Vector{Int}
    parameters = param_generator(0.0)
    remove_indexes::Vector{Int} = Int[]
    for i in constant_value_indexes
        # if parameter is apprixmately equal to zero then remove it from equations
        # use julia internal approximately equal function
        if abs(real(parameters[i])) < eps(Float64) * how_many_eps && abs(imag(parameters[i])) < eps(Float64) * how_many_eps
            push!(remove_indexes, i)
        end
    end
    return remove_indexes
end
function remove_zero_terms(terms::Vector{DE_Term_indexed}, param_generator::Function, constant_value_indexes::Vector{Int}; how_many_eps::Int=5)::Vector{DE_Term_indexed}
    # create parameters at point 
    remove_indexes = find_zero_parameters(param_generator, constant_value_indexes; how_many_eps=how_many_eps)
    new_terms::Vector{DE_Term_indexed} = Vector{DE_Term_indexed}(undef, length(terms))
    for (i, term) in enumerate(terms)
        new_terms[i] = remove_terms_with_zero_parameters(term, remove_indexes)
    end
    return new_terms
end
#function remove_zero_terms()
## Test
#param_generator, constant_value_indexes = value_generator(pulse_param, val_dict, vars_vec)
#all_eqs_indexed_reduced = remove_zero_terms(all_eqs_indexed, param_generator, constant_value_indexes; how_many_eps=5)