module qExpressions
using ..qSpace
using ..StringUtils
using ComplexRationals
import Base: show, adjoint, iterate, length, eltype, +, -, sort, *, ^, product, iszero, copy

export qExpr, qTerm, qEQ, qSum, diff_qEQ, term, simplify, base_operators, Sum, ∑, flatten, neq, d_dt

""" 
    qExpr

The abstract type `qExpr` is the base type for all quantum expressions in this module.
"""
abstract type qExpr end

Is = Union{Int,Vector{Int}}

"""
    qTerm

A `qTerm` represents a single term in a quantum expression. It contains:
    - `coeff`: The coefficient of the term, which can be a number (e.g., `Int`, `Float64`, `Rational`, etc.).
    - `var_exponents`: A vector of integers representing the exponents of the state variables, defined in a StateSpace.
    - `op_indices`: A vector of indices representing the operators in the term, which are also defined in a StateSpace.
"""
mutable struct qTerm <: qExpr
    coeff::Number                   # Coefficient (supports Complex, Rational, Int, etc.)
    var_exponents::Vector{Int}      # Exponents for each state variable (order matching StateSpace.variables).
    op_indices::Vector{Is}
end

"""
    qEQ

A `qEQ` represents a quantum equation, consisting of a Vector of quantum Expressions representing the additive terms of the equation.
It also contains a reference to the state space in which the equation is defined.
"""
mutable struct qEQ
    terms::Vector{qExpr}         # Vector of terms
    statespace::StateSpace
end

""" 
    qSum

A `qSum` represents the summation of a quantum Equation over indexes in a quantum expression.
It contains:
    - `expr`: The expression being summed over, which is a `qEQ` object.
    - `indexes`: A vector of strings representing the summation indexes (e.g., "i").
    - `subsystem_index`: The index of the subspace in which the indexes live. 
    - `element_indexes`: A vector of integers representing the position of the indexes in that subspace.
    - `neq`: A boolean indicating whether different indexes in the sum can refer to the same element in the subspace. 
            For example, the indexes i,j,k can refer to different elements in a much larger bath of elements. 
"""
mutable struct qSum <: qExpr
    expr::qEQ       # The expression being summed over.
    indexes::Vector{String}   # The summation index (e.g. "i").
    subsystem_index::Int  # The subspace index where the summation index was found.
    element_indexes::Vector{Int}    # The position in that subspace.
    neq::Bool
end

raw""" 
    diff_qEQ

A `diff_qEQ` represents a differential equation d/dt <LHS> = <RHS>, where LHS is the expectation value of a qTerm and similarly the RHS is the expectation value of a qEQ.
It contains: 
    - `left_hand_side`: The left-hand side of the equation, which is a `qTerm` object.
    - `right_hand_side`: The right-hand side of the equation, which is a `qEQ` object.
    - `statespace`: The state space in which the equation is defined.
    - `braket`: (default=true) A boolean indicating whether to use braket notation for the terms, to indicate the expectation value.
    - `do_sigma`: (default=true) A boolean indicating whether to use $\sigma_x$ notation or $x$ notation for the terms.
"""
mutable struct diff_qEQ
    left_hand_side::qTerm
    right_hand_side::qEQ
    statespace::StateSpace
    braket::Bool
    do_sigma::Bool
    function diff_qEQ(left_hand_side::qTerm, right_hand_side::qEQ, statespace::StateSpace; braket::Bool=true, do_sigma::Bool=true)
        new(left_hand_side, neq(right_hand_side), statespace, braket, do_sigma)
    end
end

"""
    Sum(index::Union{String,Symbol,Vector{String},Vector{Symbol}}, expr::qEQ; neq::Bool=false) -> qSum

Constructor of a `qSum` struct. Defines the indexes to sum over, the expressions for which to apply the sum and optionally whether the sum is only over non equal indexes. 
"""
function Sum(indexes::Union{Vector{String},Vector{Symbol}}, expr::qEQ; neq::Bool=false)::qEQ
    index_strs = [string(index) for index in indexes]
    ss = expr.statespace
    the_s_ind::Int = -1
    e_inds::Vector{Int} = []
    for index_str in index_strs
        found = false
        for (s_ind, sub) in enumerate(ss.subspaces)
            for (e_ind, key) in enumerate(sub.keys)
                if key == index_str
                    if the_s_ind == -1
                        the_s_ind = s_ind
                    else
                        if s_ind != the_s_ind
                            error("Index $index_str found in multiple subspaces. Please specify a single subspace.")
                        end
                    end
                    push!(e_inds, e_ind)
                    found = true
                    break
                end
            end
            if found
                break
            end
        end
        if !found
            error("Index $index_str not found in any subspace keys in the state space.")
        end
    end
    if length(index_strs) > 0
        return qEQ([qSum(expr, index_strs, the_s_ind, e_inds, neq)], ss)
    else
        return expr
    end
end
function Sum(index::Union{String,Symbol}, expr::qEQ; neq::Bool=false)::qEQ
    return Sum([index], expr, neq=neq)
end
""" 
    ∑(index::Union{String,Symbol}, expr::qEQ; neq::Bool=false) -> qSum

Alternative way to call the `Sum` constructor. Sum(index, expr; neq) = ∑(index, expr; neq).
"""
∑(index::Union{String,Symbol}, expr::qEQ; neq::Bool=false) = Sum(index, expr, neq=neq)
∑(indexes::Union{Vector{String},Vector{Symbol}}, expr::qEQ; neq::Bool=false) = Sum(indexes, expr, neq=neq)

# Define iteration for qEQ so that iterating over it yields its qTerm's.
function iterate(q::qEQ, state::Int=1)
    state > length(q.terms) && return nothing
    return (q.terms[state], state + 1)
end

# Optionally, define length and eltype.
length(q::qEQ) = length(q.terms)
length(q::qSum) = length(q.expr.terms)
eltype(::Type{qEQ}) = qTerm
iszero(q::qTerm) = iszero(q.coeff)
iszero(q::qEQ) = length(q.terms) == 0 || all(iszero, q.terms)
iszero(q::qSum) = iszero(q.expr)

function copy(q::qTerm)::qTerm
    return qTerm(copy(q.coeff), copy(q.var_exponents), copy(q.op_indices))
end
function copy(q::qSum)::qSum
    return qSum(copy(q.expr), q.indexes, q.subsystem_index, q.element_indexes, q.neq)
end
function copy(q::qEQ)::qEQ
    return qEQ(copy(q.terms), q.statespace)
end

function var_exponents(q::qTerm)::Vector{Int}
    return q.var_exponents
end
function var_exponents(q::qSum)::Vector{Int}
    return zeros(Int, length(q.expr.statespace.vars))
end

function get_op_inds(space::SubSpace, res_strings::Tuple{Vector{Vector{String}},Vector{Vector{String}},Vector{Vector{String}}})::Vector{Vector{Tuple{Number,Is}}}
    str_elements::Vector{Vector{String}} = res_strings[space.statespace_main_ind]
    all_results::Vector{Vector{Tuple{Number,Is}}} = []
    for sub_ind in space.statespace_inds
        curr_elements::Vector{String} = str_elements[sub_ind]
        curr_results = space.op_set.strs2ind(curr_elements)
        push!(all_results, curr_results)
    end
    return all_results
end

"""
    term(q::StateSpace, operator_str::String)
    term(q::StateSpace, coeff, operator_str::String)

Generate a quantum term (qTerm) from the StateSpace `q`. The state description is provided as a string.

- Tokens of the form `var^exp` (e.g. `"a^2"`) set the exponent for a state variable.
- Other tokens are assumed to be keys that match one of the allowed subspace keys (i.e. elements in each SubSpace.keys).
- If no coefficient is given, the default coefficient is 1.
"""
function term(qspace::StateSpace, coeff::Number, operator_str::String="")::qEQ
    # Initialize exponents for each state variable.
    var_exponents = zeros(Int, length(qspace.vars))
    res_strings = separate_terms(operator_str, qspace.vars_str, qspace.fermionic_keys, qspace.bosonic_keys)
    for var_ind in 1:length(var_exponents)
        curr_res_str = res_strings[1][var_ind]
        var_exponents[var_ind] = sum([expstr_separate(curr_res_str[res])[2] for res in 1:length(curr_res_str)])
    end
    subspace_str2inds::Vector{Vector{Tuple{Number,Is}}} = []
    for space in qspace.subspaces
        append!(subspace_str2inds, get_op_inds(space, res_strings))
    end

    # Construct terms for each combination of subspace_str2inds
    terms::Vector{qTerm} = qTerm[]
    for combo in Iterators.product(subspace_str2inds...)
        curr_inds::Vector{Is} = []
        curr_coeff::Number = coeff
        for i in 1:length(combo)
            curr_coeff *= combo[i][1]
            push!(curr_inds, combo[i][2])
        end
        push!(terms, qTerm(curr_coeff, copy(var_exponents), curr_inds))
    end
    return qEQ(terms, qspace)

end
function term(qspace::StateSpace, operator_str::String)
    return term(qspace, one(1), operator_str)
end

include("qExpressionsOps/qExpressionsAlgebra.jl")

# Use your custom_sort_key for coefficients.
function qExpr_sort_key(term::qTerm)
    # Here we convert var_exponents (a Vector{Int}) to a tuple so that it compares lexicographically.
    return (tuple(term.var_exponents...), tuple(term.op_indices...), custom_sort_key(term.coeff), tuple(0, Int[], 0))
end
function qExpr_sort_key(term::qSum)
    curr_space = term.expr.statespace
    var_exponents = zeros(Int, length(curr_space.vars))
    op_indices::Vector{Is} = []
    for subspace in curr_space.subspaces
        op = subspace.op_set
        n_ops = length(subspace.keys)
        neutral_element = op.neutral_element
        for j in 1:length(n_ops)
            push!(op_indices, neutral_element)
        end
    end
    return (tuple(var_exponents...), tuple(op_indices...), custom_sort_key(0.0), tuple(term.subsystem_index, term.element_indexes, length(term)))
end
# Sort the terms in a qEQ using the key above.
function sort(qeq::qEQ)
    sorted_terms = sort(qeq.terms)
    return qEQ(sorted_terms, qeq.statespace)
end
function sort(qterms::Vector{qExpr})
    sorted_terms = sort(qterms, by=qExpr_sort_key)
    return sorted_terms
end

function same_term_type(t1, t2)
    return false
end
function same_term_type(t1::qTerm, t2::qTerm)::Bool
    return t1.var_exponents == t2.var_exponents && t1.op_indices == t2.op_indices
end
function same_term_type(s1::qSum, s2::qSum)::Bool
    return s1.subsystem_index == s2.subsystem_index && s1.element_indexes == s2.element_indexes
end

function combine_term(t1::qTerm, t2::qTerm)::qTerm
    return qTerm(t1.coeff + t2.coeff, t1.var_exponents, t1.op_indices)
end
function combine_term(s1::qSum, s2::qSum)::qSum
    return qSum(s1.expr + s2.expr, s1.indexes, s1.subsystem_index, s1.element_indexes, s1.neq)
end
"""
    simplify(q::qEQ) -> qEQ
    simplify(q::qSum) -> qSum
    simplify(q::diff_qEQ) -> diff_qEQ

Simplify a qEQ, qSum or diff_qEQ by sorting terms and ading up terms that are equal (up to a coefficient). 
"""
function simplify(q::qEQ)::qEQ
    # If there are no terms, return an empty qEQ.
    if isempty(q.terms)
        return qEQ(qExpr[], q.statespace)
    end

    # First, sort qEQ without modifying the original.
    sorted_q = sort(q)
    sorted_terms = copy(sorted_q.terms)

    combined_terms = qExpr[]
    i = 1
    curr_term = sorted_terms[1]
    curr_i = 1
    while i < length(sorted_terms)
        # Combine adjacent like terms.
        next_term = sorted_terms[i+1]
        if same_term_type(curr_term, next_term)
            curr_term = combine_term(curr_term, next_term)
        else
            if !iszero(curr_term)
                if isa(curr_term, qSum)
                    simplified_curr_term = simplify(curr_term)
                    if !iszero(simplified_curr_term)
                        push!(combined_terms, copy(simplified_curr_term))
                    end
                else
                    push!(combined_terms, copy(curr_term))
                end
                curr_term = next_term
                curr_i = i + 1
            end
        end
        i += 1
    end
    if !iszero(curr_term)
        if isa(curr_term, qSum)
            simplified_curr_term = simplify(curr_term)
            if !iszero(simplified_curr_term)
                push!(combined_terms, copy(simplified_curr_term))
            end
        else
            push!(combined_terms, copy(curr_term))
        end
    end
    return qEQ(combined_terms, q.statespace)
end
function simplify(s::qSum)::qSum
    simplified_expr = simplify(s.expr)
    return qSum(simplified_expr, s.indexes, s.subsystem_index, s.element_indexes, s.neq)
end

include("qExpressionsOps/qExpressionsPrint.jl")


"""
    base_operators(letter::String, qspace::StateSpace) -> Union{Int,Vector{Int}}
or 
    base_operators(ss:StateSpace) -> Tuple{Dict{String,qEQ},Dict{String,qEQ}}

Returns variables and/or operators in the state space `ss`.
Specifc variables/operators can be selected by passing a string `letter`.
If no `letter` is passed, the function returns a tuple of two dictionaries:
- The first dictionary contains the variables in the state space, with their corresponding qEQ objects.
- The second dictionary contains the operators in the state space, with their corresponding qEQ objects.

If you pass "vars", it will return a tuple with elements for each variable
"""
function base_operators(letter::String, qspace::StateSpace)
    my_ops::Vector{qEQ} = []
    var_exponents = zeros(Int, length(qspace.vars))
    curr_coeff = 1
    neutral_operator = [s.op_set.neutral_element for s in qspace.subspaces for key in s.keys]
    index = 1
    if letter == "I"
        return qEQ([qTerm(1, copy(var_exponents), copy(neutral_operator))], qspace)
    end
    if letter == "vars"
        for i in 1:length(qspace.vars)
            var_exponents[i] += 1
            push!(my_ops, qEQ([qTerm(curr_coeff, copy(var_exponents), copy(neutral_operator))], qspace))
            var_exponents[i] -= 1
        end
        if length(my_ops) > 1
            return tuple(my_ops...)
        elseif length(my_ops) == 1
            return my_ops[1]
        else
            return nothing
        end
    end
    for (i, var) in enumerate(qspace.vars)
        if var == letter
            var_exponents[i] += 1
            return qEQ([qTerm(curr_coeff, copy(var_exponents), copy(neutral_operator))], qspace)
        end
    end
    for (i, var) in enumerate(qspace.vars_str)
        if var == letter
            var_exponents[i] += 1
            return qEQ([qTerm(curr_coeff, copy(var_exponents), copy(neutral_operator))], qspace)
        end
    end
    # check variables
    for sub in qspace.subspaces
        for key in sub.keys
            if key == letter
                keys = sub.keys
                op_set = sub.op_set
                base_ops = op_set.base_ops
                for base_op in base_ops
                    curr_operator = copy(neutral_operator)
                    curr_operator[index] = base_op
                    term = qTerm(1, copy(var_exponents), curr_operator)
                    push!(my_ops, qEQ([term], qspace))
                end
                if length(my_ops) > 1
                    return tuple(my_ops...)
                else
                    return my_ops[1]
                end
            end
            index += 1
        end
    end
    error("No subspace with key starting with '$letter' found in the state space.")
end
function base_operators(qspace::StateSpace)::Tuple{Dict{String,qEQ},Dict{String,qEQ}}
    # return 2 dicctionaries, one with the vars and one with the operators 
    var_dict::Dict{String,qEQ} = Dict()
    op_dict::Dict{String,qEQ} = Dict()
    neutral_operator = [s.op_set.neutral_element for s in qspace.subspaces for key in s.keys]
    var_exponents = zeros(Int, length(qspace.vars))
    curr_coeff = 1
    for (i, vars_str) in enumerate(qspace.vars_str)
        var_exponents[i] += 1
        var_dict[vars_str] = qEQ([qTerm(curr_coeff, copy(var_exponents), copy(neutral_operator))], qspace)
        var_exponents[i] -= 1
    end
    index = 1
    for sub in qspace.subspaces
        op_set = sub.op_set
        base_ops = op_set.base_ops
        for key in sub.keys
            for base_op in base_ops
                curr_operator = copy(neutral_operator)
                curr_operator[index] = base_op
                term = qTerm(1, copy(var_exponents), curr_operator)
                curr_name = op_set.op2str(base_op, key)
                op_dict[curr_name] = qEQ([term], qspace)
            end
            index += 1
        end
    end
    op_dict["I"] = qEQ([qTerm(1, copy(var_exponents), copy(neutral_operator))], qspace)
    return var_dict, op_dict
end



##### Flatten 
"""
flatten(qeq::qEQ) -> qEQ

Flattens nested Sums in quantum Equations (qEQ).
"""
function flatten(qeq::qEQ)::qEQ
    new_terms = qExpr[]
    for t in qeq.terms
        if t isa qSum
            flat_eq = flatten_qSum(t)
            append!(new_terms, flat_eq.terms)
        else
            push!(new_terms, t)
        end
    end
    return qEQ(new_terms, qeq.statespace)
end

"""
flatten_qSum(s::qSum) -> qEQ

Take one `qSum` `s`.  First do `inner = flatten(s.expr)` so that all
deeper-nested sums are already one-level.  Split `inner.terms` into
• `base_terms` (just the `qTerm`’s)  
• `nested_sums` (any `qSum`’s).

Emit up to one “parent” sum over the `base_terms` (if non-empty), then
for each nested sum `n` emit a new `qSum` whose index-list is
`vcat(s.indexes, n.indexes)`.  Any duplicate index names will error.
"""
function flatten_qSum(s::qSum)::qEQ
    # first, fully flatten the body
    inner = flatten(s.expr)

    # pull out bare terms vs. sums
    base_terms = [t for t in inner.terms if t isa qTerm]
    nested_sums = [t for t in inner.terms if t isa qSum]

    out_terms = qExpr[]

    # if there were any base qTerms directly under `s`, keep a sum
    if !isempty(base_terms)
        push!(out_terms,
            qSum(qEQ(base_terms, inner.statespace),
                s.indexes, s.subsystem_index, s.element_indexes, s.neq))
    end

    # for each nested sum, merge its indexes onto `s`'s
    for n in nested_sums
        dup = intersect(s.indexes, n.indexes)
        if !isempty(dup)
            error("Unsupported: duplicate summation indexes detected: $(dup)")
        end

        merged_idxs = vcat(s.indexes, n.indexes)
        if n.subsystem_index != s.subsystem_index
            error("Unsupported: nested sum with different subsystem index: $(n.subsystem_index)")
        end
        merged_einds = vcat(s.element_indexes, n.element_indexes)
        push!(out_terms, qSum(n.expr, merged_idxs, s.subsystem_index, merged_einds, s.neq))
    end

    return qEQ(out_terms, inner.statespace)
end


function term_equal_indexes(expr::qTerm, index1::Int, index2::Int, subspace::SubSpace, coeff_inds1::Vector{Int}, coeff_inds2::Vector{Int})::Tuple{Bool,Vector{qTerm}}  # both elements not neutral, new term
    op1 = expr.op_indices[index1]
    op2 = expr.op_indices[index2]
    op_set = subspace.op_set
    neutral_element = op_set.neutral_element

    no_op = (op1 === neutral_element || op2 === neutral_element)
    no_coeff = all(expr.var_exponents[i] == 0 for i in coeff_inds1)
    if no_op && no_coeff
        # nothing changes: return `false` plus the original term
        return false, [expr]
    end

    # if they interact 
    results = op_set.op_product(op1, op2)
    new_terms = qTerm[]
    for (coeff, op) in results
        new_term = copy(expr)
        new_term.coeff *= coeff
        new_term.op_indices[index2] = op
        new_term.op_indices[index1] = neutral_element
        for (i, j) in zip(coeff_inds1, coeff_inds2)
            new_term.var_exponents[j] += new_term.var_exponents[i]
            new_term.var_exponents[i] = 0
        end
        push!(new_terms, new_term)
    end
    return true, new_terms
end

"""
    neq(qeq::qEQ) -> qEQ

Transform sums into neq sums, where all indexes are different from each other, and returns a flattened qEQ with neq sums. 
Considers all cases of the sums, simplifying the cases in which indexes are the same, which then reduces the order of the sum (i.e. a sum_{j} x_i y_j => sum_{j} x_i y_j + im*z_i, where we used x_i*y_i=im*z_i).
"""
function neq(qeq::qEQ)::qEQ
    # flatten first 
    qeq = flatten(qeq)
    out = qEQ(qTerm[], qeq.statespace)
    for t in qeq.terms
        if t isa qSum
            # expand this sum into distinct + diag parts
            out += neq_qsum(t)
        else
            out += t
        end
    end
    for t in out.terms
        if t isa qSum
            t.neq = true
        end
    end
    return out
end

# -------------------------------------------------------------------
# handle one qSum
function neq_qsum(s::qSum, index::Int=1)::qEQ
    if s.neq
        return qEQ([s], s.expr.statespace)   # skip
    end
    n = length(s.element_indexes) # is at least 1
    if n < index
        error("neq: index $index is out of range for this qSum (with n=$n)")
    end
    # 1) the “all distinct” piece # add to the lower part and remove it here 

    # 2) fwe consider for each sum index combination all possible 
    ss = s.expr.statespace
    sub = ss.subspaces[s.subsystem_index]
    n_sub = length(sub.statespace_inds)
    # consider only one possible equality, then recursively process untill all possibilities have been checked
    curr_element = s.element_indexes[index]
    if index < n # recursively execute neq_qsum for higher possible indexes
        post_expr = neq_qsum(s, index + 1)
    else
        post_expr = qEQ(qExpr[s], ss)
    end
    pieces = copy(post_expr)

    # 3) now we assume index is equal to each of the parameters in subspace, with smaller index than the curr index of the sum 
    curr_ind_sum::Int = s.element_indexes[index]
    curr_statespace_ind::Int = sub.statespace_inds[curr_ind_sum]

    coeffs_of_subspace = ss.where_by_continuum[s.subsystem_index]
    curr_coeff_inds = [coeffs_of_subspace[i][curr_ind_sum] for i in 1:length(coeffs_of_subspace)]
    #println("\nnew run: ($index) : ", pieces) 
    for (new_ind_sum, new_statespace_sum) in zip(1:curr_ind_sum-1, sub.statespace_inds[1:curr_ind_sum-1])
        new_coeff_inds = [coeffs_of_subspace[i][new_ind_sum] for i in 1:length(coeffs_of_subspace)]
        new_statespace_ind = sub.statespace_inds[new_ind_sum]
        # check for each term in the subspace if curr_statespace_ind and new_statespace_ind are the neutral_element  
        for expr in post_expr.terms
            if isa(expr, qTerm)
                not_neutral, new_terms = term_equal_indexes(expr, curr_statespace_ind, new_statespace_ind, sub, curr_coeff_inds, new_coeff_inds)
                if not_neutral
                    error("Unsupported: Element that isn't part of a Sum should no longer contain sum indexes")
                end
                for new_term in new_terms
                    pieces += new_term
                end
            elseif isa(expr, qSum)
                #println("    ($index) - expr: ", expr)
                for t in expr.expr.terms
                    not_neutral, new_terms = term_equal_indexes(t, curr_statespace_ind, new_statespace_ind, sub, curr_coeff_inds, new_coeff_inds)
                    # add new terms as qSum(s) with corrected indexing 
                    if not_neutral
                        # remove expr.indexes[index] and similarly expr.element_indexes[index]
                        new_indexes = vcat(expr.indexes[1:index-1], expr.indexes[index+1:end])
                        new_element_indexes = vcat(expr.element_indexes[1:index-1], expr.element_indexes[index+1:end])
                        if length(new_indexes) == 0
                            for new_term in new_terms
                                pieces += new_term
                            end
                        else
                            pieces += qSum(qEQ(new_terms, ss), new_indexes, expr.subsystem_index, new_element_indexes, expr.neq)
                        end
                    else # no change to sum structure
                        pieces += qSum(qEQ(new_terms, ss), copy(expr.indexes), expr.subsystem_index, copy(expr.element_indexes), expr.neq)
                    end
                end
            end
        end
        #println("  Result for ($index => $curr_ind_sum, $new_ind_sum | $curr_statespace_ind, $new_statespace_ind):  " , pieces)
    end
    return simplify(pieces)
end


function simplify(q::diff_qEQ)::diff_qEQ
    simp_rhs = simplify(q.right_hand_side)
    return diff_qEQ(q.left_hand_side, simp_rhs, q.statespace, q.braket, q.do_sigma)
end
# Import the base addition operator.

"""
    d_dt(qspace::StateSpace, expr)

Evaluate the time derivative of an expression `expr` in the context of the given state space `ss`.

This function expects that `expr` is an equation (i.e. an Expr with an equal sign as its head),
of the form

    LHS = RHS
The function then returns a `diff_qEQ` constructed from the left-hand side qTerm and the right-hand side qEQ.
"""
function d_dt(left_hand::Union{qTerm,qEQ}, right_hand::qEQ)::diff_qEQ
    # Check if expr is an equality.
    qstate = right_hand.statespace

    if left_hand isa qEQ
        if left_hand.statespace != qstate
            error("Left and right sides of the equation must be in the same state space.")
        end
        if length(left_hand.terms) != 1
            error("Left-hand side of the equation must consist of a single qTerm.")
        end
        left_hand = left_hand.terms[1]
    end
    if abs(left_hand.coeff - 1) > 1e-10
        error("Left-hand side of the equation must be a qTerm with coeff 1.")
    end
    if !iszero(left_hand.var_exponents)
        error("Left-hand side of the equation must be a qTerm with no variable exponents.")
    end
    # Return a diff_qEQ constructed from these sides.
    return diff_qEQ(left_hand, right_hand, qstate)
end
end