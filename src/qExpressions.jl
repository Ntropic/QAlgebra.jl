module qExpressions
using ..qSpace
using ..FFunctions
using ..StringUtils
using ComplexRationals
import Base: show, adjoint, iterate, length, eltype, +, -, sort, *, ^, product, iszero, copy

export qObj, qAtom, qComposite, qTerm, qExpr, qSum, Sum, ∑, diff_qEQ, term, base_operators,simplify, flatten, neq, d_dt


# ==========================================================================================================================================================
# --------> TBase Types and Their Constructors <---------------------------------------------------------------------------------------------------------
# ==========================================================================================================================================================

""" 
    qObj

The abstract type `qObj` is the base type for all quantum expressions in this module.
"""
abstract type qObj end  # most general
""" 
    qAtom

The abstract type `qAtom` is a subtype of `qObj` and represents elementary operator definitions, such as qTerm and qAbstract. 
"""
abstract type qAtom <: qObj end # elementary operator definitions 
""" 
    qComposite

The abstract type `qComposite` is a subtype of `qObj` and represents composite expressions, 
such as qSum and qProd which consist of qAtom, qAbstract or qComposite objects themselves. 
"""
abstract type qComposite <: qObj end  # products and sums of operator definitions

Is = Union{Int,Vector{Int}}

"""
    qTerm

A `qTerm` represents a single term in a quantum expression. It contains:
    - `op_indices`: A vector of indices representing the operators in the term, which are also defined in a StateSpace.
"""
mutable struct qTerm <: qAtom
    op_indices::Vector{Is}
end

"""
    qAbstract(indices::Vector{Int})

A purely‐symbolic abstract operator
    - key_index: The index of the abstract_key in the state space.
    - sub_index: The index of the suboperator in the state_space 
    - exponent: The exponent of the operator.
    - dag: A boolean indicating whether the operator is daggered (default = `false`)
    - operator_type: A reference to the operator of which it is a type 
    - index_map: Keeps track of indexes, that are equal (for neq transformations)
Is an instance of an OperatorType 
"""
mutable struct qAbstract <: qAtom
    key_index::Int
    sub_index::Int
    exponent::Int
    dag::Bool
    operator_type::OperatorType
    index_map::Vector{Tuple{Int,Int}}
end
function qAbstract(operator_type::OperatorType, key_index::Int, sub_index::Int, exponent::Int, dag::Bool; index_map::Vector{Tuple{Int,Int}}=Tuple{Int,Int}[])
    return qAbstract(key_index, sub_index, exponent, dag, operator_type, index_map)
end

""" 
    qProd

A product of qAtom expressions, i.e. qTerms or qAbstract.
It contains:
    - `statespace`: The state space in which the product is defined.
    - `coeff_fun`: The function of parameters for the Operator product
    - `expr`: A vector of qAtoms (qTerms or qAbstract) that are multiplied together.
"""
mutable struct qProd <: qComposite
    statespace::StateSpace         # State space of the product.
    coeff_fun::FFunction            # function of scalar parameters => has +,-,*,/,^ defined 
    expr::Vector{qAtom}             # Vector of qAtoms (qTerms or qAbstract).
    function qProd(statespace::StateSpace, coeff::FFunction, expr::Vector{<:qAtom}= [])
        if isempty(expr)
           # add neutral operator
           expr = [statespace.neutral_op]
        end
        new(statespace, coeff, expr)
    end
    function qProd(statespace::StateSpace, coeff::Number, var_exponents::Vector{Int}, expr::qAtom)
        f_fun = FAtom(coeff, var_exponents) 
        return new(statespace, f_fun, [expr]) 
    end
    function qProd(statespace::StateSpace, coeff::Number, var_exponents::Vector{Int}, expr::AbstractVector{<:qAtom})
        f_fun = FAtom(coeff, var_exponents)
        return new(statespace, f_fun, expr)
    end
    function qProd(statespace::StateSpace, coeff::Number, expr::AbstractVector{<:qAtom})
        f_fun = coeff* statespace.fone
        return new(statespace, f_fun, expr)
    end
    function qProd(statespace::StateSpace, coeff::Number, expr::qAtom)
        f_fun = coeff * statespace.fone
        return new(statespace, f_fun, [expr]) 
    end
end

"""
    qExpr

A `qExpr` represents a quantum equation, consisting of a Vector of quantum Expressions representing the additive terms of the equation.
It also contains a reference to the state space in which the equation is defined.
"""
mutable struct qExpr
    statespace::StateSpace
    terms::Vector{qComposite}         # Vector of terms
    function qExpr(statespace::StateSpace, terms::Vector{<:qComposite})
        if isempty(terms) 
            # add neotral zero term
            zero_term = qProd(statespace.fone*0, qTerm(statespace.neutral_op))
            push!(terms, zero_term)
        end
        return new(statespace, terms)
    end
end
function qExpr(statespace::StateSpace, prod::qComposite)
    return qExpr(statespace, [prod])
end
function qExpr(statespace::StateSpace, terms::qAtom)
    return qExpr(statespace, qProd(statespace, statespace.fone, [terms]))
end

""" 
    qSum

A `qSum` represents the summation of a quantum Equation over indexes in a quantum expression.
It contains:
    - `expr`: The expression being summed over, which is a `qExpr` object.
    - `indexes`: A vector of strings representing the summation indexes (e.g., "i").
    - `subsystem_index`: The index of the subspace in which the indexes live. 
    - `element_indexes`: A vector of integers representing the position of the indexes in that subspace.
    - `neq`: A boolean indicating whether different indexes in the sum can refer to the same element in the subspace. 
            For example, the indexes i,j,k can refer to different elements in a much larger bath of elements. 
"""
mutable struct qSum <: qComposite
    statespace::StateSpace
    expr::qExpr       # The expression being summed over.    # use expr in other qComposites except for qProd
    indexes::Vector{String}   # The summation index (e.g. "i").
    subsystem_index::Int  # The subspace index where the summation index was found.
    element_indexes::Vector{Int}    # The position in that subspace.
    neq::Bool
end

"""
    diff_qEQ

A `diff_qEQ` represents a differential equation of the form:

    d/dt ⟨Op⟩ = RHS

It represents time evolution of operator expectation values, and wraps the symbolic structure of such an equation.

# Fields
- `left_hand_side::qTerm`: The LHS operator being differentiated.
- `expr::qExpr`: The RHS symbolic expression.
- `statespace::StateSpace`: The StateSpace in which the equation is defined.
- `braket::Bool`: Whether to use braket notation ⟨⋯⟩ (default = `true`).
- `do_sigma::Bool`: Whether to display Pauli operators as `σₓ`, etc. (default = `true`).
"""
mutable struct diff_qEQ <: qComposite
    left_hand_side::qProd
    expr::qExpr
    statespace::StateSpace
    braket::Bool
    do_sigma::Bool
end

"""
    diff_qEQ(lhs::qTerm, rhs::qExpr, statespace::StateSpace; braket=true, do_sigma=true)

Construct a [`diff_qEQ`](@ref) that represents the time derivative of ⟨lhs⟩ = rhs.

Automatically applies `neq()` to the RHS to expand sums over distinct indices.
"""
function diff_qEQ(left_hand_side::qProd, expr::qExpr, statespace::StateSpace; braket::Bool=true, do_sigma::Bool=true)
    new_rhs = neq(expr)
    return diff_qEQ(left_hand_side, new_rhs, statespace, braket, do_sigma)
end

"""
    Sum(index::Union{String,Symbol,Vector{String},Vector{Symbol}}, expr::qExpr; neq::Bool=false) -> qSum

Constructor of a `qSum` struct. Defines the indexes to sum over, the expressions for which to apply the sum and optionally whether the sum is only over non equal indexes. 
"""
function Sum(indexes::Union{Vector{String},Vector{Symbol}}, expr::qExpr; neq::Bool=false)::qExpr
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
        return qExpr(ss, [qSum(ss, expr, index_strs, the_s_ind, e_inds, neq)])
    else
        return expr
    end
end
function Sum(index::Union{String,Symbol}, expr::qExpr; neq::Bool=false)::qExpr
    return Sum([index], expr, neq=neq)
end
""" 
    ∑(index::Union{String,Symbol}, expr::qExpr; neq::Bool=false) -> qSum

Alternative way to call the `Sum` constructor. Sum(index, expr; neq) = ∑(index, expr; neq).
"""
∑(index::Union{String,Symbol}, expr::qExpr; neq::Bool=false) = Sum(index, expr, neq=neq)
∑(indexes::Union{Vector{String},Vector{Symbol}}, expr::qExpr; neq::Bool=false) = Sum(indexes, expr, neq=neq)
 
#### Helper Functions #######################################################################################
# Define iteration for qExpr so that iterating over it yields its qTerm's.
function iterate(q::qExpr, state::Int=1)
    state > length(q.terms) && return nothing
    return (q.terms[state], state + 1)
end

# Optionally, define length and eltype.
length(q::qExpr) = length(q.terms)
length(q::qSum) = length(q.expr.terms)
length(q::qProd) = length(q.expr)
iszero(q::qProd) = iszero(q.coeff_fun)
iszero(q::qExpr) = length(q.terms) == 0 || all(iszero, q.terms)
iszero(q::qSum) = iszero(q.expr)

function copy(q::qTerm)::qTerm
    return qTerm(copy(q.op_indices))
end
function copy(q::qAbstract)::qAbstract
    return qAbstract(copy(q.key_index), copy(q.sub_index), copy(q.exponent), copy(q.dag), copy(q.operator_type), copy(q.index_map))
end
function copy(q::qProd)::qProd
    return qProd(q.statespace, copy(q.coeff_fun), copy(q.expr))
end
function copy(q::qSum)::qSum
    return qSum(copy(q.statespace), copy(q.expr), copy(q.indexes), copy(q.subsystem_index), copy(q.element_indexes), copy(q.neq))
end
function copy(q::qExpr)::qExpr
    return qExpr(q.statespace, copy(q.terms))
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

# ==========================================================================================================================================================
# --------> Terms from Strings and Base Operators <---------------------------------------------------------------------------------------------------------
# ==========================================================================================================================================================

function string2qterm(statespace::StateSpace, operator_str::String="")::Tuple{Vector{qTerm}, Vector{Number}, Vector{Int}}
    # Initialize exponents for each state variable.
    var_exponents = zeros(Int, length(statespace.vars))
    res_strings = separate_terms(operator_str, statespace.vars_str, statespace.fermionic_keys, statespace.bosonic_keys)
    for var_ind in 1:length(var_exponents)
        curr_res_str = res_strings[1][var_ind]
        var_exponents[var_ind] = sum([expstr_separate(curr_res_str[res])[2] for res in 1:length(curr_res_str)])
    end
    subspace_str2inds::Vector{Vector{Tuple{Number,Is}}} = []
    for space in statespace.subspaces
        append!(subspace_str2inds, get_op_inds(space, res_strings))
    end
    # Construct terms for each combination of subspace_str2inds
    terms::Vector{qTerm} = []
    coeffs::Vector{Number} = []
    for combo in Iterators.product(subspace_str2inds...)
        curr_inds::Vector{Is} = []
        curr_coeff::Number = 1
        for i in 1:length(combo)
            curr_coeff *= combo[i][1]
            push!(curr_inds, combo[i][2])
        end
        push!(terms, qTerm(curr_inds))
        push!(coeffs, curr_coeff)
    end
    return terms, coeffs, var_exponents
end
function string2qabstract(statespace::StateSpace, operator_str::String)::qAbstract
    # find which statespace.operator_names fits the beginning of the string 
    name = nothing
    rest::String = ""
    key_index = -1
    for (i, op_name) in enumerate(statespace.operator_names)
        if startswith(operator_str, op_name)
            name = op_name
            key_index = i 
            rest = operator_str[length(op_name)+1:end]
            break
        end
    end
    if isnothing(name)
        error("Invalid string: $operator_str, must be one of $(statespace.operator_names)")
    end
    # process rest separate by separating into pre and post ^ 
    expstr, exp = expstr_separate(rest)
    # if expstr ends with ' conjugate = true 
    conjugate = occursin("'", expstr) 
    subindex = -1
    if conjugate
        a,b = split(expstr, "'")
        if length(a) > 0
            subindex = parse(Int, a)
        end
        # find the operator in statespace 
        if length(b) > 0
            error("For an abstract operator, the Dagger symbol (') can only be followed up by an exponential, nothing else.")
        end
    end
    operator_type = statespace.operatortypes[key_index]
    return qAbstract(operator_type, key_index, subindex, exp, conjugate) # subindex only printed if not -1 
end

"""
    term(operator_str::String, statespace::StateSpace)
    term(coeff, operator_str::String, statespace::StateSpace)
    term(coeff::Number, operator_str::String)
    term(operator_str::String)

Generate a quantum term (qTerm) from the StateSpace `q`. The state description is provided as a string.

- Tokens of the form `var^exp` (e.g. `"a^2"`) set the exponent for a state variable.
- Other tokens are assumed to be keys that match one of the allowed subspace keys (i.e. elements in each SubSpace.keys) or 
  are abstract operator names in the StateSpace (with potential indexes, exponents and Daggers ').
- If no coefficient is given, the default coefficient is 1.
- Daggers are given by ', and need to preceed a possible exponent
"""
function term(statespace::StateSpace, coeff::Number, operator_str::String)
    out_strings, out_type = term_pre_split(operator_str, statespace.operator_names)
    var_exponents = zeros(Int, length(statespace.vars))
    coeffs_and_terms::Vector{Vector{Tuple{Number, qAtom}}} = []
    for (out_string, type) in zip(out_strings, out_type)
        if type # qAbstract
            curr_term = string2qabstract(statespace, out_string)
            push!(coeffs_and_terms, [(1, curr_term)])
        else  # qTerm 
            curr_terms, curr_coeffs, curr_var_exponents = string2qterm(statespace, out_string)
            curr_coeffs_and_terms = []
            for (curr_coeff, curr_term) in zip(curr_coeffs, curr_terms)
                push!(curr_coeffs_and_terms, (curr_coeff, curr_term))
            end
            push!(coeffs_and_terms, curr_coeffs_and_terms)
            var_exponents .+= curr_var_exponents
        end
    end
    # Generate all combinations using Product 
    products::Vector{qProd} = []
    for comb in Iterators.product(coeffs_and_terms...)
        curr_coeff = reduce(*, [c[1] for c in comb])*coeff
        curr_atoms = [c[2] for c in comb]
        push!(products, qProd(statespace, curr_coeff, copy(var_exponents), curr_atoms))
    end
    return qExpr(statespace, products)
end
function term(statespace::StateSpace, operator_str::String)::qExpr
    return term(statespace, one(1), operator_str)
end
function term(operator_str::String)::qExpr
    return term(GLOBAL_STATE_SPACE, operator_str)
end
function term(coeff::Number, operator_str::String)::qExpr
    return term(GLOBAL_STATE_SPACE, coeff, operator_str)
end

"""
    base_operators(letter::String, statespace::StateSpace) -> Union{Int,Vector{Int}}
    base_operators(ss:StateSpace) -> Tuple{Dict{String,qExpr},Dict{String,qExpr},Dict{String,Function}}

Returns variables and/or operators in the state space `ss`.
Specifc variables/operators can be selected by passing a string `letter`.
If no `letter` is passed, the function returns a tuple of 3 dictionaries:
- The first dictionary contains the variables in the state space, with their corresponding qExpr objects.
- The second dictionary contains the operators in the state space, with their corresponding qExpr objects.
- The third dictionary contains the abstract operators in the state space.
If you pass "vars", it will return a tuple with elements for each variable
"""
function base_operators(letter::String, statespace::StateSpace)
    my_ops::Vector{qExpr} = []
    var_exponents = zeros(Int, length(statespace.vars))
    neutral_operator = [s.op_set.neutral_element for s in statespace.subspaces for key in s.keys]
    index = 1
    if letter == "I"
        return qExpr(statespace, qTerm(copy(neutral_operator)))
    end
    if letter == "vars"
        for i in 1:length(statespace.vars)
            var_exponents[i] += 1
            push!(my_ops, qExpr(qProd(statespace, 1, var_exponents, qTerm(copy(neutral_operator))), statespace))
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
    for (i, var) in enumerate(statespace.vars)
        if var == letter
            var_exponents[i] += 1
            return qExpr(statespace, qProd(statespace, 1, var_exponents, qTerm(copy(neutral_operator))))
        end
    end
    for (i, var) in enumerate(statespace.vars_str)
        if var == letter
            var_exponents[i] += 1
            return qExpr(statespace, qProd(statespace, 1, var_exponents, qTerm(copy(neutral_operator))))
        end
    end
    # check subspaces
    for sub in statespace.subspaces
        for key in sub.keys
            if key == letter
                keys = sub.keys
                op_set = sub.op_set
                base_ops = op_set.base_ops
                for base_op in base_ops
                    curr_operator = copy(neutral_operator)
                    curr_operator[index] = base_op
                    push!(my_ops, qExpr(statespace, qTerm(copy(curr_operator))))
                end

                # non base ops # in a Dict 
                non_base_ops = op_set.non_base_ops
                for (key, ops) in non_base_ops
                    curr_terms::Vector{qProd} = qProd[]
                    for op in ops
                        curr_operator = copy(neutral_operator)
                        curr_operator[index] = op[2]
                        coeff = op[1]
                        curr_prod = qProd(statespace,coeff, var_exponents, qTerm(copy(curr_operator)))
                        push!(curr_terms, curr_prod)
                    end
                    push!(my_ops, qExpr(statespace, curr_terms))
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
    # check for abstract operators
    for (key_index, name) in enumerate(statespace.operator_names)
        if name == letter
            return (subindex=-1) -> qExpr(statespace, qAbstract(key_index, subindex, 1, false))
        elseif name in letter
            abstract = string2qabstract(statespace, replace(letter, "_" => ""))
            return qExpr(statespace, abstract)
        end
    end
    error("No variable, subspace component or abstract operator with key starting with '$letter' found in the state space.")
end
function base_operators(statespace::StateSpace)::Tuple{Dict{String,qExpr},Dict{String,qExpr}, Dict{String, Function}}
    # return 2 dicctionaries, one with the vars and one with the operators 
    var_dict::Dict{String,qExpr} = Dict()
    op_dict::Dict{String,qExpr} = Dict()
    CRone = one(ComplexRational)
    var_exponents = zeros(Int, length(statespace.vars))
    #abstract_dict::Dict{String,qExpr} = Dict()
    neutral_operator = s.neutral_op
    for (i, vars_str) in enumerate(statespace.vars_str)
        var_exponents[i] += 1
        var_dict[vars_str] = qExpr(statespace, qProd(statespace, CRone, var_exponents, qTerm(copy(neutral_operator))))
        var_exponents[i] -= 1
    end
    index = 1
    for sub in statespace.subspaces
        op_set = sub.op_set
        base_ops = op_set.base_ops
        for key in sub.keys
            for base_op in base_ops
                curr_operator = copy(neutral_operator)
                curr_operator[index] = base_op
                term = qTerm(curr_operator)
                curr_name = op_set.op2str(base_op, key)
                op_dict[curr_name] = qExpr(statespace, term)
            end

            # non base ops # in a Dict 
            non_base_ops = op_set.non_base_ops
            for (key, ops) in non_base_ops
                curr_terms::Vector{qProd} = qProd[]
                for op in ops
                    curr_operator = copy(neutral_operator)
                    curr_operator[index] = op[2]
                    coeff = op[1]
                    curr_prod = qProd(statespace, coeff, var_exponents, qTerm(copy(curr_operator)))
                    push!(curr_terms, curr_prod)
                end
                op_dict[key] = qExpr(statespace, curr_terms)
            end
            index += 1
        end
    end
    op_dict["I"] = qExpr(statespace, qTerm(copy(neutral_operator)))
    # now abstract operators 
    abstract_dict::Dict{String, Function} = Dict{String, Function}()
    for (key_index, name) in enumerate(statespace.operator_names)
        abstract_dict[name] = (subindex=-1) -> qExpr(statespace, qAbstract(key_index, subindex, 1, false))
    end
    return var_dict, op_dict, abstract_dict
end
# ===========================================================================================================================================================
# --------> Sorting and Simplifying <------------------------------------------------------------------------------------------------------------------------
# ===========================================================================================================================================================

function qAtom_sort_key(term::qTerm)
    # Here we convert var_exponents (a Vector{Int}) to a tuple so that it compares lexicographically.
    return tuple(0, term.op_indices...)
end
function qAtom_sort_key(term::qAbstract)
    return tuple(1, term.key_index, term.subindex, term.exponent, Int(term.dag))
end

# Use your custom_sort_key for coefficients.
function qobj_sort_key(term::qProd; var_first::Bool=false)
    # Here we convert var_exponents (a Vector{Int}) to a tuple so that it compares lexicographically.
    if var_first
        return (tuple(term.coeff_fun.var_exponents...), tuple(term.op_indices...), custom_sort_key(term.coeff_fun.coeff), tuple(0, Int[], 0))
    else
        return (tuple(term.op_indices...), tuple(term.coeff_fun.var_exponents...), custom_sort_key(term.coeff_fun.coeff), tuple(0, Int[], 0))
    end
end
function qobj_sort_key(term::qSum; var_first::Bool=false)
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
    if var_first
        return (tuple(var_exponents...), tuple(op_indices...), custom_sort_key(0.0), tuple(term.subsystem_index, term.element_indexes, length(term)))
    else
        return (tuple(op_indices...), tuple(var_exponents...), custom_sort_key(0.0), tuple(term.subsystem_index, term.element_indexes, length(term)))
    end
end
# Sort the terms in a qExpr using the key above.
function sort(qeq::qExpr; var_first::Bool=false, kwargs...)
    sorted_terms = sort(qeq.terms, by= x-> qobj_sort_key(x; var_first=var_first), kwargs...)
    return qExpr(qeq.statespace, sorted_terms)
end
function sort(qterms::Vector{qComposite}; var_first::Bool=false, kwargs...)
    sorted_terms = sort(qterms, by= x-> qobj_sort_key(x; var_first=var_first), kwargs...)
    return sorted_terms
end


function same_term_type(t1, t2)
    return false
end
function same_term_type(t1::qProd, t2::qProd)::Bool
    return t1.coeff_fun.var_exponents == t2.coeff_fun.var_exponents && t1.op_indices == t2.op_indices
end
function same_term_type(s1::qSum, s2::qSum)::Bool
    return s1.subsystem_index == s2.subsystem_index && s1.element_indexes == s2.element_indexes
end

function combine_term(t1::qProd, t2::qProd)::qTerm
    return qProd(simplify(t1.coeff_fun + t2.coeff_fun), t1.expr)
end
function combine_term(s1::qSum, s2::qSum)::qSum
    return qSum(s1.statespace, simplify(s1.expr + s2.expr), s1.indexes, s1.subsystem_index, s1.element_indexes, s1.neq)
end

"""
    simplify(q::qProd, statespace::StateSpace) -> qProd
    simplify(q::qExpr) -> qExpr
    simplify(q::qSum) -> qSum
    simplify(q::diff_qEQ) -> diff_qEQ

Simplify a qProd, qExpr, qSum or diff_qEQ by sorting terms and ading up terms that are equal (up to a coefficient). 
"""
function simplify(q::qProd, statespace::StateSpace)::qProd   # remove this statespace, continue here
    # can'T just sort terms in qProd, we need to check for each term: if two adjacent terms are both qTerm => compute their product. for qAbstract terms, we can check ideal order
end

function simplify(q::qExpr)::qExpr
    # If there are no terms, return an empty qExpr.
    if isempty(q.terms)
        return qExpr(q.statespace, qComposite[])
    end

    # First, sort qExpr without modifying the original.
    sorted_q = sort(q)
    sorted_terms = copy(sorted_q.terms)

    combined_terms = qComposite[]
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
    return qExpr(q.statespace, combined_terms)
end
function simplify(s::qSum)::qSum
    simplified_expr = simplify(s.expr)
    return qSum(s.statespace, simplified_expr, s.indexes, s.subsystem_index, s.element_indexes, s.neq)
end


include("qExpressionsOps/qExpressionsAlgebra.jl")
include("qExpressionsOps/qExpressionsPrint.jl")

##### Flatten 
"""
flatten(qeq::qExpr) -> qExpr

Flattens nested Sums in quantum Equations (qExpr).
"""
function flatten(qeq::qExpr)::qExpr
    new_terms = qComposite[]
    for t in qeq.terms
        if t isa qSum
            flat_eq = flatten_qSum(t)
            append!(new_terms, flat_eq.terms)
        else
            push!(new_terms, t)
        end
    end
    return qExpr(qeq.statespace, new_terms)
end

"""
flatten_qSum(s::qSum) -> qExpr

Take one `qSum` `s`.  First do `inner = flatten(s.expr)` so that all
deeper-nested sums are already one-level.  Split `inner.terms` into
• `base_terms` (just the `qTerm`’s)  
• `nested_sums` (any `qSum`’s).

Emit up to one “parent” sum over the `base_terms` (if non-empty), then
for each nested sum `n` emit a new `qSum` whose index-list is
`vcat(s.indexes, n.indexes)`.  Any duplicate index names will error.
"""
# No docstring here
function flatten_qSum(s::qSum)::qExpr
    # first, fully flatten the body
    inner = flatten(s.expr)

    # pull out bare terms vs. sums
    base_terms = [t for t in inner.terms if t isa qTerm]
    nested_sums = [t for t in inner.terms if t isa qSum]

    out_terms = qComposite[]

    # if there were any base qTerms directly under `s`, keep a sum
    if !isempty(base_terms)
        push!(out_terms,
            qSum(inner.statespace, qExpr(base_terms),
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
        push!(out_terms, qSum(s.statespace, n.expr, merged_idxs, s.subsystem_index, merged_einds, s.neq))
    end

    return qExpr(inner.statespace, out_terms)
end

# change from index1 to index2
function term_equal_indexes(expr, args...)
    throw(MethodError(term_equal_indexes, (typeof(expr), args...)))
end
# multiplies from the left 
function term_equal_indexes(term::qTerm, index1::Int, index2::Int, subspace::SubSpace)::Tuple{Bool, Vector{qTerm}}
    op1 = term.op_indices[index1]
    op2 = term.op_indices[index2]
    neutral = subspace.op_set.neutral_element
    if op1 === neutral && op2 === neutral
        return false, [term]
    end
    results = subspace.op_set.op_product(op1, op2)
    new_terms = qTerm[]
    for (coeff, op) in results
        new_term = deepcopy(term)
        new_term.op_indices[index2] = op
        new_term.op_indices[index1] = neutral
        push!(new_terms, new_term)
    end

    return true, new_terms
end

function term_equal_indexes(abstract::qAbstract, index1::Int, index2::Int, subspace::SubSpace)::Tuple{Bool, Vector{qAbstract}}
    op_type = abstract.operator_type
    non_trivial_op_indices = op_type.non_trivial_op_indices
    if non_trivial_op_indices[index2]
        new_abstract = copy(abstract)
        push!(new_abstract.index_map, (index1, index2))
        return true, [new_abstract]
    end
    # append this rule to the index map 
    return false, [abstract]
end

function term_equal_indexes(prod::qProd, index1::Int, index2::Int, subspace::SubSpace)::Tuple{Bool, Vector{qProd}}
    changed_any = false
    term_variants = Vector{Vector{qAtom}}()
    for atom in prod.expr
        changed, variants = term_equal_indexes(atom, index1, index2, subspace)
        push!(term_variants, variants)
        changed_any |= changed  # Check if any term was changed
    end

    if !changed_any
        return false, [prod]
    end

    # Generate all combinations (cartesian product) of updated terms
    combinations = Iterators.product(term_variants...)

    simplified_products = qProd[]

    for combo in combinations
        new_expr = collect(combo)
        new_prod = qProd(prod.statespace, prod.coeff_fun, new_expr)
        push!(simplified_products, simplify(new_prod))
    end

    return true, simplified_products
end


"""
    neq(qeq::qExpr) -> qExpr

Transform sums into neq sums, where all indexes are different from each other, and returns a flattened qExpr with neq sums. 
Considers all cases of the sums, simplifying the cases in which indexes are the same, which then reduces the order of the sum (i.e. a sum_{j} x_i y_j => sum_{j} x_i y_j + im*z_i, where we used x_i*y_i=im*z_i).
"""
function neq(qeq::qExpr)::qExpr
    # flatten first 
    qeq = flatten(qeq)
    out = qExpr(qeq.statespace, qTerm[])
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
function neq_qsum(s::qSum, index::Int=1)::qExpr
    if s.neq
        return qExpr(s.expr.statespace, [s])   # skip
    end
    n = length(s.element_indexes) # is at least 1
    if n < index
        error("neq: index $index is out of range for this qSum (with n=$n)")
    end
    # 1) the “all distinct” piece # add to the lower part and remove it here 

    # 2) we consider for each sum index combination all possible 
    ss = s.expr.statespace
    sub = ss.subspaces[s.subsystem_index]
    n_sub = length(sub.statespace_inds)
    # consider only one possible equality, then recursively process untill all possibilities have been checked
    curr_element = s.element_indexes[index]
    if index < n # recursively execute neq_qsum for higher possible indexes
        post_expr = neq_qsum(s, index + 1)
    else
        post_expr = qExpr(ss, qComposite[s])
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
                            pieces += qSum(s.statespace, qExpr(ss, new_terms), new_indexes, expr.subsystem_index, new_element_indexes, expr.neq)
                        end
                    else # no change to sum structure
                        pieces += qSum(s.statespace, qExpr(ss, new_terms), copy(expr.indexes), expr.subsystem_index, copy(expr.element_indexes), expr.neq)
                    end
                end
            end
        end
        #println("  Result for ($index => $curr_ind_sum, $new_ind_sum | $curr_statespace_ind, $new_statespace_ind):  " , pieces)
    end
    return simplify(pieces)
end


function simplify(q::diff_qEQ)::diff_qEQ
    simp_rhs = simplify(q.expr)
    return diff_qEQ(q.left_hand_side, simp_rhs, q.statespace, q.braket, q.do_sigma)
end
function copy(q::diff_qEQ)::diff_qEQ
    return diff_qEQ(copy(q.left_hand_side), copy(q.expr), copy(q.statespace), copy(q.braket), copy(q.do_sigma))
end

"""
    d_dt(statespace::StateSpace, expr)

Evaluate the time derivative of an expression `expr` in the context of the given state space `ss`.

This function expects that `expr` is an equation (i.e. an Expr with an equal sign as its head),
of the form

    LHS = RHS
The function then returns a `diff_qEQ` constructed from the left-hand side qTerm and the right-hand side qExpr.
"""
function d_dt(left_hand::Union{qProd,qExpr}, right_hand::qExpr)::diff_qEQ
    # Check if expr is an equality.
    qstate = right_hand.statespace

    if left_hand isa qExpr
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
    if !iszero(left_hand.coeff_fun.var_exponents)
        error("Left-hand side of the equation must be a qTerm with no variable exponents.")
    end
    # Return a diff_qEQ constructed from these sides.
    return diff_qEQ(left_hand, right_hand, qstate)
end
end