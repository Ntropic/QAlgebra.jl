export term

# ==========================================================================================================================================================
# --------> Terms from Strings <----------------------------------------------------------------------------------------------------------------------------
# ==========================================================================================================================================================
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
function string2qterm(statespace::StateSpace, operator_str::String="")::Tuple{Vector{QTerm}, Vector{Number}, Vector{Int}}
    # Initialize exponents for each state variable.
    var_exponents = zeros(Int, length(statespace.vars))
    var_strs = statespace.vars_str
    alt_var_strs = [p.var_name for p in statespace.vars]
    res_strings = separate_terms(operator_str, var_strs, alt_var_strs, statespace.fermionic_keys, statespace.bosonic_keys)
    for var_ind in 1:length(var_exponents)
        curr_res_str = res_strings[1][var_ind]
        var_exponents[var_ind] = sum([expstr_separate(curr_res_str[res])[2] for res in 1:length(curr_res_str)])
    end
    subspace_str2inds::Vector{Vector{Tuple{Number,Is}}} = []
    for space in statespace.subspaces
        append!(subspace_str2inds, get_op_inds(space, res_strings))
    end
    # Construct terms for each combination of subspace_str2inds
    terms::Vector{QTerm} = []
    coeffs::Vector{Number} = []
    for combo in Iterators.product(subspace_str2inds...)
        curr_inds::Vector{Is} = []
        curr_coeff::Number = 1
        for i in 1:length(combo)
            curr_coeff *= combo[i][1]
            push!(curr_inds, combo[i][2])
        end
        push!(terms, QTerm(curr_inds))
        push!(coeffs, curr_coeff)
    end
    return terms, coeffs, var_exponents
end
function string2qabstract(statespace::StateSpace, operator_str::String)::QAbstract
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
    else 
        if length(expstr) > 0
            subindex = parse(Int, expstr)
        end
    end
    operator_type = statespace.operatortypes[key_index]
    return QAbstract(operator_type, key_index, subindex, exp, conjugate) # subindex only printed if not -1 
end

"""
    term(operator_str::String, statespace::StateSpace)
    term(coeff, operator_str::String, statespace::StateSpace)
    term(coeff::Number, operator_str::String)
    term(operator_str::String)

Generate a quantum term (QTerm) from the StateSpace `q`. The state description is provided as a string.

- Tokens of the form `var^exp` (e.g. `"a^2"`) set the exponent for a state variable.
- Other tokens are assumed to be keys that match one of the allowed subspace keys (i.e. elements in each SubSpace.keys) or 
  are abstract operator names in the StateSpace (with potential indexes, exponents and Daggers ').
- If no coefficient is given, the default coefficient is 1.
- Daggers are given by ', and need to preceed a possible exponent
"""
function term(statespace::StateSpace, coeff::Number, operator_str::String)
    out_strings, out_type = term_pre_split(operator_str, statespace.operator_names)
    var_exponents = zeros(Int, length(statespace.vars))
    coeffs_and_terms::Vector{Vector{Tuple{Number, QAtom}}} = []
    for (out_string, type) in zip(out_strings, out_type)
        if type # QAbstract
            curr_term = string2qabstract(statespace, out_string)
            push!(coeffs_and_terms, [(1, curr_term)])
        else  # QTerm 
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
    products::Vector{QAtomProduct} = []
    for comb in Iterators.product(coeffs_and_terms...)
        curr_coeff = reduce(*, [c[1] for c in comb])*coeff
        curr_atoms = [c[2] for c in comb]
        push!(products, QAtomProduct(statespace, curr_coeff, copy(var_exponents), curr_atoms))
    end
    return QExpr(statespace, products)
end
function term(statespace::StateSpace, operator_str::String)::QExpr
    return term(statespace, one(1), operator_str)
end
function term(operator_str::String)::QExpr
    return term(GLOBAL_STATE_SPACE, operator_str)
end
function term(coeff::Number, operator_str::String)::QExpr
    return term(GLOBAL_STATE_SPACE, coeff, operator_str)
end