
"""
    base_operators(ss:StateSpace; do_fun::Bool=false, formatted::Bool=false) -> Tuple{Dict{String,qExpr},Dict{String,qExpr},Dict{String,Union{Function, qExpr}}}
    base_operators(statespace::StateSpace, name::String; do_fun::Bool=false, formatted::Bool=false, do_dict::Bool=true) -> Union{Dict{String,qExpr}, Dict{String, Function}, qExpr, Vector{qExpr}}

Returns variables and/or operators in the state space `ss`.
Specifc variables/operators can be selected by passing a string `letter`.
If no `letter` is passed, the function returns a tuple of 3 dictionaries:
- The first dictionary contains the variables in the state space, with their corresponding qExpr objects.
- The second dictionary contains the operators in the state space, with their corresponding qExpr objects.
- The third dictionary contains the abstract operators in the state space either as a callable function to specify the subtype or as a standard qExpr.
If you pass "vars", "ops" or "abstract", it will return a Dictionary with elements for each variable, operator or abstruct operator
    - do_fun specifies if abstract operators are returned as functions, that can be called with no arguments or with an integer to specify the subindex of the abstract operator. 
    - formatted specifies if the dictionary keys are formatted for example true would give: γᵢ, wheras false gives: gamma_i. 
    - do_dict specifies if for non basic (not vars, ops and abstract) elements, we return a dictionary of a Vector of the operators. Requires that the user knows the order of operators for the specific subspace
"""
function base_operators(statespace::StateSpace; do_fun::Bool=false, formatted::Bool=false)::Tuple{Dict{String,qExpr},Dict{String,qExpr}, Dict{String, Union{Function, qExpr}}}
    # return 2 dicctionaries, one with the vars and one with the operators 
    var_dict::Dict{String,qExpr} = Dict()
    op_dict::Dict{String,qExpr} = Dict()
    CRone = one(ComplexRational)
    var_exponents = zeros(Int, length(statespace.vars))
    #abstract_dict::Dict{String,qExpr} = Dict()
    neutral_operator = statespace.neutral_op
    for (i, var) in enumerate(statespace.vars)
        var_exponents[i] += 1
        if formatted
            vars_str = var.var_str
        else
            vars_str = var.var_name
        end
        var_dict[vars_str] = qExpr(statespace, qAtomProduct(statespace, CRone, copy(var_exponents), qTerm[]))
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
                curr_name = op_set.op2str(base_op, key, formatted=formatted)
                op_dict[curr_name] = qExpr(statespace, term)
            end

            # non base ops # in a Dict 
            non_base_ops = op_set.non_base_ops
            for (inner_key, ops) in non_base_ops
                curr_terms::Vector{qAtomProduct} = qAtomProduct[]
                for op in ops
                    curr_operator = copy(neutral_operator)
                    curr_operator[index] = op[2]
                    coeff = op[1]
                    curr_prod = qAtomProduct(statespace, coeff, copy(var_exponents), qTerm(copy(curr_operator)))
                    push!(curr_terms, curr_prod)
                end
                if formatted
                    op_str = inner_key * str2sub(key)
                else
                    op_str = inner_key * "_" * key 
                end 
                op_dict[op_str] = qExpr(statespace, curr_terms)
            end
            index += 1
        end
    end
    op_dict["I"] = qExpr(statespace, qAtomProduct(statespace, CRone, qTerm[]))
    # now abstract operators 
    abstract_dict::Dict{String, Union{Function, qExpr}} = Dict()
    for (key_index, (name, operatortype)) in enumerate(zip(statespace.operator_names, statespace.operatortypes))
        if do_fun
            abstract_dict[name] = (subindex=-1) -> qExpr(statespace, qAbstract(operatortype, key_index, subindex))
        else
            abstract_dict[name] = qExpr(statespace, qAbstract(operatortype, key_index))
        end
    end
    return var_dict, op_dict, abstract_dict
end

function base_operators(statespace::StateSpace, name::String; do_fun::Bool=false, formatted::Bool=false, do_dict::Bool=true)::Union{Dict{String,qExpr}, Dict{String, Function}, qExpr, Vector{qExpr}, Function}
    # return 2 dicctionaries, one with the vars and one with the operators 
    CRone = one(ComplexRational)
    var_exponents = zeros(Int, length(statespace.vars))
    #abstract_dict::Dict{String,qExpr} = Dict()
    neutral_operator = statespace.neutral_op
    if name == "vars"
        var_dict::Dict{String,qExpr} = Dict()
        for (i, var) in enumerate(statespace.vars)
            var_exponents[i] += 1
            if formatted
                vars_str = var.var_str
            else
                vars_str = var.var_name
            end
            var_dict[vars_str] = qExpr(statespace, qAtomProduct(statespace, CRone, copy(var_exponents), qTerm[]))
            var_exponents[i] -= 1
        end
        return var_dict
    elseif name == "ops"
        op_dict::Dict{String,qExpr} = Dict()
        index = 1
        for sub in statespace.subspaces
            op_set = sub.op_set
            base_ops = op_set.base_ops
            for key in sub.keys
                for base_op in base_ops
                    curr_operator = copy(neutral_operator)
                    curr_operator[index] = base_op
                    term = qTerm(curr_operator)
                    curr_name = op_set.op2str(base_op, key, formatted=formatted)
                    op_dict[curr_name] = qExpr(statespace, term)
                end

                # non base ops # in a Dict 
                non_base_ops = op_set.non_base_ops
                for (inner_key, ops) in non_base_ops
                    curr_terms::Vector{qAtomProduct} = qAtomProduct[]
                    for op in ops
                        curr_operator = copy(neutral_operator)
                        curr_operator[index] = op[2]
                        coeff = op[1]
                        curr_prod = qAtomProduct(statespace, coeff, copy(var_exponents), qTerm(copy(curr_operator)))
                        push!(curr_terms, curr_prod)
                    end
                    if formatted
                        op_str = inner_key * str2sub(key)
                    else
                        op_str = inner_key * "_" * key 
                    end 
                    op_dict[op_str] = qExpr(statespace, curr_terms)
                end
                index += 1
            end
        end
        op_dict["I"] = qExpr(statespace, qAtomProduct(statespace, CRone, qTerm[]))
        return op_dict
    elseif name == "abstract"
        # now abstract operators 
        if do_fun 
            abstract_dict::Dict{String, Function} = Dict()  
            for (key_index, (name, operatortype)) in enumerate(zip(statespace.operator_names, statespace.operatortypes))
                abstract_dict[name] = (subindex=-1) -> qExpr(statespace, qAbstract(operatortype, key_index, subindex))
            end
            return abstract_dict
        else
            abstract_dict2::Dict{String, qExpr} = Dict()
            for (key_index, (name, operatortype)) in enumerate(zip(statespace.operator_names, statespace.operatortypes))
                abstract_dict2[name] = qExpr(statespace, qAbstract(operatortype, key_index))
            end
            return abstract_dict2
        end
    elseif name == "I"
        return qExpr(statespace, qAtomProduct(statespace, CRone, qTerm[]))
    end

    # vars
    for (i, var) in enumerate(statespace.vars)
        if name == var.var_name || name == var.var_str
            var_exponents[i] += 1
            return qExpr(statespace, qAtomProduct(statespace, 1, copy(var_exponents), qTerm[]))
        end
    end
    # Option 2 for vars: gather all the ones for which it occurs in the name, colelct those return if vector is not empty
    curr_ops::Dict{String,qExpr} = Dict()
    ops_vec::Vector{qExpr} = []
    for (i, var) in enumerate(statespace.vars)
        pref = string(split(var.var_name, "_")[1])
        pref_subs = pref
        if haskey(var_substitution, pref)
            pref_subs = var_substitution[pref]
        end
        if name == pref || name == pref_subs
            if formatted
                vars_str = var.var_str
            else
                vars_str = var.var_name
            end
            var_exponents[i] += 1
            if do_dict
                curr_ops[vars_str] = qExpr(statespace, qAtomProduct(statespace, 1, copy(var_exponents), qTerm[]))
            else 
                push!(ops_vec, qExpr(statespace, qAtomProduct(statespace, 1, copy(var_exponents), qTerm[])))
            end
            var_exponents[i] -= 1
        end
    end
    if length(curr_ops) > 0 
        return curr_ops
    end
    if length(ops_vec) > 0
        return ops_vec
    end

    # check Operators (subspaces)
    index = 1
    for sub in statespace.subspaces
        for key in sub.keys
            if key == name
                keys = sub.keys
                op_set = sub.op_set
                base_ops = op_set.base_ops
                base_strs = op_set.ops
                for (inner_key, base_op) in zip(base_strs, base_ops)
                    curr_operator = copy(neutral_operator)
                    curr_operator[index] = base_op
                    if formatted
                        op_str = inner_key * str2sub(key)
                    else
                        op_str = inner_key * "_" * key 
                    end 
                    if do_dict
                        curr_ops[op_str] = qExpr(statespace, qTerm(copy(curr_operator)))
                    else
                        push!(ops_vec, qExpr(statespace, qTerm(copy(curr_operator))))
                    end
                end

                # non base ops # in a Dict 
                non_base_ops = op_set.non_base_ops
                for (inner_key, ops) in non_base_ops
                    curr_terms::Vector{qAtomProduct} = qAtomProduct[]
                    for op in ops
                        curr_operator = copy(neutral_operator)
                        curr_operator[index] = op[2]
                        coeff = op[1]
                        curr_prod = qAtomProduct(statespace,coeff, copy(var_exponents), qTerm(copy(curr_operator)))
                        push!(curr_terms, curr_prod)
                    end
                    if formatted
                        op_str = inner_key * str2sub(key)
                    else
                        op_str = inner_key * "_" * key 
                    end 
                    if do_dict
                        curr_ops[op_str] = qExpr(statespace, curr_terms)
                    else
                        push!(ops_vec, qExpr(statespace, curr_terms))
                    end
                end
                if do_dict 
                    return curr_ops
                else 
                    return ops_vec
                end
            end
            index += 1
        end
    end
    # check for abstract operators
    for (key_index, (curr_name, operatortype)) in enumerate(zip(statespace.operator_names, statespace.operatortypes))
        if name == curr_name
            if do_fun 
                return (subindex=-1) -> qExpr(statespace, qAbstract(operatortype, key_index, subindex))
            else
                return qExpr(statespace, qAbstract(operatortype, key_index))
            end
        else
            if contains(name, "_")
                reduced_name = string(split(name, "_")[1])
                if reduced_name == curr_name
                    return qExpr(statespace, string2qabstract(statespace, replace(name, "_" => "")))
                end
            end
        end
    end
    return term(statespace, name)
    # error("No variable, subspace component or abstract operator with key starting with '$letter' found in the state space.")
end