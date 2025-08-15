
"""
    base_operators(ss:StateSpace; do_fun::Bool=false, formatted::Bool=false) -> Tuple{Dict{String,QExpr},Dict{String,QExpr},Dict{String,Union{Function, QExpr}}}
    base_operators(statespace::StateSpace, name::String; do_fun::Bool=false, formatted::Bool=false, do_dict::Bool=true) -> Union{Dict{String,QExpr}, Dict{String, Function}, QExpr, Vector{QExpr}}

Returns variables and/or operators in the state space `ss`.
Specifc variables/operators can be selected by passing a string `letter`.
If no `letter` is passed, the function returns a tuple of 3 dictionaries:
- The first dictionary contains the variables in the state space, with their corresponding QExpr objects.
- The second dictionary contains the operators in the state space, with their corresponding QExpr objects.
- The third dictionary contains the abstract operators in the state space either as a callable function to specify the subtype or as a standard QExpr.
If you pass "vars", "ops" or "abstract", it will return a Dictionary with elements for each variable, operator or abstruct operator
    - do_fun specifies if abstract operators are returned as functions, that can be called with no arguments or with an integer to specify the subindex of the abstract operator. 
    - formatted specifies if the dictionary keys are formatted for example true would give: γᵢ, wheras false gives: gamma_i. 
    - do_dict specifies if for non basic (not vars, ops and abstract) elements, we return a dictionary of a Vector of the operators. Requires that the user knows the order of operators for the specific subspace
"""
function base_operators(statespace::StateSpace; do_fun::Bool=false, formatted::Bool=false)::Tuple{Dict{String,QExpr},Dict{String,QExpr}, Dict{String, Union{Function, QExpr}}}
    # return 2 dicctionaries, one with the vars and one with the operators 
    var_dict::Dict{String,QExpr} = Dict()
    op_dict::Dict{String,QExpr} = Dict()
    CRone = one(ComplexRational)
    var_exponents = zeros(Int, length(statespace.vars))
    #abstract_dict::Dict{String,QExpr} = Dict()
    I_operator = statespace.I_op
    for (i, var) in enumerate(statespace.vars)
        var_exponents[i] += 1
        if formatted
            vars_str = var.var_str
        else
            vars_str = var.var_name
        end
        var_dict[vars_str] = QExpr(statespace, QAtomProduct(statespace, CAtom(var_exponents), QTerm[]))
        var_exponents[i] -= 1
    end
    index = 1
    for sub in statespace.subspaces
        op_set = sub.op_set
        base_ops = op_set.base_ops
        for key in sub.keys
            for base_op in base_ops
                curr_operator = copy(I_operator)
                curr_operator[index] = base_op
                term = QTerm(curr_operator)
                curr_name = op_set.op2str(base_op, key, formatted=formatted)
                op_dict[curr_name] = QExpr(statespace, term)
            end

            # non base ops # in a Dict 
            non_base_ops = op_set.non_base_ops
            for (inner_key, ops) in non_base_ops
                curr_terms::Vector{QAtomProduct} = QAtomProduct[]
                for op in ops
                    curr_operator = copy(I_operator)
                    curr_operator[index] = op[2]
                    coeff = op[1]
                    curr_prod = QAtomProduct(statespace, CAtom(coeff, var_exponents), [QTerm(curr_operator)])
                    push!(curr_terms, curr_prod)
                end
                if formatted
                    op_str = inner_key * str2sub(key)
                else
                    op_str = inner_key * "_" * key 
                end 
                op_dict[op_str] = QExpr(statespace, curr_terms)
            end
            index += 1
        end
    end
    op_dict["I"] = QExpr(statespace, QAtomProduct(statespace, CAtom(var_exponents), QTerm[]))
    # now abstract operators 
    abstract_dict::Dict{String, Union{Function, QExpr}} = Dict()
    for (key_index, (name, operatortype)) in enumerate(zip(statespace.operator_names, statespace.operatortypes))
        if do_fun
            abstract_dict[name] = (subindex=-1) -> QExpr(statespace, QAbstract(operatortype, key_index, subindex))
        else
            abstract_dict[name] = QExpr(statespace, QAbstract(operatortype, key_index))
        end
    end
    return var_dict, op_dict, abstract_dict
end

function base_operators(statespace::StateSpace, name::String; do_fun::Bool=false, formatted::Bool=false, do_dict::Bool=true)::Union{Dict{String,QExpr}, Dict{String, Function}, QExpr, Vector{QExpr}, Function}
    # return 2 dicctionaries, one with the vars and one with the operators 
    CRone = one(ComplexRational)
    var_exponents = zeros(Int, length(statespace.vars))
    #abstract_dict::Dict{String,QExpr} = Dict()
    I_operator = statespace.I_op
    if name == "vars"
        var_dict::Dict{String,QExpr} = Dict()
        for (i, var) in enumerate(statespace.vars)
            var_exponents[i] += 1
            if formatted
                vars_str = var.var_str
            else
                vars_str = var.var_name
            end
            var_dict[vars_str] = QExpr(statespace, QAtomProduct(statespace, CAtom(var_exponents), QTerm[]))
            var_exponents[i] -= 1
        end
        return var_dict
    elseif name == "ops"
        op_dict::Dict{String,QExpr} = Dict()
        index = 1
        for sub in statespace.subspaces
            op_set = sub.op_set
            base_ops = op_set.base_ops
            for key in sub.keys
                for base_op in base_ops
                    curr_operator = copy(I_operator)
                    curr_operator[index] = base_op
                    term = QTerm(curr_operator)
                    curr_name = op_set.op2str(base_op, key, formatted=formatted)
                    op_dict[curr_name] = QExpr(statespace, term)
                end

                # non base ops # in a Dict 
                non_base_ops = op_set.non_base_ops
                for (inner_key, ops) in non_base_ops
                    curr_terms::Vector{QAtomProduct} = QAtomProduct[]
                    for op in ops
                        curr_operator = copy(I_operator)
                        curr_operator[index] = op[2]
                        coeff = op[1]
                        curr_prod = QAtomProduct(statespace, CAtom(coeff, var_exponents), [QTerm(curr_operator)])
                        push!(curr_terms, curr_prod)
                    end
                    if formatted
                        op_str = inner_key * str2sub(key)
                    else
                        op_str = inner_key * "_" * key 
                    end 
                    op_dict[op_str] = QExpr(statespace, curr_terms)
                end
                index += 1
            end
        end
        op_dict["I"] = QExpr(statespace, QAtomProduct(statespace, CAtom(var_exponents), QTerm[]))
        return op_dict
    elseif name == "abstract"
        # now abstract operators 
        if do_fun 
            abstract_dict::Dict{String, Function} = Dict()  
            for (key_index, (name, operatortype)) in enumerate(zip(statespace.operator_names, statespace.operatortypes))
                abstract_dict[name] = (subindex=-1) -> QExpr(statespace, QAbstract(operatortype, key_index, subindex))
            end
            return abstract_dict
        else
            abstract_dict2::Dict{String, QExpr} = Dict()
            for (key_index, (name, operatortype)) in enumerate(zip(statespace.operator_names, statespace.operatortypes))
                abstract_dict2[name] = QExpr(statespace, QAbstract(operatortype, key_index))
            end
            return abstract_dict2
        end
    elseif name == "I"
        return QExpr(statespace, QAtomProduct(statespace, CAtom(var_exponents), QTerm[]))
    end

    # vars
    for (i, var) in enumerate(statespace.vars)
        if name == var.var_name || name == var.var_str
            var_exponents[i] += 1
            return QExpr(statespace, QAtomProduct(statespace, CAtom(var_exponents), QTerm[]))
        end
    end
    # Option 2 for vars: gather all the ones for which it occurs in the name, colelct those return if vector is not empty
    curr_ops::Dict{String,QExpr} = Dict()
    ops_vec::Vector{QExpr} = []
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
                curr_ops[vars_str] = QExpr(statespace, QAtomProduct(statespace, CAtom(var_exponents), QTerm[]))
            else 
                push!(ops_vec, QExpr(statespace, QAtomProduct(statespace, CAtom(var_exponents), QTerm[])))
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
                    curr_operator = copy(I_operator)
                    curr_operator[index] = base_op
                    if formatted
                        op_str = inner_key * str2sub(key)
                    else
                        op_str = inner_key * "_" * key 
                    end 
                    if do_dict
                        curr_ops[op_str] = QExpr(statespace, QTerm(curr_operator))
                    else
                        push!(ops_vec, QExpr(statespace, QTerm(curr_operator)))
                    end
                end

                # non base ops # in a Dict 
                non_base_ops = op_set.non_base_ops
                for (inner_key, ops) in non_base_ops
                    curr_terms::Vector{QAtomProduct} = QAtomProduct[]
                    for op in ops
                        curr_operator = copy(I_operator)
                        curr_operator[index] = op[2]
                        coeff = op[1]
                        curr_prod = QAtomProduct(statespace, CAtom(coeff, var_exponents), [QTerm(curr_operator)])
                        push!(curr_terms, curr_prod)
                    end
                    if formatted
                        op_str = inner_key * str2sub(key)
                    else
                        op_str = inner_key * "_" * key 
                    end 
                    if do_dict
                        curr_ops[op_str] = QExpr(statespace, curr_terms)
                    else
                        push!(ops_vec, QExpr(statespace, curr_terms))
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
                return (subindex=-1) -> QExpr(statespace, QAbstract(operatortype, key_index, subindex))
            else
                return QExpr(statespace, QAbstract(operatortype, key_index))
            end
        else
            if contains(name, "_")
                reduced_name = string(split(name, "_")[1])
                if reduced_name == curr_name
                    return QExpr(statespace, string2qabstract(statespace, replace(name, "_" => "")))
                end
            end
        end
    end
    return term(statespace, name)
    # error("No variable, subspace component or abstract operator with key starting with '$letter' found in the state space.")
end