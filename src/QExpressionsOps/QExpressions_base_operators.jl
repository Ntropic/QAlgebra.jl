
"""
    base_operators(statespace::StateSpace, name::String; do_fun::Bool=false) -> Tuples of QExpr, Function, QExprLookup

Returns variables and/or operators in the state space `ss`.
Specifc variables/operators can be selected by passing a string `letter`.
If no `letter` is passed, the function returns a tuple of 3 dictionaries:
- The first dictionary contains the variables in the state space, with their corresponding QExpr objects.
- The second dictionary contains the operators in the state space, with their corresponding QExpr objects.
- The third dictionary contains the abstract operators in the state space either as a callable function to specify the subtype or as a standard QExpr.
If you pass "vars", "ops" or "abstract", it will return a Dictionary with elements for each variable, operator or abstruct operator
    - do_fun specifies if abstract operators are returned as functions, that can be called with no arguments or with an integer to specify the subindex of the abstract operator. 
    - by_ensemble specifies if theoperator subspaces are checked by ensemble key or by sub key.
"""
function base_operators(statespace::StateSpace, name::String; do_fun::Bool=false, by_ensemble::Bool=false)#::Union{QExpr, Vector{QExpr}, Function, QExprLookup}
    # return 2 dicctionaries, one with the vars and one with the operators 
    var_exponents = zeros(Int, length(statespace.vars))
    I_operator = statespace.I_op
    if name == "I"
        return QExpr(statespace, QAtomProduct(statespace, CAtom(var_exponents), QTerm[]))
    end

    # ========>  Variables / Parameters  <========================================================================
    # Direct hit 
    if !do_fun # do first 
        for (i, var) in enumerate(statespace.vars)
            if name == var.var_name || name == var.var_str || name == var.var_name_no_t || name == var.var_str_no_t
                var_exponents[i] += 1
                return QExpr(statespace, QAtomProduct(statespace, CAtom(var_exponents), QTerm[]))
            end
        end
    end
    # Option 2 for vars: gather all the ones for which it occurs in the name, collect those return if vector is not empty
    ops_vec::Vector{QExpr} = []
    ops_comb::Vector{Vector{Symbol}} = []
    for (i, var) in enumerate(statespace.vars)
        pref = string(var.var_symbol)
        pref_formatted = symbol2formatted(pref)
        if name == pref || name == pref_formatted
            #vars_str = var.var_name   # no longer needed -> but can be added to extended do_fun
            var_exponents[i] += 1
            push!(ops_comb, copy(var.index_comb_symbol))
            push!(ops_vec, QExpr(statespace, QAtomProduct(statespace, CAtom(var_exponents), QTerm[])))
            var_exponents[i] -= 1
        end
    end
    if length(ops_vec) > 0
        if length(ops_vec) == 1
            return ops_vec[1]
        end
        if !do_fun
            
            return (ops_vec...,)
        else # make a function that allows indexing variables with Symbols or Strings
            return QExprLookup(ops_comb, ops_vec)
        end
    end
    if !do_fun # do last 
        for (i, var) in enumerate(statespace.vars)
            if name == var.var_name || name == var.var_str || name == var.var_name_no_t || name == var.var_str_no_t
                var_exponents[i] += 1
                return QExpr(statespace, QAtomProduct(statespace, CAtom(var_exponents), QTerm[]))
            end
        end
    end

    # check Operators (subspaces)
    index = 1
    has_underscore = occursin("_", name)
    if has_underscore
        name_inner_key, name = string.(split(name, "_", limit=2))
    end
    for sub in statespace.subspaces
        ensemble_key = false
        if by_ensemble && sub.key == name
            ensemble_key = true
        end
        continue_ = false 
        for key in sub.keys   # key for outer subspace
            do_it = false 
            if ensemble_key 
                do_it = true 
                continue_ = true
            elseif key == name 
                do_it = true 
                continue_ = false 
            end
            if do_it 
                keys = sub.keys
                op_set = sub.op_set
                base_ops = op_set.base_ops
                base_strs = op_set.ops
                for (inner_key, base_op) in zip(base_strs, base_ops)
                    curr_operator = copy(I_operator)
                    curr_operator[index] = base_op
                    if has_underscore && inner_key == name_inner_key    #### Specific key not general key return directly
                        return QExpr(statespace, QTerm(curr_operator))
                    end
                    if by_ensemble
                        push!(ops_comb, [Symbol(inner_key), Symbol(key)])
                    else
                        push!(ops_comb, [Symbol(inner_key)])
                    end
                    push!(ops_vec, QExpr(statespace, QTerm(curr_operator)))
                end
                if !continue_
                    if length(ops_vec) == 1
                        return ops_vec[1]
                    end
                    if !do_fun
                        return (ops_vec...,)
                    else
                        return QExprLookup(ops_comb, ops_vec)
                    end
                end
            end
            index += 1
        end
    end
    if length(ops_vec) > 0 
        if length(ops_vec) == 1
            return ops_vec[1]
        end
        if !do_fun
            return (ops_vec...,)
        else
            return QExprLookup(ops_comb, ops_vec)
        end
    end

    # check for abstract operators
    for (key_index, operatortype) in enumerate(statespace.operatortypes)
        if name == operatortype.name
            if do_fun 
                return (subindex=-1) -> QExpr(statespace, QAbstract(operatortype, key_index, subindex))
            else
                return QExpr(statespace, QAbstract(operatortype, key_index))
            end
        else
            if contains(name, "_")   # specific index 
                reduced_name_s, subindex_s = split(name, "_") 
                reduced_name = string(reduced_name_s)
                subindex = parse(Int, subindex_s) 
                if reduced_name == operatortype.name
                    return QExpr(statespace, QAbstract(operatortype, key_index, subindex)) 
                end
            end
        end
    end
    error("No variable, subspace component or abstract operator found for key='$name'.")
end
function base_operators(statespace::StateSpace, names::Vector{String}; do_fun::Bool=false, by_ensemble::Bool=false)
    return [base_operators(statespace, name; do_fun=do_fun, by_ensemble=by_ensemble) for name in names]
end 