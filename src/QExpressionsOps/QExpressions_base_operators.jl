
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
function base_operators(statespace::StateSpace, name::String; do_fun::Bool=false, by_ensemble::Bool=false)
    # -------------------- helpers --------------------
    _I_expr() = QExpr(statespace, QAtomProduct(statespace, CAtom(zeros(Int, length(statespace.vars))), QTerm[]))
    _qexpr_for_var(i::Int) = begin
        vexp = zeros(Int, length(statespace.vars))
        vexp[i] += 1
        QExpr(statespace, QAtomProduct(statespace, CAtom(vexp), QTerm[]))
    end
    # parse "x_tN"
    function _parse_time_suffix(s::String)
        m = match(r"^(.*)_t(-?\d+)$", s)
        m === nothing ? (s, nothing) : (String(m.captures[1]), parse(Int, m.captures[2]))
    end

    # -------------------- identity --------------------
    if name == "I"
        return _I_expr()
    end

    # ==================> Variables / Parameters <===============================
    name_base, t_spec = _parse_time_suffix(name)

    ops_comb = Vector{Vector{Symbol}}()
    ops_vec  = QExpr[]

    for (i, var) in enumerate(statespace.vars)
        pref = string(var.var_symbol)
        pref_fmt = symbol2formatted(pref)

        if name_base == var.var_name_no_t || name_base == var.var_str_no_t ||
           name_base == pref || name_base == pref_fmt

            if t_spec === nothing || t_spec == var.t_index
                comb = copy(var.index_comb_symbol)
                push!(comb, Symbol("t$(var.t_index)"))
                push!(ops_comb, comb)
                push!(ops_vec, _qexpr_for_var(i))
            end
        end
    end

    if !isempty(ops_vec)
        if do_fun
            return QExprLookup(ops_comb, ops_vec)
        else
            return length(ops_vec) == 1 ? ops_vec[1] : (ops_vec...,)
        end
    end

    # ===================> Subspace Operators <=================================
    index = 1
    has_us = occursin("_", name)
    name_inner_key = nothing
    outer_name = name
    if has_us
        name_inner_key, outer_name = string.(split(name, "_", limit=2))
    end

    ops_comb = Vector{Vector{Symbol}}()
    ops_vec  = QExpr[]

    for sub in statespace.subspaces
        ensemble_key = by_ensemble && sub.key == outer_name
        for key in sub.keys
            do_it = ensemble_key || (key == outer_name)
            if do_it
                base_ops = sub.op_set.base_ops
                base_strs = sub.op_set.ops
                for (inner_key, base_op) in zip(base_strs, base_ops)
                    curr = copy(statespace.I_op)
                    curr[index] = base_op
                    if has_us && inner_key == name_inner_key
                        return QExpr(statespace, QTerm(curr))
                    end
                    if by_ensemble
                        push!(ops_comb, [Symbol(inner_key), Symbol(key)])
                    else
                        push!(ops_comb, [Symbol(inner_key)])
                    end
                    push!(ops_vec, QExpr(statespace, QTerm(curr)))
                end
                if !ensemble_key
                    return do_fun ? QExprLookup(ops_comb, ops_vec) :
                           (length(ops_vec) == 1 ? ops_vec[1] : (ops_vec...,))
                end
            end
            index += 1
        end
    end

    if !isempty(ops_vec)
        return do_fun ? QExprLookup(ops_comb, ops_vec) :
               (length(ops_vec) == 1 ? ops_vec[1] : (ops_vec...,))
    end

    # ===================> Abstract Operators <=================================
    for (key_index, operatortype) in enumerate(statespace.operatortypes)
        if name == operatortype.name
            return do_fun ?
                ((subindex=-1) -> QExpr(statespace, QAbstract(operatortype, key_index, subindex))) :
                QExpr(statespace, QAbstract(operatortype, key_index))
        elseif contains(name, "_")
            reduced_name_s, subindex_s = split(name, "_")
            reduced_name = String(reduced_name_s)
            subindex = parse(Int, subindex_s)
            if reduced_name == operatortype.name
                return QExpr(statespace, QAbstract(operatortype, key_index, subindex))
            end
        end
    end

    error("No variable, subspace component or abstract operator found for key='$name'.")
end

function base_operators(statespace::StateSpace, names::Vector{String}; do_fun::Bool=false, by_ensemble::Bool=false)
    return [base_operators(statespace, name; do_fun=do_fun, by_ensemble=by_ensemble) for name in names]
end 