"""
    base_operators(qspace::QSpace, name::String; do_fun::Bool=false, by_ensemble::Bool=false)

Returns variables and/or operators in the state space `qspace`.

Time handling:
- If an operator (subspace or abstract) is not time-dependent (`.of_time == false`), its `time_index` defaults to `-1`.
- If it is time-dependent, `time_index` defaults to `0`, unless `name` ends with `_tN` (e.g. `"A_t3"`), in which case `time_index = N`.
- For QExprLookup results of subspace operators, the time argument is optional: you can query with or without a `:tN` key; both map to the same expression.
"""
function base_operators(qspace::QSpace, name::String; do_fun::Bool=false, by_ensemble::Bool=false)
    # -------------------- helpers --------------------
    _I_expr() = QExpr(qspace, QAtomProduct(qspace, CAtom(qspace.param_info, spzeros(Int, length(qspace.params))), QTerm[]))

    _qexpr_for_var(i::Int) = begin
        vexp = spzeros(Int, length(qspace.params))
        vexp[i] += 1
        QExpr(qspace, QAtomProduct(qspace, CAtom(qspace.param_info, vexp), QTerm[]))
    end

    # parse "x_tN" once and reuse
    function _parse_time_suffix(s::String)
        m = match(r"^(.*)_t(-?\d+)$", s)
        m === nothing ? (s, nothing) : (String(m.captures[1]), parse(Int, m.captures[2]))
    end

    _default_time_index(of_time::Bool, t_spec::Union{Nothing, Int}) =
        of_time ? (t_spec === nothing ? 0 : t_spec) : -1

    _check_t_bounds(t::Int, max_t::Int) = (t < 0 || t > max_t) && error("time_index=$t out of range 0:$max_t")

    # -------------------- identity --------------------
    name_base, t_spec = _parse_time_suffix(name)
    if name_base == "I"
        return _I_expr()
    end

    # ==================> Variables / Parameters <===============================
    ops_comb = Vector{Vector{Symbol}}()
    ops_vec  = QExpr[]

    for (i, var) in enumerate(qspace.params)
        pref     = string(var.param_symbol)
        pref_fmt = symbol2formatted(pref)

        if name_base == var.param_name_no_t || name_base == var.param_str_no_t ||
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

    # ==================> CAbstract Parameters <===============================
    if name_base == name 
        abstract_definitions = qspace.param_info.abstract_definitions
        for abstract_defintion in abstract_definitions
            if abstract_defintion.name == name
                cfun = CAbstract(qspace.param_info, ComplexRational(1,0,1), abstract_defintion.index)
                return QExpr(qspace, [QAtomProduct(qspace, cfun)])
            end
        end
    end
    if name == "CAbstract"
        abstract_definitions = qspace.param_info.abstract_definitions
        exprs = QExpr[]
        for abstract_defintion in abstract_definitions
            if abstract_defintion.name == name
                cfun = CAbstract(qspace.param_info, ComplexRational(1,0,1), abstract_defintion.index)
                push!(exprs, QExpr(qspace, [QAtomProduct(qspace, cfun)]))
            end
        end
        return exprs
    end

    # ===================> Subspace Operators <=================================
    index = 1
    has_us = occursin("_", name_base)
    name_inner_key = nothing
    outer_name = name_base
    if has_us
        name_inner_key, outer_name = string.(split(name_base, "_", limit=2))
    end

    ops_comb = Vector{Vector{Symbol}}()
    ops_vec  = QExpr[]

    sub_of_time = qspace.subspace_info.of_time
    max_t = qspace.max_t_ind
    ti = _default_time_index(sub_of_time, t_spec)  # used for non-iterating cases

    for sub in qspace.subspaces
        ensemble_key = by_ensemble && sub.key == outer_name
        for key_symbol in sub.keys_symbols
            key = String(key_symbol)
            do_it = ensemble_key || (key == outer_name)
            if do_it
                base_ops  = sub.op_set.base_ops
                base_strs = sub.op_set.ops

                for (inner_key, base_op) in zip(base_strs, base_ops)
                    curr = copy(qspace.I_op)
                    curr[index] = base_op

                    base_comb = by_ensemble ? [Symbol(inner_key), Symbol(key)] : [Symbol(inner_key)]

                    # If an exact inner op was requested (e.g. "a_b"), special-case:
                    if has_us && inner_key == name_inner_key
                        if do_fun && sub_of_time && t_spec === nothing
                            # Build a lookup over all times + unsuffixed->t0
                            # unsuffixed key -> t=0
                            q0 = QExpr(qspace, QTerm(curr, 0))
                            push!(ops_comb, base_comb); push!(ops_vec, q0)
                            # all explicit times 0..max_t
                            for t in 0:max_t
                                qt = t == 0 ? q0 : QExpr(qspace, QTerm(curr, t))
                                comb_t = copy(base_comb); push!(comb_t, Symbol("t$(t)"))
                                push!(ops_comb, comb_t); push!(ops_vec, qt)
                            end
                            return QExprLookup(ops_comb, ops_vec)
                        else
                            # single QExpr path (either timeless, or time specified)
                            if sub_of_time && t_spec !== nothing
                                _check_t_bounds(t_spec, max_t)
                            end
                            return QExpr(qspace, QTerm(curr, ti))
                        end
                    end

                    # General accumulation path
                    if sub_of_time
                        if t_spec === nothing
                            if do_fun
                                # Add unsuffixed -> t0
                                q0 = QExpr(qspace, QTerm(curr, 0))
                                push!(ops_comb, base_comb); push!(ops_vec, q0)
                                # Add all explicit times 0..max_t
                                for t in 0:max_t
                                    qt = t == 0 ? q0 : QExpr(qspace, QTerm(curr, t))
                                    comb_t = copy(base_comb); push!(comb_t, Symbol("t$(t)"))
                                    push!(ops_comb, comb_t); push!(ops_vec, qt)
                                end
                            else
                                # no lookup requested: default to t=0
                                push!(ops_vec, QExpr(qspace, QTerm(curr, 0)))
                            end
                        else
                            _check_t_bounds(t_spec, max_t)
                            if do_fun
                                comb_t = copy(base_comb); push!(comb_t, Symbol("t$(t_spec)"))
                                push!(ops_comb, comb_t)
                            else
                                # ignore combs when not returning a lookup
                            end
                            push!(ops_vec, QExpr(qspace, QTerm(curr, t_spec)))
                        end
                    else
                        # timeless
                        push!(ops_comb, base_comb)  # harmless if do_fun=false
                        push!(ops_vec, QExpr(qspace, QTerm(curr, -1)))
                    end
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

    # ===================> QAbstract Operators <=================================
    for (key_index, operatortype) in enumerate(qspace.operatortypes)
        ti_default = _default_time_index(operatortype.of_time, t_spec)
        max_t = qspace.max_t_ind

        if name_base == operatortype.name
            if do_fun
                if operatortype.of_time
                    # time arg present with default 0; enforce 0..max_t
                    return (subindex::Int=-1, time_index::Int=0) -> begin
                        _check_t_bounds(time_index, max_t)
                        QExpr(qspace, QAbstract(operatortype, key_index, subindex, 1, false, time_index))
                    end
                else
                    # no time arg when not of_time
                    return (subindex::Int=-1) -> begin
                        QExpr(qspace, QAbstract(operatortype, key_index, subindex, 1, false, -1))
                    end
                end
            else
                # honor optional _tN in name; if provided and of_time, also bound-check it
                if operatortype.of_time && t_spec !== nothing
                    _check_t_bounds(t_spec, max_t)
                end
                return QExpr(qspace, QAbstract(operatortype, key_index, -1, 1, false, ti_default))
            end
        elseif occursin("_", name_base)
            reduced_name_s, subindex_s = split(name_base, "_", limit=2)
            reduced_name = String(reduced_name_s)
            subindex = parse(Int, subindex_s)
            if reduced_name == operatortype.name
                if operatortype.of_time && t_spec !== nothing
                    _check_t_bounds(t_spec, max_t)
                end
                return QExpr(qspace, QAbstract(operatortype, key_index, subindex, 1, false, ti_default))
            end
        end
    end

    error("No variable, subspace component or abstract operator found for key='$name'.")
end

function base_operators(qspace::QSpace, name::Symbol; do_fun::Bool=false, by_ensemble::Bool=false)
    return base_operators(qspace, String(name); do_fun=do_fun, by_ensemble=by_ensemble)
end

function base_operators(qspace::QSpace, names::Vector{String}; do_fun::Bool=false, by_ensemble::Bool=false)
    return [base_operators(qspace, name; do_fun=do_fun, by_ensemble=by_ensemble) for name in names]
end 
function base_operators(qspace::QSpace, names::Vector{Symbol}; do_fun::Bool=false, by_ensemble::Bool=false)
    return [base_operators(qspace, name; do_fun=do_fun, by_ensemble=by_ensemble) for name in names]
end 
