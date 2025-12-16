@inline function _identity_expr(qspace::QSpace)
    atom = CAtom(qspace.param_info, spzeros(Int, length(qspace.params)))
    return QExpr(qspace, QAtomProduct(qspace, atom, QTerm[]))
end

@inline function _qexpr_for_var(qspace::QSpace, index::Int)
    exponents = spzeros(Int, length(qspace.params))
    exponents[index] += 1
    atom = CAtom(qspace.param_info, exponents)
    return QExpr(qspace, QAtomProduct(qspace, atom, QTerm[]))
end

function _parse_time_suffix(name::String)
    m = match(r"^(.*)_t(-?\d+)$", name)
    return m === nothing ? (name, nothing) : (String(m.captures[1]), parse(Int, m.captures[2]))
end

@inline _default_time_index(of_time::Bool, t_spec::Union{Nothing,Int}) =
    of_time ? (t_spec === nothing ? 0 : t_spec) : -1

@inline function _check_t_bounds(t::Int, max_t::Int)
    (t < 0 || t > max_t) && error("time_index=$t out of range 0:$max_t")
end

function _finalize_lookup_result(ops_vec::Vector{QExpr}, ops_comb::Vector{Vector{Symbol}}, do_fun::Bool)
    isempty(ops_vec) && return nothing
    if do_fun
        return QExprLookup(ops_comb, ops_vec)
    else
        return length(ops_vec) == 1 ? ops_vec[1] : (ops_vec...,)
    end
end

function get_parameter(qspace::QSpace, name::String;
                       do_fun::Bool=false,
                       name_base::String=name,
                       t_spec::Union{Nothing,Int}=nothing,
                       error_if_missing::Bool=true)
    ops_comb = Vector{Vector{Symbol}}()
    ops_vec = QExpr[]

    for (i, var) in enumerate(qspace.params)
        pref = string(var.param_symbol)
        pref_fmt = symbol2formatted(pref)

        if name_base == var.param_name_no_t || name_base == var.param_str_no_t ||
           name_base == pref || name_base == pref_fmt

            if t_spec === nothing || t_spec == var.t_index
                comb = copy(var.index_comb_symbol)
                push!(comb, Symbol("t$(var.t_index)"))
                push!(ops_comb, comb)
                push!(ops_vec, _qexpr_for_var(qspace, i))
            end
        end
    end

    result = _finalize_lookup_result(ops_vec, ops_comb, do_fun)
    result !== nothing && return result

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

    if error_if_missing
        error("No parameter found for key='$name'.")
    end
    return nothing
end

function get_operator(qspace::QSpace, name::String;
                      do_fun::Bool=false,
                      by_ensemble::Bool=false,
                      name_base::String=name,
                      t_spec::Union{Nothing,Int}=nothing,
                      error_if_missing::Bool=true)
    if name_base == "I"
        return _identity_expr(qspace)
    end
    index = 1
    has_us = occursin("_", name_base)
    name_inner_key = nothing
    outer_name = name_base
    if has_us
        name_inner_key, outer_name = string.(split(name_base, "_", limit=2))
    end

    ops_comb = Vector{Vector{Symbol}}()
    ops_vec = QExpr[]

    sub_of_time = qspace.subspace_info.of_time
    max_t = qspace.max_t_ind
    ti = _default_time_index(sub_of_time, t_spec)

    for sub in qspace.subspaces
        ensemble_key = by_ensemble && sub.key == outer_name
        match_outer = sub.key == outer_name
        for (inner_idx, key_symbol) in enumerate(subspace_symbols(sub))
            key = String(key_symbol)
            do_it = ensemble_key || (key == outer_name) || (match_outer && !has_us && !by_ensemble)
            if do_it
                base_ops = sub.op_set.base_ops
                base_strs = sub.op_set.ops

                for (inner_key, base_op) in zip(base_strs, base_ops)
                    curr = copy(qspace.I_op)
                    curr[index] = base_op

                    base_comb = by_ensemble ? [Symbol(inner_key), Symbol(key)] : [Symbol(inner_key)]

                    if has_us && inner_key == name_inner_key
                        if do_fun && sub_of_time && t_spec === nothing
                            q0 = QExpr(qspace, QTerm(curr, 0))
                            push!(ops_comb, base_comb); push!(ops_vec, q0)
                            for t in 0:max_t
                                qt = t == 0 ? q0 : QExpr(qspace, QTerm(curr, t))
                                comb_t = copy(base_comb); push!(comb_t, Symbol("t$(t)"))
                                push!(ops_comb, comb_t); push!(ops_vec, qt)
                            end
                            return QExprLookup(ops_comb, ops_vec)
                        else
                            if sub_of_time && t_spec !== nothing
                                _check_t_bounds(t_spec, max_t)
                            end
                            return QExpr(qspace, QTerm(curr, ti))
                        end
                    end

                    if sub_of_time
                        if t_spec === nothing
                            if do_fun
                                q0 = QExpr(qspace, QTerm(curr, 0))
                                push!(ops_comb, base_comb); push!(ops_vec, q0)
                                for t in 0:max_t
                                    qt = t == 0 ? q0 : QExpr(qspace, QTerm(curr, t))
                                    comb_t = copy(base_comb); push!(comb_t, Symbol("t$(t)"))
                                    push!(ops_comb, comb_t); push!(ops_vec, qt)
                                end
                            else
                                push!(ops_vec, QExpr(qspace, QTerm(curr, 0)))
                            end
                        else
                            _check_t_bounds(t_spec, max_t)
                            if do_fun
                                comb_t = copy(base_comb); push!(comb_t, Symbol("t$(t_spec)"))
                                push!(ops_comb, comb_t)
                            end
                            push!(ops_vec, QExpr(qspace, QTerm(curr, t_spec)))
                        end
                    else
                        push!(ops_comb, base_comb)
                        push!(ops_vec, QExpr(qspace, QTerm(curr, -1)))
                    end
                end

                if !ensemble_key
                    result = _finalize_lookup_result(ops_vec, ops_comb, do_fun)
                    result !== nothing && return result
                end
            end
            index += 1
        end
    end

    result = _finalize_lookup_result(ops_vec, ops_comb, do_fun)
    result !== nothing && return result

    if error_if_missing
        error("No subspace operator found for key='$name'.")
    end
    return nothing
end

function get_abstract(qspace::QSpace, name::String;
                      do_fun::Bool=false,
                      name_base::String=name,
                      t_spec::Union{Nothing,Int}=nothing,
                      error_if_missing::Bool=true)
    for (key_index, operatortype) in enumerate(qspace.operatortypes)
        ti_default = _default_time_index(operatortype.of_time, t_spec)
        max_t = qspace.max_t_ind

        if name_base == operatortype.name
            if do_fun
                if operatortype.of_time
                    return (subindex::Int=-1, time_index::Int=0) -> begin
                        _check_t_bounds(time_index, max_t)
                        QExpr(qspace, QAbstract(operatortype, key_index, subindex, 1, false, time_index))
                    end
                else
                    return (subindex::Int=-1) -> begin
                        QExpr(qspace, QAbstract(operatortype, key_index, subindex, 1, false, -1))
                    end
                end
            else
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

    if error_if_missing
        error("No abstract operator found for key='$name'.")
    end
    return nothing
end

"""
    base_operators(qspace::QSpace, name::String; do_fun::Bool=false, by_ensemble::Bool=false)

Returns variables and/or operators in the state space `qspace`.

Time handling:
- If an operator (subspace or abstract) is not time-dependent (`.of_time == false`), its `time_index` defaults to `-1`.
- If it is time-dependent, `time_index` defaults to `0`, unless `name` ends with `_tN` (e.g. `"A_t3"`), in which case `time_index = N`.
- For QExprLookup results of subspace operators, the time argument is optional: you can query with or without a `:tN` key; both map to the same expression.
"""
function base_operators(qspace::QSpace, name::String; do_fun::Bool=false, by_ensemble::Bool=false)
    name_base, t_spec = _parse_time_suffix(name)
    if name_base == "I"
        return _identity_expr(qspace)
    end

    result = get_parameter(qspace, name; do_fun=do_fun, name_base=name_base, t_spec=t_spec, error_if_missing=false)
    result !== nothing && return result

    result = get_operator(qspace, name; do_fun=do_fun, by_ensemble=by_ensemble, name_base=name_base, t_spec=t_spec, error_if_missing=false)
    result !== nothing && return result

    result = get_abstract(qspace, name; do_fun=do_fun, name_base=name_base, t_spec=t_spec, error_if_missing=false)
    result !== nothing && return result

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

function Base.getindex(qspace::QSpace, key::Union{String,Symbol})
    name = key isa Symbol ? String(key) : key
    name_base, t_spec = _parse_time_suffix(name)
    if name_base == "I"
        return _identity_expr(qspace)
    end

    result = get_parameter(qspace, name; name_base=name_base, t_spec=t_spec, error_if_missing=false)
    result !== nothing && return result

    result = get_operator(qspace, name; name_base=name_base, t_spec=t_spec, error_if_missing=false)
    result !== nothing && return result

    result = get_abstract(qspace, name; name_base=name_base, t_spec=t_spec, error_if_missing=false)
    result !== nothing && return result

    error("No variable, subspace component or abstract operator found for key='$name'.")
end

function Base.getindex(qspace::QSpace, keys::AbstractVector{<:Union{String,Symbol}})
    return [qspace[key] for key in keys]
end
