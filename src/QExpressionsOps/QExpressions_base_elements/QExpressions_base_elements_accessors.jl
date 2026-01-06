using ..QIndexes: AbstractIndex, TimeIndex
using ..StringUtils: normalize_underscore_indices, split_index
using ..QSpaces: CR_ONE

"""
    AbstractOperatorAccessor

Accessor returned for abstract operators with unspecified indices.

Usage supports `()` and `[]`:
- Index: `Int` or numeric `String`/`Symbol` (e.g., `1`, `"2"`, `:3`).
- Time (if applicable): `Int`, `String`, `Symbol`, or `TimeIndex` (e.g., `2`, `"t2"`, `:t2`, `TimeIndex(2)`).

Example mixing types:
`A(1, :t2)` or `A["2", 0]`.
"""
struct AbstractOperatorAccessor
    qspace::QSpace
    operator_type::OperatorType
    key_index::Int
    fixed_time::Union{Nothing, Int}
end

"""
    ParameterGroupAccessor

Accessor returned for parameter groups with unspecified indices.

Usage supports `()` and `[]`:
- Indices: `Int`, `String`, or `Symbol` (e.g., `1`, `"i2"`, `:j3`).
- Time (if applicable): `Int`, `String`, `Symbol`, or `TimeIndex` (e.g., `2`, `"t2"`, `:t2`, `TimeIndex(2)`).

Example mixing types:
`gamma(1, "j2", :t3)` or `gamma[:i1, 2, "t0"]`.
"""
struct ParameterGroupAccessor
    qspace::QSpace
    group::ParameterGroupLike
    group_index::Int
    subspace_indices::Vector{Int}
    ensemble_indices::Vector{Int}
    index_name_pairs::Vector{Tuple{String,String}}
    fixed_time::Union{Nothing, Int}
end

function Base.show(io::IO, accessor::AbstractOperatorAccessor)
    label = accessor.operator_type.name
    has_time = accessor.operator_type.of_time
    fixed = accessor.fixed_time === nothing ? "" : " (fixed t_$(accessor.fixed_time))"
    if has_time
        print(io, "AbstractOperatorAccessor: call as ", label, "(n) or ", label, "(n, t)",
              " {or index as ", label, "[n] or ", label, "[n, t]}.", fixed)
    else
        print(io, "AbstractOperatorAccessor: call as ", label, "(n)",
              " {or index as ", label, "[n]}.", fixed)
    end
end

function Base.show(io::IO, accessor::ParameterGroupAccessor)
    label = string(accessor.group.param_symbol)
    n_indices = length(accessor.group.indices)
    idx_hint = n_indices == 0 ? "" : (n_indices == 1 ? "i" : "i1,...,i$(n_indices)")
    call_args = n_indices == 0 ? "" : idx_hint
    accessor.group.of_t && (call_args = isempty(call_args) ? "t" : call_args * ", t")
    idx_args = n_indices == 0 ? "" : idx_hint
    accessor.group.of_t && (idx_args = isempty(idx_args) ? "t" : idx_args * ", t")
    call_usage = label * "(" * call_args * ")"
    idx_usage = label * "[" * idx_args * "]"
    fixed = accessor.fixed_time === nothing ? "" : " (fixed t_$(accessor.fixed_time))"
    print(io, "ParameterGroupAccessor: call as ", call_usage,
          " {or index as ", idx_usage, "}.", fixed)
end

@inline function _parse_time_arg(arg)::Int
    if arg isa TimeIndex
        return arg.order
    elseif arg isa Integer
        return Int(arg)
    end
    base, comps = normalize_underscore_indices(String(arg))
    base == "t" || error("Time argument must use symbol t, got $(arg).")
    if isempty(comps)
        return 0
    elseif length(comps) == 1
        return parse(Int, comps[1])
    else
        error("Time argument can only have a single index, got $(arg).")
    end
end

@inline function _resolve_time(fixed::Union{Nothing, Int}, provided::Union{Nothing, Int}, of_time::Bool, label::String)::Int
    provided !== nothing && fixed !== nothing && provided != fixed &&
        error("Time index for $(label) already fixed to t_$(fixed), got t_$(provided).")
    if of_time
        return provided !== nothing ? provided : (fixed !== nothing ? fixed : 0)
    else
        provided === nothing || error("$(label) is not time dependent, but a time index was provided.")
        return -1
    end
end

@inline function _parse_param_index(subspace::Int, ensemble::Int, names::Tuple{String,String}, arg)::AbstractIndex
    base = ""
    number = 0
    if arg isa Integer
        base, number = "", Int(arg)
    else
        base, number = split_index(String(arg))
    end
    summation = false
    if isempty(base)
        summation = false
    elseif base == names[1]
        summation = false
    elseif base == names[2]
        summation = true
    else
        error("Index $(arg) does not match expected labels $(names).")
    end
    return AbstractIndex(subspace, ensemble, number, summation)
end

@inline function _resolve_op_index(op_set::OperatorSet, arg)::Int
    if arg isa Integer
        idx = Int(arg)
        1 <= idx <= length(op_set.ops) || error("Operator index $(idx) out of bounds, available operators: $(op_set.ops).")
        return idx
    end
    name = arg isa Symbol ? String(arg) : (arg isa AbstractString ? String(arg) :
        error("Operator name must be a String, Symbol, or Int, got $(arg)."))
    idx = findfirst(==(name), op_set.ops)
    idx === nothing && error("Unrecognized operator $(name), subspace supports: $(op_set.ops).")
    return idx
end

@inline function _abstract_index_from_arg(arg)::Int
    if arg isa Integer
        return Int(arg)
    end
    str = arg isa Symbol ? String(arg) : (arg isa AbstractString ? String(arg) :
        error("Abstract operator index must be numeric, got $(arg)."))
    try
        return parse(Int, str)
    catch
        error("Abstract operator index must be numeric, got $(arg).")
    end
end

function (accessor::ParameterGroupAccessor)(params::Vararg{Any})::QExpr
    group = accessor.group
    n_indices = length(group.indices)
    n_params = length(params)
    idx_args = params
    provided_time = nothing

    if n_params == n_indices + 1
        provided_time = _parse_time_arg(params[end])
        idx_args = params[1:n_indices]
    elseif n_params != n_indices
        suffix = group.of_t ? " (and optional time)" : ""
        error("Expected $(n_indices) indices$(suffix), got $(n_params).")
    end

    time_idx = _resolve_time(accessor.fixed_time, provided_time, group.of_t, string(group.param_symbol))

    abstract_indices = Vector{AbstractIndex}(undef, n_indices)
    @inbounds for pos in 1:n_indices
        abstract_indices[pos] = _parse_param_index(accessor.subspace_indices[pos],
                                                   accessor.ensemble_indices[pos],
                                                   accessor.index_name_pairs[pos],
                                                   idx_args[pos])
    end

    particle = CParticle(accessor.group_index, 1, abstract_indices, TimeIndex(time_idx))
    atom = CAtom(accessor.qspace.param_info, CR_ONE, CParticle{AbstractIndex}[particle])
    return QExpr(accessor.qspace, QComposite[QAtomProduct(accessor.qspace, atom, QAtom[])])
end

function Base.getindex(accessor::ParameterGroupAccessor, params::Vararg{Any})::QExpr
    return accessor(params...)
end

function (accessor::AbstractOperatorAccessor)(params::Vararg{Any})::QExpr
    length(params) in (1, 2) ||
        error("Expected index or (index, time) for abstract operator $(accessor.operator_type.name), got $(length(params)) arguments.")

    sub_index = _abstract_index_from_arg(params[1])
    provided_time = length(params) == 2 ? _parse_time_arg(params[2]) : nothing
    time_idx = _resolve_time(accessor.fixed_time, provided_time, accessor.operator_type.of_time, accessor.operator_type.name)

    abstrac_op = QAbstract(accessor.operator_type, accessor.key_index, sub_index, 1, false, time_idx)
    return QExpr(accessor.qspace, QComposite[QAtomProduct(accessor.qspace, accessor.qspace.c_one, QAtom[abstrac_op])])
end

function Base.getindex(accessor::AbstractOperatorAccessor, params::Vararg{Any})::QExpr
    return accessor(params...)
end

"""
    SubSpaceAccessor

Accessor returned for subspaces with unspecified operators.

Usage supports `()` and `[]`:
- Operator: `String`, `Symbol`, or `Int` (operator name or position).
- Index (ensemble only): `Int`, `String`, or `Symbol` (e.g., `1`, `"i2"`, `:j3`).
- Time (if `qspace.of_time`): `Int`, `String`, `Symbol`, or `TimeIndex` (e.g., `2`, `"t2"`, `:t2`, `TimeIndex(2)`).

Example mixing types:
`i(:x, "j2", :t0)` or `i["y", :i1, 1]`.
"""
struct SubSpaceAccessor
    qspace::QSpace
    subspace::SubSpace
    subspace_index::Int
    ensemble_index::Int
    index_name_pair::Tuple{String,String}
end

"""
    SubSpaceOperatorsAccessor

Accessor returned for subspace operator sets with unspecified ensemble indices.

Usage supports `()` and `[]`:
- Index (ensemble only): `Int`, `String`, or `Symbol` (e.g., `1`, `"i2"`, `:j3`).
- Time (if `qspace.of_time`): `Int`, `String`, `Symbol`, or `TimeIndex` (e.g., `2`, `"t2"`, `:t2`, `TimeIndex(2)`).

Example mixing types:
`ops("i1", "t0")` or `ops[:j3, 1]`.
"""
struct SubSpaceOperatorsAccessor
    qspace::QSpace
    subspace::SubSpace
    subspace_index::Int
    ensemble_index::Int
    index_name_pair::Tuple{String,String}
end

function Base.show(io::IO, accessor::SubSpaceAccessor)
    label = accessor.subspace.key
    time_hint = accessor.qspace.of_time ? ", t" : ""
    if accessor.subspace.is_ensemble_ss
        print(io, "SubSpaceAccessor: call as ", label, "(op, i", time_hint, ")",
              " {or index as ", label, "[op, i", time_hint, "]}.")
    else
        print(io, "SubSpaceAccessor: call as ", label, "(op", time_hint, ")",
              " {or index as ", label, "[op", time_hint, "]}.")
    end
end

function Base.show(io::IO, accessor::SubSpaceOperatorsAccessor)
    label = accessor.subspace.key
    time_hint = accessor.qspace.of_time ? ", t" : ""
    print(io, "SubSpaceOperatorsAccessor: call as ", label, "(i", time_hint, ")",
          " {or index as ", label, "[i", time_hint, "]}.")
end

@inline function _subspace_ops_exprs(qspace::QSpace, sub::SubSpace, abstract_index::AbstractIndex, time_index::TimeIndex)::Vector{QExpr}
    op_set = sub.op_set
    ops = Vector{QExpr}(undef, length(op_set.base_ops))
    @inbounds for idx in eachindex(op_set.base_ops)
        operator = QParticle(op_set.base_ops[idx], abstract_index)
        term = QTerm(QParticle[operator], time_index)
        ops[idx] = QExpr(qspace, QComposite[QAtomProduct(qspace, qspace.c_one, QAtom[term])], Val(:nosimp))
    end
    return ops
end

function (accessor::SubSpaceAccessor)(params::Vararg{Any})::QExpr
    sub = accessor.subspace
    n_params = length(params)
    n_params == 0 && error("Expected operator argument for subspace $(sub.key).")
    provided_time = nothing
    if accessor.qspace.of_time
        if sub.is_ensemble_ss
            if n_params == 3
                provided_time = _parse_time_arg(params[3])
                n_params = 2
            elseif n_params != 2
                error("Expected (op, index) or (op, index, t) for ensemble subspace $(sub.key), got $(length(params)) arguments.")
            end
        else
            if n_params == 2
                provided_time = _parse_time_arg(params[2])
                n_params = 1
            elseif n_params != 1
                error("Expected (op) or (op, t) for subspace $(sub.key), got $(length(params)) arguments.")
            end
        end
    else
        if sub.is_ensemble_ss
            n_params == 2 ||
                error("Expected (op, index) for ensemble subspace $(sub.key), got $(length(params)) arguments.")
        else
            n_params == 1 ||
                error("Expected (op) for subspace $(sub.key), got $(length(params)) arguments.")
        end
    end

    op_idx = _resolve_op_index(sub.op_set, params[1])
    abstract_index = if sub.is_ensemble_ss
        _parse_param_index(accessor.subspace_index, accessor.ensemble_index, accessor.index_name_pair, params[2])
    else
        AbstractIndex(accessor.subspace_index, accessor.ensemble_index, 0, false)
    end
    operator = QParticle(sub.op_set.base_ops[op_idx], abstract_index)
    term = QTerm(QParticle[operator], resolve_time_index(accessor.qspace, provided_time))
    return QExpr(accessor.qspace, QComposite[QAtomProduct(accessor.qspace, accessor.qspace.c_one, QAtom[term])], Val(:nosimp))
end

function Base.getindex(accessor::SubSpaceAccessor, params::Vararg{Any})::QExpr
    return accessor(params...)
end

function (accessor::SubSpaceOperatorsAccessor)(params::Vararg{Any})::Vector{QExpr}
    sub = accessor.subspace
    sub.is_ensemble_ss ||
        error("SubSpaceOperatorsAccessor only applies to ensemble subspaces, got $(sub.key).")
    n_params = length(params)
    provided_time = nothing
    if accessor.qspace.of_time
        if n_params == 2
            provided_time = _parse_time_arg(params[2])
        elseif n_params != 1
            error("Expected (index) or (index, t) for ensemble subspace $(sub.key), got $(length(params)) arguments.")
        end
    else
        n_params == 1 ||
            error("Expected (index) for ensemble subspace $(sub.key), got $(length(params)) arguments.")
    end
    abstract_index = _parse_param_index(accessor.subspace_index, accessor.ensemble_index,
                                        accessor.index_name_pair, params[1])
    return _subspace_ops_exprs(accessor.qspace, sub, abstract_index,
                               resolve_time_index(accessor.qspace, provided_time))
end

function Base.getindex(accessor::SubSpaceOperatorsAccessor, params::Vararg{Any})::Vector{QExpr}
    return accessor(params...)
end
