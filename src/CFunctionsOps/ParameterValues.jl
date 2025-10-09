const _GroupStorage = Union{ComplexF64, Array{ComplexF64}, Vector{Float64}}

export attach_samples!, register_ensemble_sample_size!, resolve_param!

using Base: WeakRef
using ..Sampler: QEnsembleFunction, QDistribution, build_discrete_samples, build_continuous_samples
using ..EnsembleSamples: AbstractEnsembleSample
using ..ParameterGroups: ParameterGroup, ParameterGroupKind,
                          ParameterGroupScalar, ParameterGroupTimeFunction,
                          ParameterGroupDistribution, ParameterGroupEnsembleFunction

"""
    ParameterValues(param_info::ParameterInfo; qspace_ref=WeakRef())

Concrete value store that mirrors the parameter layout described by
`param_info`. Each parameter group owns either a scalar, a time vector, or a
dense array over `(time, indexes...)`. Payloads declared on the shared
`ParameterGroup` structures (distributions, scalar functions, ensemble
functions, literals) drive the storage layout and evaluation behaviour. The
constructor allocates storage and marks which groups are ready based on those
payloads.

`ParameterValues` keeps lightweight bookkeeping so function-driven groups update
in a deterministic order: `time_group` points to the explicit time parameter and
`update_order` captures the remaining function groups sorted by their declared
dependencies. Scalar groups update first via [`update_functions`](@ref), followed
by ensemble-backed groups handled by [`update_ensemble_group_functions`](@ref).
"""
mutable struct ParameterValues
    qspace::WeakRef
    param_info::ParameterInfo
    group_values::Vector{_GroupStorage}
    group_definition_initialized::BitVector
    group_initialized::BitVector
    group_time_initialized::Vector{Vector{Bool}}
    got_all_definitions::Bool
    update_order::Vector{Int}
    time_group::Int
end

function ParameterValues(param_info::ParameterInfo; qspace_ref::WeakRef=WeakRef())
    groups = param_info.param_groups
    group_count = length(groups)

    group_values = Vector{_GroupStorage}(undef, group_count)
    group_definition_initialized = falses(group_count)
    group_initialized = falses(group_count)
    group_time_initialized = Vector{Vector{Bool}}(undef, group_count)

    @inbounds for g in 1:group_count
        group = groups[g]
        time_count = max(group.time_count, 1)
        group_time_initialized[g] = fill(false, time_count)
        index_sizes = _effective_index_sizes(group)
        if group.kind == ParameterGroupDistribution
            count = length(group.parameter_indices)
            group_values[g] = Vector{Float64}(undef, count)
            group_definition_initialized[g] = group.payload isa QDistribution
            continue
        end

        if isempty(index_sizes) && time_count == 1
            group_values[g] = ComplexF64(NaN)
        else
            dims = isempty(index_sizes) ? (time_count,) : (time_count, index_sizes...)
            group_values[g] = Array{ComplexF64}(undef, dims...)
        end

        if group.kind == ParameterGroupTimeFunction
            group_definition_initialized[g] = group.payload isa Function
        elseif group.kind == ParameterGroupEnsembleFunction
            group_definition_initialized[g] = group.payload isa QEnsembleFunction
        elseif group.kind == ParameterGroupScalar && group.payload !== nothing
            group_definition_initialized[g] = true
        end

        if group.is_time_group
            group_definition_initialized[g] = true
        end
    end

    got_all_definitions = all(group_definition_initialized)
    time_group_idx = findfirst(group -> group.is_time_group, groups)
    time_group = time_group_idx === nothing ? 0 : time_group_idx
    update_order = _compute_update_order(groups, time_group)

    pv = ParameterValues(qspace_ref, param_info, group_values, group_definition_initialized,
                         group_initialized, group_time_initialized, got_all_definitions, update_order, time_group)

    @inbounds for g in 1:group_count
        group = groups[g]
        if group.kind == ParameterGroupScalar && group.payload !== nothing
            _set_group!(pv, g, group.payload; allow_function=true)
        end
    end

    recompute_functions!(pv, 1)

    return pv
end

function _compute_update_order(groups::Vector{ParameterGroup}, time_group::Int)
    function_groups = Int[]
    ensemble_groups = Int[]
    @inbounds for (idx, group) in enumerate(groups)
        idx == time_group && continue
        if group.kind == ParameterGroupTimeFunction
            push!(function_groups, idx)
        elseif group.kind == ParameterGroupEnsembleFunction
            push!(ensemble_groups, idx)
        end
    end
    return isempty(function_groups) ?
        ensemble_groups :
        isempty(ensemble_groups) ? function_groups : vcat(function_groups, ensemble_groups)
end

@inline function _effective_index_sizes(group::ParameterGroup)
    sizes = group.sample_sizes
    if !isempty(sizes) && length(sizes) == length(group.indexes) && all(>(0), sizes)
        return sizes
    end
    return group.index_sizes
end

@inline function _coords_tuple(coords::AbstractVector{Int})
    return Tuple(coords)
end

@inline function _reset_group_flags!(pv::ParameterValues, group_idx::Int)
    flags = pv.group_time_initialized[group_idx]
    fill!(flags, false)
    pv.group_initialized[group_idx] = false
end

@inline function _mark_slot_initialized!(pv::ParameterValues, group_idx::Int, slot::Int)
    flags = pv.group_time_initialized[group_idx]
    slot = clamp(slot, 1, length(flags))
    flags[slot] = true
    pv.group_initialized[group_idx] = flags[1]
end

@inline function _clamp_slot(pv::ParameterValues, group_idx::Int, slot::Int)
    flags = pv.group_time_initialized[group_idx]
    return clamp(slot, 1, length(flags))
end

@inline function _ensure_group_slot!(pv::ParameterValues, group_idx::Int, slot::Int)
    slot_idx = _clamp_slot(pv, group_idx, slot)
    flags = pv.group_time_initialized[group_idx]
    flags[slot_idx] && return
    group = pv.param_info.param_groups[group_idx]
    if group.kind == ParameterGroupTimeFunction
        _evaluate_scalar_function_group!(pv, group_idx, slot_idx)
    elseif group.kind == ParameterGroupEnsembleFunction
        _evaluate_qensemble_group!(pv, group_idx, slot_idx)
    else
        flags[slot_idx] = true
        pv.group_initialized[group_idx] = flags[1]
    end
end

@inline function _distribution_slot(info::ParameterInfo, group_idx::Int, param_idx::Int)
    params = info.param_groups[group_idx].parameter_indices
    pos = findfirst(==(param_idx), params)
    pos === nothing && error("Parameter $(info.params_str[param_idx]) does not belong to group $(info.param_groups[group_idx].name).")
    return pos
end

@inline function _distribution_slot(info::ParameterInfo, group_idx::Int, coords::Vector{Int})
    params = info.param_groups[group_idx].parameter_indices
    for (pos, idx) in enumerate(params)
        info.param_coords[idx] == coords && return pos
    end
    error("No parameter with coordinates $(coords) in group $(info.param_groups[group_idx].name).")
end

function _store_value!(pv::ParameterValues, param_idx::Int, value; allow_function::Bool=false)
    info = pv.param_info
    group_idx = info.param_group_by_index[param_idx]
    group = info.param_groups[group_idx]
    if !allow_function && ((group.kind == ParameterGroupTimeFunction && group.payload !== nothing) ||
                            (group.kind == ParameterGroupEnsembleFunction && group.payload !== nothing))
        error("Cannot assign values to function-defined parameter group $(group.name).")
    end
    coords = info.param_coords[param_idx]
    storage = pv.group_values[group_idx]
    if storage isa ComplexF64
        length(coords) == 1 ||
            error("Expected scalar storage for parameter $(info.params_name[param_idx]).")
        pv.group_values[group_idx] = ComplexF64(value)
    elseif storage isa Array{ComplexF64}
        storage_tuple = _coords_tuple(coords)
        storage[storage_tuple...] = ComplexF64(value)
    elseif storage isa Vector{Float64}
        position = _distribution_slot(info, group_idx, param_idx)
        storage[position] = Float64(value)
    else
        error("Unsupported storage type $(typeof(storage)) for group $(group_idx).")
    end
    slot_idx = length(coords) >= 1 ? coords[1] : 1
    _mark_slot_initialized!(pv, group_idx, slot_idx)
    pv.group_definition_initialized[group_idx] = true
    pv.got_all_definitions = all(pv.group_definition_initialized)
    return value
end

function get_parameter_index(info::ParameterInfo, dicts::ParameterDicts, name::Symbol)
    matches = get(dicts.param_name_to_indices, name, nothing)
    matches === nothing && error("No parameter named $(name) registered in ParameterInfo.")
    length(matches) == 1 && return matches[1]
    labels = info.params_str[matches]
    error("Parameter name $(name) is ambiguous. Matches: $(join(labels, ", ")).")
end

function get_parameter_index(info::ParameterInfo, name::Symbol)
    dicts = build_parameter_dicts(info)
    return get_parameter_index(info, dicts, name)
end

get_parameter_index(info::ParameterInfo, dicts::ParameterDicts, name::String) =
    get_parameter_index(info, dicts, Symbol(name))

get_parameter_index(info::ParameterInfo, name::String) = get_parameter_index(info, Symbol(name))
get_parameter_index(pv::ParameterValues, name) = get_parameter_index(pv.param_info, name)

function _fill_group!(pv::ParameterValues, group_idx::Int, value)
    storage = pv.group_values[group_idx]
    if storage isa ComplexF64
        pv.group_values[group_idx] = ComplexF64(value)
    elseif storage isa Array{ComplexF64}
        storage .= ComplexF64(value)
    elseif storage isa Vector{Float64}
        storage .= Float64(value)
    else
        error("Unsupported storage type $(typeof(storage)) for group $(group_idx).")
    end
    flags = pv.group_time_initialized[group_idx]
    fill!(flags, true)
    pv.group_initialized[group_idx] = flags[1]
    return value
end

function _set_group_array!(pv::ParameterValues, group_idx::Int, value::AbstractArray)
    info = pv.param_info
    group = info.param_groups[group_idx]
    time_count = max(group.time_count, 1)
    index_sizes = group.index_sizes
    expected_dims = isempty(index_sizes) ? (time_count,) : (time_count, index_sizes...)
    storage = pv.group_values[group_idx]
    if storage isa ComplexF64
        length(value) == 1 ||
            error("Value shape $(size(value)) does not match expected scalar for group $(group.name).")
        pv.group_values[group_idx] = ComplexF64(value[1])
    elseif storage isa Array{ComplexF64}
        isempty(index_sizes) && time_count == 1 && return _fill_group!(pv, group_idx, value[1])
        size(value) == expected_dims ||
            error("Value shape $(size(value)) does not match expected $(expected_dims) for group $(group.name).")
        storage .= ComplexF64.(value)
    elseif storage isa Vector{Float64}
        length(value) == length(storage) ||
            error("Value length $(length(value)) does not match expected $(length(storage)) for group $(group.name).")
        storage .= Float64.(value)
    else
        error("Group $(group.name) does not accept array assignments.")
    end
    flags = pv.group_time_initialized[group_idx]
    fill!(flags, true)
    pv.group_initialized[group_idx] = flags[1]
    return value
end

function set_param!(pv::ParameterValues, param_idx::Int, value)
    _store_value!(pv, param_idx, value)
    return value
end

function _set_group!(pv::ParameterValues, group_idx::Int, value; allow_function::Bool=false)
    info = pv.param_info
    group = info.param_groups[group_idx]
    if !allow_function && ((group.kind == ParameterGroupTimeFunction && group.payload !== nothing) ||
                           (group.kind == ParameterGroupEnsembleFunction && group.payload !== nothing))
        error("Cannot assign values to function-defined parameter group $(group.name).")
    end
    if value isa Number
        _fill_group!(pv, group_idx, value)
    elseif value isa AbstractArray
        _set_group_array!(pv, group_idx, value)
    else
        error("Unsupported value type $(typeof(value)) for group assignment.")
    end
    pv.group_definition_initialized[group_idx] = true
    pv.got_all_definitions = all(pv.group_definition_initialized)
    return value
end

function set_param!(pv::ParameterValues, name::Symbol, value)
    info = pv.param_info
    group_idx = findfirst(==(name), info.outer_labels_symbols)
    if group_idx !== nothing
        return _set_group!(pv, group_idx, value)
    end
    param_idx = get_parameter_index(info, name)
    return set_param!(pv, param_idx, value)
end

set_param!(pv::ParameterValues, name::String, value) = set_param!(pv, Symbol(name), value)

function update_t!(pv::ParameterValues, value, slot::Int=0)
    info = pv.param_info
    idx = nothing
    for (param_idx, is_t) in enumerate(info.param_is_t)
        is_t || continue
        coords = info.param_coords[param_idx]
        if coords[1] - 1 == slot
            idx = param_idx
            break
        end
    end
    idx === nothing && error("No time parameter t$(slot) registered in ParameterValues.")
    set_param!(pv, idx, value)
    slot_idx = slot + 1
    pv.got_all_definitions && recompute_functions!(pv, slot_idx)
    return value
end

set_time!(pv::ParameterValues, value) = update_t!(pv, value, 0)
set_time!(pv::ParameterValues, value, slot::Int) = update_t!(pv, value, slot)

function _lookup_storage_value(info::ParameterInfo, storage, group_idx::Int, coords::Vector{Int})
    if storage isa ComplexF64
        return storage
    elseif storage isa Array{ComplexF64}
        group = info.param_groups[group_idx]
        idxs = if !isempty(group.sample_sizes) && ndims(storage) == length(coords) - 1
            _coords_tuple(@view coords[2:end])
        else
            _coords_tuple(coords)
        end
        return @inbounds storage[idxs...]
    elseif storage isa Vector{Float64}
        position = _distribution_slot(info, group_idx, coords)
        return @inbounds storage[position]
    end
    error("Unsupported storage type $(typeof(storage)) for group $(group_idx).")
end

@inline function _raw_value(pv::ParameterValues, param_idx::Int)
    info = pv.param_info
    group_idx = info.param_group_by_index[param_idx]
    storage = pv.group_values[group_idx]
    if storage isa Vector{Float64}
        position = _distribution_slot(info, group_idx, param_idx)
        return @inbounds storage[position]
    end
    coords = info.param_coords[param_idx]
    slot_idx = length(coords) >= 1 ? coords[1] : 1
    slot_idx = _clamp_slot(pv, group_idx, slot_idx)
    _ensure_group_slot!(pv, group_idx, slot_idx)
    coords_adj = copy(coords)
    if !isempty(coords_adj)
        coords_adj[1] = slot_idx
    end
    return _lookup_storage_value(info, storage, group_idx, coords_adj)
end

function value(pv::ParameterValues, param_idx::Int)
    coords = pv.param_info.param_coords[param_idx]
    slot_idx = length(coords) >= 1 ? coords[1] : 1
    recompute_functions!(pv, slot_idx)
    return _raw_value(pv, param_idx)
end

function _resolve_index_coords(info::ParameterInfo, param_idx::Int, indexes::ConcreteIndexes)
    tuples = info.param_index_tuples[param_idx]
    coords = copy(info.param_coords[param_idx])
    @inbounds for (offset, (ensemble, inner)) in enumerate(tuples)
        coords[offset + 1] = indexes.indexes[ensemble][inner]
    end
    return coords
end

function value(pv::ParameterValues, param_idx::Int, indexes::ConcreteIndexes)
    info = pv.param_info
    coords = _resolve_index_coords(info, param_idx, indexes)
    slot_idx = length(coords) >= 1 ? coords[1] : 1
    recompute_functions!(pv, slot_idx)
    group_idx = info.param_group_by_index[param_idx]
    storage = pv.group_values[group_idx]
    if storage isa Vector{Float64}
        position = _distribution_slot(info, group_idx, coords)
        return @inbounds storage[position]
    end
    slot_idx = _clamp_slot(pv, group_idx, slot_idx)
    coords_adj = copy(coords)
    if !isempty(coords_adj)
        coords_adj[1] = slot_idx
    end
    return _lookup_storage_value(info, storage, group_idx, coords_adj)
end

function value(pv::ParameterValues, name::Symbol, indexes::Union{Nothing,ConcreteIndexes}=nothing)
    idx = get_parameter_index(pv, name)
    if indexes === nothing
        return value(pv, idx)
    else
        return value(pv, idx, indexes)
    end
end

value(pv::ParameterValues, name::String, indexes::Union{Nothing,ConcreteIndexes}=nothing) =
    value(pv, Symbol(name), indexes)

function Base.show(io::IO, pv::ParameterValues)
    info = pv.param_info
    groups = info.param_groups
    group_count = length(groups)
    if get(io, :compact, false)
        print(io, "ParameterValues(", group_count, " groups)")
        return
    end
    println(io, "ParameterValues:")
    for idx in 1:group_count
        status = pv.group_initialized[idx] ? "✓" : "x"
        group = groups[idx]
        base_name = group.display_signature
        pdf_hint = group.kind == ParameterGroupDistribution ? " (pdf)" : ""
        storage = pv.group_values[idx]
        size_str = if storage isa ComplexF64
            "1"
        elseif storage isa Vector{Float64}
            string(length(storage))
        elseif storage isa Array{ComplexF64}
            dims = size(storage)
            isempty(dims) ? "1" : join(string.(dims), "×")
        else
            string(typeof(storage))
        end
        println(io, "  ", status, " ", base_name, pdf_hint, " (size=", size_str, ")")
    end
end

@inline function _evaluate_scalar_function_group!(pv::ParameterValues, group_idx::Int, slot::Int)
    info = pv.param_info
    group = info.param_groups[group_idx]
    func = group.payload
    func isa Function || return
    slot_idx = _clamp_slot(pv, group_idx, slot)
    for idx in group.parameter_indices
        coords = info.param_coords[idx]
        if group.of_t && coords[1] != slot_idx
            continue
        end
        refs = info.function_param_refs[idx]
        if refs === nothing
            result = func()
        else
            args = map(refs) do ref
                dep_group = info.param_group_by_index[ref]
                dep_coords = info.param_coords[ref]
                _ensure_group_slot!(pv, dep_group, dep_coords[1])
                _raw_value(pv, ref)
            end
            result = func(args...)
        end
        _store_value!(pv, idx, result; allow_function=true)
    end
    if group.of_t
        _mark_slot_initialized!(pv, group_idx, slot_idx)
    else
        flags = pv.group_time_initialized[group_idx]
        fill!(flags, true)
        pv.group_initialized[group_idx] = true
    end
end

# Fallback evaluator for ensemble groups without full sample data.
@inline function _evaluate_qensemble_group_paramwise!(pv::ParameterValues, group_idx::Int, slot::Int)
    info = pv.param_info
    group = info.param_groups[group_idx]
    func = group.payload
    func isa QEnsembleFunction || return
    slot_idx = _clamp_slot(pv, group_idx, slot)
    for idx in group.parameter_indices
        refs = info.function_param_refs[idx]
        refs === nothing && continue
        coords = info.param_coords[idx]
        if group.of_t && coords[1] != slot_idx
            continue
        end
        args = map(refs) do ref
            dep_group = info.param_group_by_index[ref]
            dep_coords = info.param_coords[ref]
            _ensure_group_slot!(pv, dep_group, dep_coords[1])
            _raw_value(pv, ref)
        end
        result = func.func(args...)
        _store_value!(pv, idx, result; allow_function=true)
    end
    if group.of_t
        _mark_slot_initialized!(pv, group_idx, slot_idx)
    else
        flags = pv.group_time_initialized[group_idx]
        fill!(flags, true)
        pv.group_initialized[group_idx] = true
    end
end

@inline function _time_value(pv::ParameterValues, time_group_idx::Int, t_slot::Int)
    if time_group_idx == 0
        return t_slot - 1
    end
    slot_idx = _clamp_slot(pv, time_group_idx, t_slot)
    _ensure_group_slot!(pv, time_group_idx, slot_idx)
    storage = pv.group_values[time_group_idx]
    if storage isa AbstractArray
        return storage[slot_idx]
    elseif storage isa AbstractVector
        return storage[slot_idx]
    else
        return storage
    end
end

@inline function _extract_storage_value(storage, indices::Tuple)
    if storage isa AbstractArray
        return isempty(indices) ? storage[] : storage[indices...]
    elseif storage isa AbstractVector
        isempty(indices) && error("Vector storage expects at least one index.")
        length(indices) == 1 || error("Vector storage expects one index, got $(length(indices)).")
        return storage[first(indices)]
    else
        return storage
    end
end

@inline function _argument_value(pv::ParameterValues,
                                 arg_group_idx::Int,
                                 t_slot::Int,
                                 sample_coords::Vector{Int},
                                 self_positions::Vector{Int})
    if arg_group_idx == 0
        return t_slot - 1
    end
    info = pv.param_info
    group = info.param_groups[arg_group_idx]
    slot_idx = _clamp_slot(pv, arg_group_idx, t_slot)
    _ensure_group_slot!(pv, arg_group_idx, slot_idx)
    storage = pv.group_values[arg_group_idx]
    time_indices = group.of_t ? (slot_idx,) : ()
    if isempty(group.indexes)
        return _extract_storage_value(storage, time_indices)
    end
    sub_indices = isempty(self_positions) ? Int[] : [sample_coords[pos] for pos in self_positions]
    length(sub_indices) == length(group.indexes) ||
        error("Mismatch between argument sample positions and target group indexes for group $(group.name).")
    indices_tuple = group.of_t ? (time_indices[1], sub_indices...) : Tuple(sub_indices)
    return _extract_storage_value(storage, indices_tuple)
end

# Strip spurious zero-imaginary parts before passing user payloads.
@inline function _normalize_argument(val)
    if val isa Complex
        imag(val) == 0 ? real(val) : val
    else
        val
    end
end

# Call QEnsembleFunction safely, retrying with truncated arguments if needed.
@inline function _call_ensemble_function(func::QEnsembleFunction, args::Vector)
    try
        return func.func(args...)
    catch err
        if err isa MethodError && err.f === func.func
            arities = map(m -> m.nargs - 1, methods(func.func))
            if !isempty(arities)
                target = maximum(arities)
                if target < length(args)
                    return func.func(args[1:target]...)
                end
            end
        end
        rethrow(err)
    end
end

# Compute an ensemble-function group when sample grids are available.
@inline function _evaluate_qensemble_group_samples!(pv::ParameterValues, group_idx::Int, slot::Int)
    info = pv.param_info
    group = info.param_groups[group_idx]
    func = group.payload::QEnsembleFunction
    sample_sizes = group.sample_sizes
    length(sample_sizes) == length(group.indexes) && all(>(0), sample_sizes) || return false
    slot_idx = _clamp_slot(pv, group_idx, slot)
    for (arg_pos, dep_group_idx) in pairs(func.argument_group_indices)
        dep_group_idx == 0 && continue
        dep_group = info.param_groups[dep_group_idx]
        dep_slot = dep_group.of_t ? slot_idx : 1
        _ensure_group_slot!(pv, dep_group_idx, dep_slot)
    end
    _ensure_ensemble_storage!(pv, group_idx)
    storage = pv.group_values[group_idx]
    time_count = max(group.time_count, 1)
    storage isa AbstractArray{ComplexF64} || return false
    size(storage, 1) == time_count || return false
    args = Vector{Any}(undef, length(func.argument_symbols))
    sample_axis_count = length(sample_sizes)
    if sample_axis_count == 0
        sample_coords = Int[]
        @inbounds for (arg_pos, arg_name) in pairs(func.argument_group_names)
            dep_group_idx = func.argument_group_indices[arg_pos]
            dep_group = info.param_groups[dep_group_idx]
            dep_slot = dep_group.of_t ? slot_idx : 1
            if arg_name == "t"
                args[arg_pos] = _normalize_argument(_time_value(pv, dep_group_idx, slot_idx))
            else
                arg_val = _argument_value(pv, dep_group_idx, dep_slot, sample_coords, func.argument_self_index_positions[arg_pos])
                args[arg_pos] = _normalize_argument(arg_val)
            end
        end
        result = _call_ensemble_function(func, args)
        storage[slot_idx] = ComplexF64(result)
        _mark_slot_initialized!(pv, group_idx, slot_idx)
        return true
    end
    dims = Tuple(sample_sizes)
    size(storage) == (time_count, dims...) || return false
    sample_coords = Vector{Int}(undef, sample_axis_count)
    for CI in CartesianIndices(dims)
        @inbounds for k in 1:sample_axis_count
            sample_coords[k] = CI[k]
        end
        @inbounds for (arg_pos, arg_name) in pairs(func.argument_group_names)
            dep_group_idx = func.argument_group_indices[arg_pos]
            dep_group = info.param_groups[dep_group_idx]
            dep_slot = dep_group.of_t ? slot_idx : 1
            if arg_name == "t"
                args[arg_pos] = _normalize_argument(_time_value(pv, dep_group_idx, slot_idx))
            else
                arg_val = _argument_value(pv, dep_group_idx, dep_slot, sample_coords, func.argument_self_index_positions[arg_pos])
                args[arg_pos] = _normalize_argument(arg_val)
            end
        end
        result = _call_ensemble_function(func, args)
        idxs = (slot_idx, sample_coords...)
        storage[idxs...] = ComplexF64(result)
    end
    _mark_slot_initialized!(pv, group_idx, slot_idx)
    return true
end

@inline function _evaluate_qensemble_group!(pv::ParameterValues, group_idx::Int, slot::Int)
    info = pv.param_info
    group = info.param_groups[group_idx]
    func = group.payload
    func isa QEnsembleFunction || return
    _evaluate_qensemble_group_samples!(pv, group_idx, slot) && return
    _evaluate_qensemble_group_paramwise!(pv, group_idx, slot)
end

# Execute scalar time and function groups in dependency order.
function update_time_group(pv::ParameterValues, slot::Int=1)
    group_count = length(pv.group_values)
    time_group = pv.time_group
    (0 < time_group <= group_count) || return pv
    group = pv.param_info.param_groups[time_group]
    slot_idx = _clamp_slot(pv, time_group, slot)
    if group.kind == ParameterGroupTimeFunction && group.payload isa Function
        _evaluate_scalar_function_group!(pv, time_group, slot_idx)
    elseif group.kind == ParameterGroupEnsembleFunction && group.payload isa QEnsembleFunction
        _evaluate_qensemble_group!(pv, time_group, slot_idx)
    end
    return pv
end

# Recompute all pure function groups (no ensemble coupling).
function update_functions(pv::ParameterValues, slot::Int=1)
    isempty(pv.update_order) && return pv
    for g in pv.update_order
        g == pv.time_group && continue
        group = pv.param_info.param_groups[g]
        group.kind == ParameterGroupTimeFunction || continue
        group.payload isa Function || continue
        slot_idx = _clamp_slot(pv, g, slot)
        _evaluate_scalar_function_group!(pv, g, slot_idx)
    end
    return pv
end

# Recompute ensemble-driven groups, respecting dependency readiness.
function update_ensemble_group_functions(pv::ParameterValues, slot::Int=1)
    isempty(pv.update_order) && return pv
    for g in pv.update_order
        group = pv.param_info.param_groups[g]
        group.kind == ParameterGroupEnsembleFunction || continue
        group.payload isa QEnsembleFunction || continue
        slot_idx = _clamp_slot(pv, g, slot)
        _evaluate_qensemble_group!(pv, g, slot_idx)
    end
    return pv
end

"""
    recompute_functions!(pv::ParameterValues)

Evaluate all function-defined parameter groups attached to `pv`, respecting time
dependencies and ensemble ordering. Scalar groups driven by pure functions are
updated before time-dependent ensembles so that downstream evaluations see the
latest values.
"""
# Ensure every function-backed group reflects the latest inputs.
function recompute_functions!(pv::ParameterValues, slot::Int=1)
    update_time_group(pv, slot)
    update_functions(pv, slot)
    update_ensemble_group_functions(pv, slot)
    return pv
end

ensure_functions!(pv::ParameterValues) = recompute_functions!(pv)

function _assign_group_storage!(storage, values::AbstractVector{<:Real})
    if storage isa Vector{Float64}
        resize!(storage, length(values))
        storage .= Float64.(values)
        return storage
    elseif storage isa Array{ComplexF64}
        storage .= ComplexF64.(values)
        return storage
    elseif storage isa ComplexF64
        return ComplexF64(first(values))
    else
        return Float64.(values)
    end
end

# Resize ensemble storage to match the available sample counts.
function _ensure_ensemble_storage!(pv::ParameterValues, group_idx::Int)
    info = pv.param_info
    group = info.param_groups[group_idx]
    sizes = group.sample_sizes
    isempty(sizes) && return
    length(sizes) == length(group.indexes) || return
    all(>(0), sizes) || return
    time_count = max(group.time_count, 1)
    dims = isempty(group.indexes) ? (time_count,) : (time_count, sizes...)
    storage = pv.group_values[group_idx]
    if storage isa Array{ComplexF64} && size(storage) == dims
        return
    end
    pv.group_values[group_idx] = Array{ComplexF64}(undef, dims)
    _reset_group_flags!(pv, group_idx)
end

# Record sample counts for each ensemble index of a group.
function _update_group_sample_sizes!(pv::ParameterValues, group_idx::Int, positions::Vector{Int}, sample_count::Int)
    info = pv.param_info
    group = info.param_groups[group_idx]
    isempty(group.indexes) && return
    isempty(positions) && return
    n = length(group.indexes)
    if isempty(group.sample_sizes) || length(group.sample_sizes) != n
        group.sample_sizes = fill(0, n)
    end
    for pos in positions
        1 <= pos <= n || error("Sample size position $(pos) out of bounds for group $(group.name).")
        group.sample_sizes[pos] = sample_count
    end
    if group.kind == ParameterGroupEnsembleFunction && all(>(0), group.sample_sizes)
        _ensure_ensemble_storage!(pv, group_idx)
    end
    _reset_group_flags!(pv, group_idx)
end

# Broadcast a new ensemble sample size to every dependent group.
function register_ensemble_sample_size!(pv::ParameterValues, outer_idx::Int, sample_count::Int)
    info = pv.param_info
    for (group_idx, group) in enumerate(info.param_groups)
        positions = findall(==(outer_idx), group.index_outer_subspaces)
        isempty(positions) && continue
        _update_group_sample_sizes!(pv, group_idx, positions, sample_count)
    end
    return pv
end

# Assign raw sample vectors to their owning groups.
function attach_samples!(pv::ParameterValues, sample::AbstractEnsembleSample)
    info = pv.param_info
    for (col, group_idx) in enumerate(sample.group_indices)
        values = sample.samples[:, col]
        group = info.param_groups[group_idx]
        if !isempty(group.indexes)
            positions = collect(eachindex(group.indexes))
            _update_group_sample_sizes!(pv, group_idx, positions, length(values))
        else
            _reset_group_flags!(pv, group_idx)
        end
        storage = pv.group_values[group_idx]
        pv.group_values[group_idx] = _assign_group_storage!(storage, values)
        flags = pv.group_time_initialized[group_idx]
        fill!(flags, true)
        pv.group_initialized[group_idx] = flags[1]
        pv.group_definition_initialized[group_idx] = true
    end
    pv.got_all_definitions = all(pv.group_definition_initialized)
    return sample
end

"""
    resolve_param!(qspace, group_name, payload)

Attach or update the payload for the parameter group identified by
`group_name` on an existing `qspace`.  The helper validates the payload type
against the group's declared kind, updates the shared `ParameterGroup` record,
marks the corresponding storage in `qspace.param_values` as needing refresh,
and, when applicable, rebuilds ensemble samplers so newly-specified
distributions take effect.  Calling [`ensure_functions!`](@ref) is handled
internally.

`group_name` may be a `Symbol` or `String`.  Supported payloads:

  * `Number` – literal value for scalar groups.
  * `Function` – definition for time/ensemble function groups (the latter are
    wrapped into a `QEnsembleFunction`).
  * `QDistribution` – distribution backing ensemble sampling.

Passing `nothing` is not supported; construct a fresh `QSpace` if you need to
remove a definition. Returns the modified `qspace` for convenience.
"""
function resolve_param!(qspace, group_name, payload)
    info = qspace.param_info
    groups = info.param_groups
    idx = get(qspace.parameter_dicts.group_name_to_index, Symbol(group_name)) do
        error("Unknown parameter group $(group_name).")
    end
    group = groups[idx]
    kind = group.kind
    old_payload = group.payload
    payload === nothing &&
        error("resolve_param! cannot remove payloads; construct a new QSpace if you need to clear definitions.")

    if kind == ParameterGroupScalar
        payload isa Number ||
            error("Scalar group $(group.name) expects a literal number payload.")
    elseif kind == ParameterGroupTimeFunction
        payload === nothing || payload isa Function ||
            error("Time-function group $(group.name) expects a Function payload.")
    elseif kind == ParameterGroupDistribution
        payload === nothing || payload isa QDistribution ||
            error("Distribution group $(group.name) expects a QDistribution payload.")
    elseif kind == ParameterGroupEnsembleFunction
        payload === nothing || payload isa QEnsembleFunction || payload isa Function ||
            error("Ensemble-function group $(group.name) expects a QEnsembleFunction or plain Function payload.")
    end

    assigned_payload = payload
    if kind == ParameterGroupEnsembleFunction
        if payload isa Function
            assigned_payload = QEnsembleFunction(String(group.name), group.function_args, payload)
        elseif payload isa QEnsembleFunction
            assigned_payload = QEnsembleFunction(String(group.name), group.function_args, payload.func)
        end
    end
    if kind in (ParameterGroupDistribution, ParameterGroupEnsembleFunction) &&
       old_payload !== nothing && assigned_payload !== nothing && old_payload !== assigned_payload
        @warn "Overwriting payload for parameter group $(group.name)."
    end

    group.payload = assigned_payload

    pv = qspace.param_values
    if kind == ParameterGroupDistribution
        pv.group_definition_initialized[idx] = assigned_payload isa QDistribution
        _reset_group_flags!(pv, idx)
    elseif kind == ParameterGroupScalar
        pv.group_definition_initialized[idx] = true
        pv.group_time_initialized[idx] .= true
        pv.group_initialized[idx] = true
    elseif kind == ParameterGroupTimeFunction
        pv.group_definition_initialized[idx] = assigned_payload isa Function
        _reset_group_flags!(pv, idx)
    elseif kind == ParameterGroupEnsembleFunction
        pv.group_definition_initialized[idx] = assigned_payload isa QEnsembleFunction
        _reset_group_flags!(pv, idx)
        if assigned_payload isa QEnsembleFunction
            time_group_idx = findfirst(==(Symbol("t")), info.outer_labels_symbols)
            for arg_pos in eachindex(assigned_payload.argument_symbols)
                arg_name = assigned_payload.argument_group_names[arg_pos]
                if arg_name == "t"
                    assigned_payload.argument_group_indices[arg_pos] = time_group_idx === nothing ? 0 : time_group_idx
                    continue
                end
                target_idx = findfirst(==(Symbol(arg_name)), info.outer_labels_symbols)
                target_idx === nothing &&
                    error("resolve_param!: ensemble function $(group.name) references unknown parameter group \"$arg_name\".")
                assigned_payload.argument_group_indices[arg_pos] = target_idx
                target_group = groups[target_idx]
                expected_len = length(target_group.indexes)
                self_len = length(assigned_payload.argument_self_index_positions[arg_pos])
                expected_len == self_len ||
                    error("resolve_param!: argument $(assigned_payload.argument_symbols[arg_pos]) expects $(self_len) index references but parameter group $(target_group.name) declares $(expected_len).")
            end
        end
    end

    pv.group_initialized[idx] = false
    pv.got_all_definitions = all(pv.group_definition_initialized)

    group_sym = Symbol(group.name)
    ensembles = qspace.ensembles
    ensemble_outer_lookup = IdDict{Any,Int}()
    for ss in qspace.subspaces
        ens = ss.ensemble
        ens === nothing && continue
        ensemble_outer_lookup[ens] = ss.ss_outer_ind
    end

    if kind == ParameterGroupDistribution
        for ens in ensembles
            idx ∈ ens.distribution_group_indices || continue
            dist_indices = ens.distribution_group_indices
            all_assigned = true
            dists = Vector{QDistribution}(undef, length(dist_indices))
            for (pos, gidx) in enumerate(dist_indices)
                payload = groups[gidx].payload
                if payload isa QDistribution
                    dists[pos] = payload
                else
                    all_assigned = false
                    break
                end
            end
            if all_assigned
                outer_symbols = info.outer_labels_symbols
                outer_names = info.outer_labels
                method = ens.sample_method === :default ?
                    (ens.as_continuum ? :chebychev : :random) : ens.sample_method
                sample = if ens.as_continuum
                    build_continuous_samples(ens, dist_indices, outer_symbols[dist_indices], outer_names[dist_indices], dists;
                        method=method)
                else
                    build_discrete_samples(ens, dist_indices, outer_symbols[dist_indices], outer_names[dist_indices], dists;
                        method=method,
                        num_nodes=ens.sample_num_nodes,
                        atol=ens.sample_atol,
                        rtol=ens.sample_rtol,
                        max_iter=ens.sample_max_iter)
                end
                ens.sampler = sample
                attach_samples!(pv, sample)
                outer_idx = get(ensemble_outer_lookup, ens, 0)
                if outer_idx != 0
                    sample_count = size(sample.samples, 1)
                    register_ensemble_sample_size!(pv, outer_idx, sample_count)
                end
            else
                ens.sampler = nothing
                for gidx in dist_indices
                    pv.group_initialized[gidx] = false
                    if !(groups[gidx].payload isa QDistribution)
                        pv.group_definition_initialized[gidx] = false
                    end
                end
            end
        end
    elseif kind == ParameterGroupEnsembleFunction
        for ens in ensembles
            idx ∈ ens.ensemble_function_group_indices || continue
            # nothing specific to update beyond payload assignment
        end
    elseif kind == ParameterGroupScalar && assigned_payload !== nothing
        _set_group!(pv, idx, assigned_payload; allow_function=true)
    end

    ensure_functions!(pv)
    return qspace
end
