const _GroupStorage = Union{Nothing, Number, AbstractArray}

export attach_samples!

using ..QDistributions: QEnsembleFunction
using ..EnsembleSamples: AbstractEnsembleSample

"""
    ParameterValues(param_info::ParameterInfo; ensemble_group_samples,
                    ensemble_group_functions, group_functions)

Concrete value store that mirrors the parameter layout described by
`param_info`. Each group owns either a scalar, a time vector, or a dense array
over `(time, indexes...)`.  The optional keyword arguments allow callers to pass
in precomputed ensemble sample descriptors (`ensemble_group_samples`), the
ensemble functions that produce indexed values (`ensemble_group_functions`), and
plain scalar functions for non-indexed groups (`group_functions`)—the same
values populated by `QSpace` during construction.  If omitted, entries default
to `nothing`.

`ParameterValues` keeps lightweight bookkeeping so time updates run in a
deterministic order: `time_group` points to the explicit time parameter,
`no_index_function_of_t` lists scalar time-driven groups, and
`indexed_function_of_t` collects the time-dependent ensemble groups.  All
ensemble-backed groups appear in `ensemble_group_functions`, and are recomputed
by [`recompute_functions!`](@ref) when dependencies change.
"""
struct ParameterValues
    param_info::ParameterInfo
    group_values::Vector{_GroupStorage}
    # Info about which group has which properties -> necessary to process time updates in the correct order. 
    group_initialized::BitVector
    ensemble_group_samples::Vector{Union{Nothing,AbstractEnsembleSample}}
    ensemble_group_functions::Vector{Union{Nothing,QEnsembleFunction}}
    group_functions::Vector{Union{Nothing,Function}}
    no_index_function_of_t::Vector{Int}
    indexed_function_of_t::Vector{Int}
    time_group::Int
end

@inline function _init_group_storage(time_count::Int, index_sizes::Vector{Int})
    if isempty(index_sizes)
        if time_count == 1
            return nothing
        else
            values = Vector{Any}(undef, time_count)
            fill!(values, nothing)
            return values
        end
    else
        dims = (time_count, index_sizes...)
        arr = Array{Any}(undef, dims...)
        fill!(arr, nothing)
        return arr
    end
end

@inline function _coords_tuple(coords::Vector{Int})
    return Tuple(coords)
end

@inline function _check_initialized(value, param_name::String)
    value === nothing && error("Parameter $(param_name) is unset.")
    return value
end

function ParameterValues(param_info::ParameterInfo;
                         ensemble_group_samples::Vector{Union{Nothing,AbstractEnsembleSample}}=fill!(Vector{Union{Nothing,AbstractEnsembleSample}}(undef, length(param_info.outer_labels_symbols)), nothing),
                         ensemble_group_functions::Vector{Union{Nothing,QEnsembleFunction}}=fill!(Vector{Union{Nothing,QEnsembleFunction}}(undef, length(param_info.outer_labels_symbols)), nothing),
                         group_functions::Vector{Union{Nothing,Function}}=fill!(Vector{Union{Nothing,Function}}(undef, length(param_info.outer_labels_symbols)), nothing))
    group_count = length(param_info.outer_labels_symbols)

    group_values = Vector{_GroupStorage}(undef, group_count)
    for g in 1:group_count
        time_count = param_info.group_time_counts[g]
        index_sizes = param_info.group_index_sizes[g]
        storage = _init_group_storage(time_count, index_sizes)
        group_values[g] = storage
    end

    ensemble_group_functions_map = BitVector(map(!isnothing, ensemble_group_functions))
    no_index_function_of_t = Int[]
    indexed_function_of_t = Int[]
    time_group = findfirst(param_info.group_is_t)
    for g in 1:group_count
        is_time_group = g == time_group
        has_indexes = !isempty(param_info.group_index_sizes[g])
        if param_info.group_of_t[g] && !is_time_group
            if !has_indexes && group_functions[g] !== nothing
                push!(no_index_function_of_t, g)
            elseif has_indexes && ensemble_group_functions_map[g]
                push!(indexed_function_of_t, g)
            end
        end
    end
    group_initialized = falses(group_count)

    time_group_index = time_group === nothing ? 0 : time_group

    return ParameterValues(param_info, group_values, group_initialized,
                           ensemble_group_samples, ensemble_group_functions, group_functions,
                           no_index_function_of_t, indexed_function_of_t,
                           time_group_index)
end

function _store_value!(pv::ParameterValues, param_idx::Int, value; allow_function::Bool=false)
    info = pv.param_info
    group_idx = info.param_group_by_index[param_idx]
    if !allow_function && (pv.group_functions[group_idx] !== nothing || pv.ensemble_group_functions[group_idx] !== nothing)
        name = info.outer_labels_symbols[group_idx]
        error("Cannot assign values to function-defined parameter group $(name).")
    end
    coords = info.param_coords[param_idx]
    storage = pv.group_values[group_idx]
    if storage === nothing || storage isa Number
        length(coords) == 1 || error("Expected scalar storage for parameter $(info.params_name[param_idx]).")
        pv.group_values[group_idx] = value
    else
        storage_tuple = _coords_tuple(coords)
        storage[storage_tuple...] = value
    end
    pv.group_initialized[group_idx] = true
    return value
end

function get_parameter_index(info::ParameterInfo, name::Symbol)
    matches = get(info.param_dicts.param_name_to_indices, name, nothing)
    matches === nothing && error("No parameter named $(name) registered in ParameterInfo.")
    length(matches) == 1 && return matches[1]
    labels = info.params_str[matches]
    error("Parameter name $(name) is ambiguous. Matches: $(join(labels, ", ")).")
end

get_parameter_index(info::ParameterInfo, name::String) = get_parameter_index(info, Symbol(name))
get_parameter_index(pv::ParameterValues, name) = get_parameter_index(pv.param_info, name)

function _fill_group!(pv::ParameterValues, group_idx::Int, value)
    storage = pv.group_values[group_idx]
    if storage === nothing || storage isa Number
        pv.group_values[group_idx] = value
    elseif storage isa AbstractArray
        storage .= value
    else
        error("Unsupported storage type $(typeof(storage)) for group $(group_idx).")
    end
    pv.group_initialized[group_idx] = true
    return value
end

function _set_group_array!(pv::ParameterValues, group_idx::Int, value::AbstractArray)
    info = pv.param_info
    time_count = info.group_time_counts[group_idx]
    index_sizes = info.group_index_sizes[group_idx]
    expected_dims = isempty(index_sizes) ? (time_count,) : (time_count, index_sizes...)
    if isempty(index_sizes) && time_count == 1
        length(value) == 1 ||
            error("Value shape $(size(value)) does not match expected scalar for group $(info.outer_labels_symbols[group_idx]).")
        pv.group_values[group_idx] = value[1]
    else
        size(value) == expected_dims ||
            error("Value shape $(size(value)) does not match expected $(expected_dims) for group $(info.outer_labels_symbols[group_idx]).")
        storage = pv.group_values[group_idx]
        storage isa AbstractArray || error("Group $(info.outer_labels_symbols[group_idx]) does not accept array assignments.")
        storage .= value
    end
    pv.group_initialized[group_idx] = true
    return value
end

function set_param!(pv::ParameterValues, param_idx::Int, value)
    _store_value!(pv, param_idx, value)
    return value
end

function _set_group!(pv::ParameterValues, group_idx::Int, value)
    info = pv.param_info
    (pv.group_functions[group_idx] === nothing && pv.ensemble_group_functions[group_idx] === nothing) ||
        error("Cannot assign values to function-defined parameter group $(info.outer_labels_symbols[group_idx]).")
    if value isa Number
        _fill_group!(pv, group_idx, value)
    elseif value isa AbstractArray
        _set_group_array!(pv, group_idx, value)
    else
        error("Unsupported value type $(typeof(value)) for group assignment.")
    end
    return value
end

function set_param!(pv::ParameterValues, name::Symbol, value)
    dicts = pv.param_info.param_dicts
    if haskey(dicts.group_name_to_index, name)
        return _set_group!(pv, dicts.group_name_to_index[name], value)
    end
    matches = get(dicts.param_name_to_indices, name, nothing)
    matches === nothing && error("No parameter named $(name) registered in ParameterValues.")
    length(matches) == 1 || error("Parameter name $(name) is ambiguous; matches indices $(join(string.(matches), ", ")).")
    return set_param!(pv, matches[1], value)
end

set_param!(pv::ParameterValues, name::String, value) = set_param!(pv, Symbol(name), value)

function update_t!(pv::ParameterValues, value, slot::Int=0)
    dicts = pv.param_info.param_dicts
    idx = get(dicts.time_slot_to_param, slot, nothing)
    idx === nothing && error("No time parameter t$(slot) registered in ParameterValues.")
    set_param!(pv, idx, value)
    recompute_functions!(pv)
    return value
end

set_time!(pv::ParameterValues, value) = update_t!(pv, value, 0)

function _lookup_storage_value(storage, coords::Vector{Int}, param_label::String)
    if storage === nothing || storage isa Number
        storage === nothing && error("Parameter $(param_label) is unset.")
        return storage
    else
        result = storage[_coords_tuple(coords)...]
        result === nothing && error("Parameter $(param_label) is unset.")
        return result
    end
end

@inline function _raw_value(pv::ParameterValues, param_idx::Int)
    info = pv.param_info
    group_idx = info.param_group_by_index[param_idx]
    storage = pv.group_values[group_idx]
    coords = info.param_coords[param_idx]
    return _lookup_storage_value(storage, coords, info.params_str[param_idx])
end

function value(pv::ParameterValues, param_idx::Int)
    recompute_functions!(pv)
    return _raw_value(pv, param_idx)
end

value(pv::ParameterValues, param_idx::Int, ::Nothing) = value(pv, param_idx)

function _resolve_index_coords(info::ParameterInfo, param_idx::Int, indexes::ConcreteIndexes)
    tuples = info.param_index_tuples[param_idx]
    coords = copy(info.param_coords[param_idx])
    for (offset, (ensemble, inner)) in enumerate(tuples)
        ensemble <= length(indexes.indexes) ||
            error("Concrete indexes missing ensemble $(ensemble) needed for $(info.params_str[param_idx]).")
        entries = indexes.indexes[ensemble]
        inner <= length(entries) ||
            error("Concrete indexes missing entry $(inner) for ensemble $(ensemble) in parameter $(info.params_str[param_idx]).")
        idx_val = entries[inner]
        idx_val > 0 ||
            error("Concrete index for ensemble $(ensemble) entry $(inner) is unset (<=0) for parameter $(info.params_str[param_idx]).")
        coords[offset + 1] = idx_val
    end
    return coords
end

function value(pv::ParameterValues, param_idx::Int, indexes::ConcreteIndexes)
    recompute_functions!(pv)
    info = pv.param_info
    coords = _resolve_index_coords(info, param_idx, indexes)
    group_idx = info.param_group_by_index[param_idx]
    storage = pv.group_values[group_idx]
    return _lookup_storage_value(storage, coords, info.params_str[param_idx])
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

@inline function _evaluate_scalar_function_group!(pv::ParameterValues, group_idx::Int)
    func = pv.group_functions[group_idx]
    func === nothing && return
    info = pv.param_info
    for idx in info.params_by_group[group_idx]
        refs = info.function_param_refs[idx]
        if refs === nothing
            result = func()
        else
            args = map(refs) do ref
                _raw_value(pv, ref)
            end
            result = func(args...)
        end
        _store_value!(pv, idx, result; allow_function=true)
    end
    pv.group_initialized[group_idx] = true
end

@inline function _evaluate_qensemble_group!(pv::ParameterValues, group_idx::Int)
    info = pv.param_info
    func = pv.ensemble_group_functions[group_idx]
    func === nothing && return
    for idx in info.params_by_group[group_idx]
        refs = info.function_param_refs[idx]
        refs === nothing && continue
        args = map(refs) do ref
            _raw_value(pv, ref)
        end
        result = func.func(args...)
        _store_value!(pv, idx, result; allow_function=true)
    end
    pv.group_initialized[group_idx] = true
end

"""
    recompute_functions!(pv::ParameterValues)

Evaluate all function-defined parameter groups attached to `pv`, respecting time
dependencies and ensemble ordering. Scalar groups driven by pure functions are
updated before time-dependent ensembles so that downstream evaluations see the
latest values.
"""
function recompute_functions!(pv::ParameterValues)
    group_count = length(pv.group_values)
    has_scalar = any(!isnothing, pv.group_functions)
    has_ensemble = any(!isnothing, pv.ensemble_group_functions)
    (has_scalar || has_ensemble) || return pv

    processed = falses(group_count)
    time_group = pv.time_group
    if 0 < time_group <= group_count
        processed[time_group] = true
    end

    # update time-dependent, non-ensemble groups first (if driven by a function)
    for g in pv.no_index_function_of_t
        pv.group_functions[g] === nothing && continue
        _evaluate_scalar_function_group!(pv, g)
        processed[g] = true
    end

    # then time-dependent ensemble groups
    for g in pv.indexed_function_of_t
        pv.ensemble_group_functions[g] === nothing && continue
        _evaluate_qensemble_group!(pv, g)
        processed[g] = true
    end

    # remaining scalar function groups
    for g in 1:group_count
        processed[g] && continue
        pv.group_functions[g] === nothing && continue
        _evaluate_scalar_function_group!(pv, g)
        processed[g] = true
    end

    # remaining ensemble function groups
    for g in 1:group_count
        processed[g] && continue
        pv.ensemble_group_functions[g] === nothing && continue
        if g == time_group
            continue
        end
        _evaluate_qensemble_group!(pv, g)
        processed[g] = true
    end

    return pv
end

ensure_functions!(pv::ParameterValues) = recompute_functions!(pv)

param_value(pv::ParameterValues, args...) = value(pv, args...)


function _assign_group_storage!(storage, values)
    if storage === nothing || storage isa Number
        return values
    elseif storage isa AbstractArray
        storage .= Ref(values)
        return storage
    else
        error("Unsupported storage type $(typeof(storage)) for ensemble sample assignment.")
    end
end

function attach_samples!(pv::ParameterValues, sample::AbstractEnsembleSample)
    for (col, group_idx) in enumerate(sample.group_indices)
        values = sample.samples[:, col]
        storage = pv.group_values[group_idx]
        pv.group_values[group_idx] = _assign_group_storage!(storage, values)
        pv.group_initialized[group_idx] = true
        pv.ensemble_group_samples[group_idx] = sample
    end
    return sample
end
