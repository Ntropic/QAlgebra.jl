const _GroupStorage = Union{ComplexF64, Array{ComplexF64}, Vector{Float64}}

export attach_samples!

using Base: WeakRef
using ..Sampler: QEnsembleFunction
using ..EnsembleSamples: AbstractEnsembleSample

const KIND_SCALAR = UInt8(1)
const KIND_TIME_FUNCTION = UInt8(2)
const KIND_DISTRIBUTION = UInt8(3)
const KIND_ENSEMBLE_FUNCTION = UInt8(4)

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
    got_all_definitions::Bool
    ensemble_group_samples::Vector{Union{Nothing,AbstractEnsembleSample}}
    ensemble_group_functions::Vector{Union{Nothing,QEnsembleFunction}}
    group_functions::Vector{Union{Nothing,Function}}
    update_order::Vector{Int}
    time_group::Int
end

function ParameterValues(param_info::ParameterInfo, ensemble_group_samples::Vector{Union{Nothing,AbstractEnsembleSample}}, ensemble_group_functions::Vector{Union{Nothing,QEnsembleFunction}}, group_functions::Vector{Union{Nothing,Function}}; qspace_ref::WeakRef=WeakRef())
    group_count = length(param_info.outer_labels_symbols)
    ensemble_group_functions = copy(ensemble_group_functions)
    length(ensemble_group_functions) == group_count ||
        error("Expected $(group_count) ensemble group functions, got $(length(ensemble_group_functions)).")
    group_functions = copy(group_functions)
    length(group_functions) == group_count ||
        error("Expected $(group_count) scalar group functions, got $(length(group_functions)).")
    length(ensemble_group_samples) == group_count ||
        error("Expected $(group_count) ensemble group samples, got $(length(ensemble_group_samples)).")
    ensemble_group_samples = copy(ensemble_group_samples)

    group_values = Vector{_GroupStorage}(undef, group_count)
    group_definition_initialized = falses(group_count)
    @inbounds for g in 1:group_count
        kind_code = param_info.group_kind_codes[g]
        time_count = param_info.group_time_counts[g]
        index_sizes = param_info.group_index_sizes[g]
        if kind_code == KIND_DISTRIBUTION
            count = length(param_info.params_by_group[g])
            group_values[g] = Vector{Float64}(undef, count)
            group_definition_initialized[g] = param_info.group_distributions[g] !== nothing
            continue
        end
        if isempty(index_sizes) && time_count == 1
            group_values[g] = ComplexF64(NaN)
            group_definition_initialized[g] = (group_functions[g] !== nothing) ||
                                              (param_info.group_initial_values[g] !== nothing)
            continue
        end
        dims = isempty(index_sizes) ? (time_count,) : (time_count, index_sizes...)
        group_values[g] = Array{ComplexF64}(undef, dims...)
        if kind_code == KIND_TIME_FUNCTION
            group_definition_initialized[g] = group_functions[g] !== nothing
        elseif kind_code == KIND_ENSEMBLE_FUNCTION
            group_definition_initialized[g] = ensemble_group_functions[g] !== nothing
        elseif param_info.group_initial_values[g] !== nothing
            group_definition_initialized[g] = true
        end
        if param_info.group_is_t[g]
            group_definition_initialized[g] = true
        end
    end

    group_initialized = falses(group_count)
    got_all_definitions = all(group_definition_initialized)
    time_group_idx = findfirst(param_info.group_is_t)
    time_group = time_group_idx === nothing ? 0 : time_group_idx
    update_order = _compute_update_order(param_info.group_kind_codes, time_group)

    pv = ParameterValues(qspace_ref, param_info, group_values, group_definition_initialized,
                         group_initialized, got_all_definitions, ensemble_group_samples,
                         ensemble_group_functions, group_functions, update_order, time_group)

    @inbounds for g in 1:group_count
        init_val = param_info.group_initial_values[g]
        init_val === nothing && continue
        _set_group!(pv, g, init_val; allow_function=true)
    end

    update_time_group(pv)
    update_functions(pv)
    update_ensemble_group_functions(pv)

    return pv
end

function _compute_update_order(group_kind_codes::Vector{UInt8}, time_group::Int)
    function_groups = Int[]
    ensemble_groups = Int[]
    @inbounds for (idx, code) in enumerate(group_kind_codes)
        if idx == time_group
            continue
        elseif code == KIND_TIME_FUNCTION
            push!(function_groups, idx)
        elseif code == KIND_ENSEMBLE_FUNCTION
            push!(ensemble_groups, idx)
        end
    end
    return isempty(function_groups) ?
        ensemble_groups :
        isempty(ensemble_groups) ? function_groups : vcat(function_groups, ensemble_groups)
end

@inline function _coords_tuple(coords::Vector{Int})
    return Tuple(coords)
end

@inline function _distribution_slot(info::ParameterInfo, group_idx::Int, param_idx::Int)
    params = info.params_by_group[group_idx]
    pos = findfirst(==(param_idx), params)
    pos === nothing && error("Parameter $(info.params_str[param_idx]) does not belong to group $(info.outer_labels_symbols[group_idx]).")
    return pos
end

@inline function _distribution_slot(info::ParameterInfo, group_idx::Int, coords::Vector{Int})
    params = info.params_by_group[group_idx]
    for (pos, idx) in enumerate(params)
        info.param_coords[idx] == coords && return pos
    end
    error("No parameter with coordinates $(coords) in group $(info.outer_labels_symbols[group_idx]).")
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
    pv.group_initialized[group_idx] = true
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
    pv.group_initialized[group_idx] = true
    return value
end

function _set_group_array!(pv::ParameterValues, group_idx::Int, value::AbstractArray)
    info = pv.param_info
    time_count = info.group_time_counts[group_idx]
    index_sizes = info.group_index_sizes[group_idx]
    expected_dims = isempty(index_sizes) ? (time_count,) : (time_count, index_sizes...)
    storage = pv.group_values[group_idx]
    if storage isa ComplexF64
        length(value) == 1 ||
            error("Value shape $(size(value)) does not match expected scalar for group $(info.outer_labels_symbols[group_idx]).")
        pv.group_values[group_idx] = ComplexF64(value[1])
    elseif storage isa Array{ComplexF64}
        isempty(index_sizes) && time_count == 1 && return _fill_group!(pv, group_idx, value[1])
        size(value) == expected_dims ||
            error("Value shape $(size(value)) does not match expected $(expected_dims) for group $(info.outer_labels_symbols[group_idx]).")
        storage .= ComplexF64.(value)
    elseif storage isa Vector{Float64}
        length(value) == length(storage) ||
            error("Value length $(length(value)) does not match expected $(length(storage)) for group $(info.outer_labels_symbols[group_idx]).")
        storage .= Float64.(value)
    else
        error("Group $(info.outer_labels_symbols[group_idx]) does not accept array assignments.")
    end
    pv.group_initialized[group_idx] = true
    return value
end

function set_param!(pv::ParameterValues, param_idx::Int, value)
    _store_value!(pv, param_idx, value)
    return value
end

function _set_group!(pv::ParameterValues, group_idx::Int, value; allow_function::Bool=false)
    info = pv.param_info
    if !allow_function &&
       (pv.group_functions[group_idx] !== nothing || pv.ensemble_group_functions[group_idx] !== nothing)
        error("Cannot assign values to function-defined parameter group $(info.outer_labels_symbols[group_idx]).")
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
    recompute_functions!(pv)
    return value
end

set_time!(pv::ParameterValues, value) = update_t!(pv, value, 0)

function _lookup_storage_value(info::ParameterInfo, storage, group_idx::Int, coords::Vector{Int}, param_label::String, initialized::Bool)
    if storage isa ComplexF64
        initialized || error("Parameter $(param_label) is unset.")
        return storage
    elseif storage isa Array{ComplexF64}
        idxs = _coords_tuple(coords)
        isassigned(storage, idxs...) ||
            error("Parameter $(param_label) is unset.")
        return storage[idxs...]
    elseif storage isa Vector{Float64}
        position = _distribution_slot(info, group_idx, coords)
        isassigned(storage, position) ||
            error("Parameter $(param_label) is unset.")
        return storage[position]
    else
        error("Unsupported storage type $(typeof(storage)) for $(param_label).")
    end
end

@inline function _raw_value(pv::ParameterValues, param_idx::Int)
    info = pv.param_info
    group_idx = info.param_group_by_index[param_idx]
    storage = pv.group_values[group_idx]
    if storage isa Vector{Float64}
        position = _distribution_slot(info, group_idx, param_idx)
        isassigned(storage, position) ||
            error("Parameter $(info.params_str[param_idx]) is unset.")
        return storage[position]
    end
    coords = info.param_coords[param_idx]
    return _lookup_storage_value(info, storage, group_idx, coords, info.params_str[param_idx], pv.group_initialized[group_idx])
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
    if storage isa Vector{Float64}
        position = _distribution_slot(info, group_idx, coords)
        isassigned(storage, position) ||
            error("Parameter $(info.params_str[param_idx]) is unset.")
        return storage[position]
    end
    return _lookup_storage_value(info, storage, group_idx, coords, info.params_str[param_idx], pv.group_initialized[group_idx])
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
    group_count = length(info.outer_labels_symbols)
    if get(io, :compact, false)
        print(io, "ParameterValues(", group_count, " groups)")
        return
    end
    println(io, "ParameterValues:")
    for idx in 1:group_count
        status = pv.group_initialized[idx] ? "✓" : "x"
        base_name = info.group_display_signatures[idx]
        kind_code = info.group_kind_codes[idx]
        pdf_hint = kind_code == KIND_DISTRIBUTION ? " (pdf)" : ""
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

@inline function _evaluate_scalar_function_group!(pv::ParameterValues, group_idx::Int)
    func = pv.group_functions[group_idx]
    func === nothing && return
    info = pv.param_info
    for idx in info.params_by_group[group_idx]
        refs = info.function_param_refs[idx]
        if refs === nothing
            result = func()
        else
            if any(!pv.group_initialized[info.param_group_by_index[ref]] for ref in refs)
                continue
            end
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
        if any(!pv.group_initialized[info.param_group_by_index[ref]] for ref in refs)
            continue
        end
        args = map(refs) do ref
            _raw_value(pv, ref)
        end
        result = func.func(args...)
        _store_value!(pv, idx, result; allow_function=true)
    end
    pv.group_initialized[group_idx] = true
end

function update_time_group(pv::ParameterValues)
    group_count = length(pv.group_values)
    time_group = pv.time_group
    (0 < time_group <= group_count) || return pv
    if pv.group_functions[time_group] !== nothing
        _evaluate_scalar_function_group!(pv, time_group)
    elseif pv.ensemble_group_functions[time_group] !== nothing
        _evaluate_qensemble_group!(pv, time_group)
    end
    return pv
end

function update_functions(pv::ParameterValues)
    isempty(pv.update_order) && return pv
    for g in pv.update_order
        g == pv.time_group && continue
        pv.group_functions[g] === nothing && continue
        _evaluate_scalar_function_group!(pv, g)
    end
    return pv
end

function update_ensemble_group_functions(pv::ParameterValues)
    isempty(pv.update_order) && return pv
    for g in pv.update_order
        pv.ensemble_group_functions[g] === nothing && continue
        _evaluate_qensemble_group!(pv, g)
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
function recompute_functions!(pv::ParameterValues)
    update_time_group(pv)
    update_functions(pv)
    update_ensemble_group_functions(pv)
    return pv
end

ensure_functions!(pv::ParameterValues) = recompute_functions!(pv)

param_value(pv::ParameterValues, args...) = value(pv, args...)


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

function attach_samples!(pv::ParameterValues, sample::AbstractEnsembleSample)
    for (col, group_idx) in enumerate(sample.group_indices)
        values = sample.samples[:, col]
        storage = pv.group_values[group_idx]
        pv.group_values[group_idx] = _assign_group_storage!(storage, values)
        pv.group_initialized[group_idx] = true
        pv.ensemble_group_samples[group_idx] = sample
        pv.group_definition_initialized[group_idx] = true
    end
    pv.got_all_definitions = all(pv.group_definition_initialized)
    return sample
end
