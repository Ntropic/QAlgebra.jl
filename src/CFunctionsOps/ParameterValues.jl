const _GroupStorage = Union{ComplexF64, Array{ComplexF64}, Vector{Float64}}

export attach_samples!, register_ensemble_sample_size!, resolve_param!

using Base: WeakRef
using ..Sampler: QEnsembleFunction, QDistribution, build_discrete_samples, build_continuous_samples
using ..EnsembleSamples: AbstractEnsembleSample
import ..ParameterGroups
using ..ParameterGroups: ParameterGroup, ParameterGroupKind,
                          ParameterGroupScalar, ParameterGroupTimeScalar, ParameterGroupTimeFunction,
                          ParameterGroupDistribution, ParameterGroupEnsembleFunction,
                          ParameterGroupPayload

const _cached_qspaces = Ref{Union{Nothing,Module}}(nothing)

@inline function _qspaces_module()
    cached = _cached_qspaces[]
    if cached === nothing
        parent = parentmodule(@__MODULE__)
        isdefined(parent, :QSpaces) ||
            error("QSpaces module is not available; include QSpace.jl before calling resolve_param!")
        cached = getfield(parent, :QSpaces)
        _cached_qspaces[] = cached
    end
    return cached
end

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
    group_dependencies::Vector{Vector{Int}}
    group_update_waves::Vector{Vector{Int}}
    time_update_wave::Vector{Int}
end

function ParameterValues(param_info::ParameterInfo; qspace_ref::WeakRef=WeakRef())
    groups = param_info.param_groups
    where_which = ParameterGroups.WhereWhichParamGroup(groups)
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
        elseif group.kind == ParameterGroupTimeScalar
            group_values[g] = fill(Float64(NaN), time_count)
            group_definition_initialized[g] = true
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
    time_group = where_which.time_group
    if 0 < time_group <= group_count
        group_initialized[time_group] = true
    end
    update_order = _compute_update_order(groups, time_group)
    group_dependencies = _compute_group_dependencies(param_info)
    group_update_waves = _compute_group_update_waves(group_dependencies, update_order, group_count)
    time_update_wave = (0 < time_group <= group_count) ? group_update_waves[time_group] : Int[]

    pv = ParameterValues(qspace_ref, param_info, group_values, group_definition_initialized,
                         group_initialized, group_time_initialized, got_all_definitions,
                         update_order, time_group, group_dependencies, group_update_waves, time_update_wave)

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

function _compute_group_dependencies(info::ParameterInfo)
    groups = info.param_groups
    group_count = length(groups)
    deps = [Int[] for _ in 1:group_count]
    @inbounds for group_idx in 1:group_count
        group = groups[group_idx]
        for param_idx in group.parameter_indices
            refs = info.function_param_refs[param_idx]
            refs === nothing && continue
            for ref in refs
                dep_group = info.param_group_by_index[ref]
                dep_group == group_idx && continue
                push!(deps[group_idx], dep_group)
            end
        end
        if !isempty(deps[group_idx])
            sort!(deps[group_idx])
            unique!(deps[group_idx])
        end
    end
    return deps
end

function _compute_group_update_waves(group_dependencies::Vector{Vector{Int}},
                                     update_order::Vector{Int},
                                     group_count::Int)
    dependents = [Int[] for _ in 1:group_count]
    @inbounds for group_idx in 1:group_count
        for dep in group_dependencies[group_idx]
            push!(dependents[dep], group_idx)
        end
    end
    @inbounds for idx in 1:group_count
        if !isempty(dependents[idx])
            sort!(dependents[idx])
            unique!(dependents[idx])
        end
    end
    waves = Vector{Vector{Int}}(undef, group_count)
    queue = Vector{Int}()
    visited = BitVector()
    for start in 1:group_count
        empty!(queue)
        append!(queue, dependents[start])
        resize!(visited, group_count)
        fill!(visited, false)
        idx = 1
        while idx <= length(queue)
            current = queue[idx]
            idx += 1
            visited[current] && continue
            visited[current] = true
            append!(queue, dependents[current])
        end
        wave = Int[]
        if !isempty(update_order)
            for g in update_order
                if visited[g]
                    push!(wave, g)
                end
            end
        end
        waves[start] = wave
    end
    return waves
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

# Clear cached time-slot flags after the backing storage changes.
@inline function _reset_group_flags!(pv::ParameterValues, group_idx::Int)
    flags = pv.group_time_initialized[group_idx]
    fill!(flags, false)
    pv.group_initialized[group_idx] = false
end

# Mark a specific time slot as populated and update the group-level flag.
@inline function _mark_slot_initialized!(pv::ParameterValues, group_idx::Int, slot::Int)
    flags = pv.group_time_initialized[group_idx]
    slot = clamp(slot, 1, length(flags))
    flags[slot] = true
    pv.group_initialized[group_idx] = flags[1]
end

# Clamp time-slot requests to the storage bounds recorded for a group.
@inline function _clamp_slot(pv::ParameterValues, group_idx::Int, slot::Int)
    flags = pv.group_time_initialized[group_idx]
    return clamp(slot, 1, length(flags))
end

# Ensure that a group's storage contains an evaluated value for the requested slot.
@inline function _ensure_group_slot!(pv::ParameterValues, group_idx::Int, slot::Int)::Bool
    slot_idx = _clamp_slot(pv, group_idx, slot)
    flags = pv.group_time_initialized[group_idx]
    flags[slot_idx] && return true
    pv.group_definition_initialized[group_idx] || return false
    group = pv.param_info.param_groups[group_idx]
    if group.kind == ParameterGroupTimeFunction
        _evaluate_scalar_function_group!(pv, group_idx, slot_idx)
    elseif group.kind == ParameterGroupEnsembleFunction
        _evaluate_qensemble_group!(pv, group_idx, slot_idx)
    elseif group.kind == ParameterGroupDistribution
        return false
    else
        flags[slot_idx] = true
        pv.group_initialized[group_idx] = flags[1]
    end
    return flags[slot_idx]
end

# Locate the flat-vector position used by a distribution parameter within its group.
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

# Write a scalar parameter value into the appropriate storage buffer.
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
        if group.kind == ParameterGroupTimeScalar
            slot_idx = length(coords) >= 1 ? coords[1] : 1
            storage[slot_idx] = Float64(value)
        else
            position = _distribution_slot(info, group_idx, param_idx)
            storage[position] = Float64(value)
        end
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

function update_t!(pv::ParameterValues, value::Float64; slot::Int=0)
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
    _store_value!(pv, idx, value)
    slot_idx = slot + 1
    _refresh_time_dependents!(pv, slot_idx)
    return value
end
update_t!(qspace::QSpace, value::Float64; slot::Int=0) = update_t!(qspace.sample_index_param_values, value, slot)

set_time!(pv::ParameterValues, value::Float64; slot::Int=0) = update_t!(pv, value; slot=slot)

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
    if !_ensure_group_slot!(pv, group_idx, slot_idx)
        group = info.param_groups[group_idx]
        error("Parameter group $(group.display_signature) is not initialized; missing definition or samples.")
    end
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
    labels = Vector{String}(undef, group_count)
    sizes = Vector{String}(undef, group_count)
    for idx in 1:group_count
        group = groups[idx]
        pdf_hint = group.kind == ParameterGroupDistribution ? " (pdf)" : ""
        labels[idx] = group.display_signature * pdf_hint
        storage = pv.group_values[idx]
        sizes[idx] = storage isa ComplexF64 ? "1" :
                     storage isa Vector{Float64} ? string(length(storage)) :
                     storage isa Array{ComplexF64} ? (isempty(size(storage)) ? "1" : join(string.(size(storage)), "×")) :
                     string(typeof(storage))
    end
    name_width = isempty(labels) ? length("group") : max(length("group"), maximum(length, labels))
    size_width = isempty(sizes) ? length("size") : max(length("size"), maximum(length, sizes))
    def_hdr, init_hdr = "def", "init"
    order = collect(1:group_count)
    qspace = pv.qspace.value
    where_which = qspace === nothing ? ParameterGroups.WhereWhichParamGroup(groups) : qspace.where_which_param_groups
    time_group = where_which.time_group
    if 0 < time_group <= group_count
        order = vcat([time_group], filter(!=(time_group), order))
    end
    println(io, "ParameterValues:")
    println(io, "  ", rpad(def_hdr, 3), " ", rpad(init_hdr, 4), " ", rpad("group", name_width), "  ", rpad("size", size_width))
    for idx in order
        def_flag = pv.group_definition_initialized[idx] ? "✓" : "x"
        init_flag = pv.group_initialized[idx] ? "✓" : "x"
        println(io, "  ", rpad(def_flag, 3), " ", rpad(init_flag, 4), " ", rpad(labels[idx], name_width), "  ", rpad(sizes[idx], size_width))
    end
    if group_count > 0
        println(io, "  (def → payload ready, init → values populated)")
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
        args = Vector{Any}(undef, length(refs))
        for (pos, ref) in pairs(refs)
            dep_group = info.param_group_by_index[ref]
            dep_coords = info.param_coords[ref]
            _ensure_group_slot!(pv, dep_group, dep_coords[1]) || return
            args[pos] = _raw_value(pv, ref)
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
        args = Vector{Any}(undef, length(refs))
        for (pos, ref) in pairs(refs)
            dep_group = info.param_group_by_index[ref]
            dep_coords = info.param_coords[ref]
            _ensure_group_slot!(pv, dep_group, dep_coords[1]) || return
            args[pos] = _raw_value(pv, ref)
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
    _ensure_group_slot!(pv, time_group_idx, slot_idx) || error("Time parameter group is not initialized.")
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
    _ensure_group_slot!(pv, arg_group_idx, slot_idx) ||
        error("Parameter group $(group.display_signature) is not initialized; cannot evaluate ensemble function argument.")
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
        _ensure_group_slot!(pv, dep_group_idx, dep_slot) || return false
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

# Evaluate a time group (if present) for the requested slot.
@inline function _evaluate_time_group!(pv::ParameterValues, slot::Int)
    time_group = pv.time_group
    (0 < time_group <= length(pv.group_values)) || return pv
    group = pv.param_info.param_groups[time_group]
    slot_idx = _clamp_slot(pv, time_group, slot)
    if group.kind == ParameterGroupTimeFunction && group.payload isa Function
        _evaluate_scalar_function_group!(pv, time_group, slot_idx)
    elseif group.kind == ParameterGroupEnsembleFunction && group.payload isa QEnsembleFunction
        _evaluate_qensemble_group!(pv, time_group, slot_idx)
    end
    return pv
end

@inline function _evaluate_group_slot!(pv::ParameterValues, group_idx::Int, slot::Int)
    group = pv.param_info.param_groups[group_idx]
    slot_idx = _clamp_slot(pv, group_idx, slot)
    if group.kind == ParameterGroupTimeFunction && group.payload isa Function
        _evaluate_scalar_function_group!(pv, group_idx, slot_idx)
    elseif group.kind == ParameterGroupEnsembleFunction && group.payload isa QEnsembleFunction
        _evaluate_qensemble_group!(pv, group_idx, slot_idx)
    end
    return pv
end

function _run_update_wave!(pv::ParameterValues, wave::Vector{Int}, slot::Int)
    isempty(wave) && return pv
    for group_idx in wave
        pv.group_definition_initialized[group_idx] || continue
        deps = pv.group_dependencies[group_idx]
        all(pv.group_initialized[dep] for dep in deps) || continue
        _evaluate_group_slot!(pv, group_idx, slot)
    end
    pv.got_all_definitions = all(pv.group_definition_initialized)
    return pv
end

function _run_update_wave_all_slots!(pv::ParameterValues, wave::Vector{Int})
    isempty(wave) && return pv
    for group_idx in wave
        pv.group_definition_initialized[group_idx] || continue
        deps = pv.group_dependencies[group_idx]
        all(pv.group_initialized[dep] for dep in deps) || continue
        flags = pv.group_time_initialized[group_idx]
        if isempty(flags) || length(flags) == 1
            _evaluate_group_slot!(pv, group_idx, 1)
        else
            for slot_idx in eachindex(flags)
                _evaluate_group_slot!(pv, group_idx, slot_idx)
            end
        end
    end
    pv.got_all_definitions = all(pv.group_definition_initialized)
    return pv
end

@inline function _refresh_time_dependents!(pv::ParameterValues, slot::Int)
    _evaluate_time_group!(pv, slot)
    wave = pv.time_update_wave
    isempty(wave) || _run_update_wave!(pv, wave, slot)
    return pv
end

@inline function _refresh_group_dependents!(pv::ParameterValues, group_idx::Int; slot::Union{Nothing,Int}=nothing)
    wave = pv.group_update_waves[group_idx]
    isempty(wave) && return pv
    if slot === nothing
        _run_update_wave_all_slots!(pv, wave)
    else
        _run_update_wave!(pv, wave, slot)
    end
    return pv
end

function _refresh_group_and_dependents!(pv::ParameterValues, group_idx::Int; slot::Union{Nothing,Int}=nothing)
    pv.group_definition_initialized[group_idx] || return pv
    if slot === nothing
        flags = pv.group_time_initialized[group_idx]
        if isempty(flags) || length(flags) == 1
            _evaluate_group_slot!(pv, group_idx, 1)
        else
            for slot_idx in eachindex(flags)
                _evaluate_group_slot!(pv, group_idx, slot_idx)
            end
        end
    else
        _evaluate_group_slot!(pv, group_idx, slot)
    end
    _refresh_group_dependents!(pv, group_idx; slot=slot)
    return pv
end

"""
    recompute_functions!(pv::ParameterValues, slot::Int=1)

Refresh function-defined parameter groups so that every dependent sees the
latest inputs. The time group (if any) is evaluated first, followed by function
and ensemble groups in dependency order.
"""
function recompute_functions!(pv::ParameterValues, slot::Int=1)
    _refresh_time_dependents!(pv, slot)
    isempty(pv.update_order) || _run_update_wave!(pv, pv.update_order, slot)
    return pv
end

ensure_functions!(pv::ParameterValues) = recompute_functions!(pv)

# Adapt existing storage buffers to hold freshly attached ensemble samples.
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

# Resize array-backed ensemble storage once sample counts are known.
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

# Track how many samples belong to each index position of an ensemble group.
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
    if group.kind != ParameterGroupDistribution
        _reset_group_flags!(pv, group_idx)
    end

    if group.kind == ParameterGroupDistribution && all(>(0), group.sample_sizes)
        pv.group_definition_initialized[group_idx] = true
    end
end

"""
    register_ensemble_sample_size!(pv::ParameterValues, outer_idx::Int, sample_count::Int)

Propagate the sample count recorded for ensemble subspace `outer_idx` to every
parameter group that depends on it. Stored sample-size metadata is refreshed,
function-backed groups are marked for recomputation, and distribution groups
are marked as defined once all slots report a positive size. Returns `pv`.
"""
function register_ensemble_sample_size!(pv::ParameterValues, outer_idx::Int, sample_count::Int)
    info = pv.param_info
    for (group_idx, group) in enumerate(info.param_groups)
        positions = findall(==(outer_idx), group.index_outer_subspaces)
        isempty(positions) && continue
        _update_group_sample_sizes!(pv, group_idx, positions, sample_count)
        if group.kind == ParameterGroupDistribution && all(>(0), group.sample_sizes)
            pv.group_definition_initialized[group_idx] = true
        end
    end
    pv.got_all_definitions = all(pv.group_definition_initialized)
    return pv
end

"""
    attach_samples!(pv::ParameterValues, sample::AbstractEnsembleSample)

Attach the ensemble `sample` to its target parameter groups. Each group's
storage is resized or converted as needed, sample sizes are refreshed, and the
definition/initialisation flags are updated so downstream evaluations can
consume the data. Returns the input `sample`.
"""
function attach_samples!(pv::ParameterValues, sample::AbstractEnsembleSample)
    info = pv.param_info
    touched_groups = Int[]
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
        push!(touched_groups, group_idx)
    end
    if !isempty(touched_groups)
        sort!(touched_groups)
        unique!(touched_groups)
        for group_idx in touched_groups
            _refresh_group_dependents!(pv, group_idx)
        end
    end
    pv.got_all_definitions = all(pv.group_definition_initialized)
    return sample
end

"""
    resolve_param!(qspace, group_idx, payload)
    resolve_param!(qspace, name, payload)

Attach or update the payload for a parameter group on an existing `qspace`.
Both methods validate the payload against the group's declared kind, refresh
the corresponding storage in `qspace.sample_index_param_values`, and rebuild
ensemble samplers when distributions change. Calling [`ensure_functions!`](@ref)
is handled internally. Returns the modified `qspace` for convenience.
"""
function resolve_param!(qspace, group_idx::Integer, payload::ParameterGroupPayload)
    info = qspace.param_info
    groups = info.param_groups
    idx = Int(group_idx)
    1 <= idx <= length(groups) ||
        error("Parameter group index $(group_idx) out of bounds (1:$(length(groups))).")
    return _resolve_param_core!(qspace, idx, payload)
end

function resolve_param!(qspace, name::Union{Symbol,String}, payload::ParameterGroupPayload)
    qspaces = _qspaces_module()
    param_idx = qspaces.get_parameter_index(qspace, name)
    group_idx = qspace.param_info.param_group_by_index[param_idx]
    return resolve_param!(qspace, group_idx, payload)
end

function _resolve_param_core!(qspace, idx::Int, payload::ParameterGroupPayload)
    info = qspace.param_info
    groups = info.param_groups
    group = groups[idx]
    kind = group.kind
    old_payload = group.payload
    payload === nothing &&
        error("resolve_param! cannot remove payloads; construct a new QSpace if you need to clear definitions.")

    groups_to_refresh = Int[]

    if kind == ParameterGroupTimeScalar
        error("Time parameter group $(group.name) cannot be reassigned via resolve_param!; use set_time!/update_t! instead.")
    elseif kind == ParameterGroupScalar
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

    pv = qspace.sample_index_param_values
    if kind == ParameterGroupDistribution
        pv.group_definition_initialized[idx] = assigned_payload isa QDistribution
        _reset_group_flags!(pv, idx)
        push!(groups_to_refresh, idx)
    elseif kind == ParameterGroupScalar
        pv.group_definition_initialized[idx] = true
        pv.group_time_initialized[idx] .= true
        pv.group_initialized[idx] = true
    elseif kind == ParameterGroupTimeFunction
        pv.group_definition_initialized[idx] = assigned_payload isa Function
        _reset_group_flags!(pv, idx)
        push!(groups_to_refresh, idx)
    elseif kind == ParameterGroupEnsembleFunction
        pv.group_definition_initialized[idx] = assigned_payload isa QEnsembleFunction
        _reset_group_flags!(pv, idx)
        push!(groups_to_refresh, idx)
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

    if kind != ParameterGroupScalar
        pv.group_initialized[idx] = false
    end
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
                append!(groups_to_refresh, dist_indices)
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
        push!(groups_to_refresh, idx)
    end
    if !isempty(groups_to_refresh)
        sort!(groups_to_refresh)
        unique!(groups_to_refresh)
        for g in groups_to_refresh
            _refresh_group_and_dependents!(pv, g)
        end
    end
    return qspace
end
