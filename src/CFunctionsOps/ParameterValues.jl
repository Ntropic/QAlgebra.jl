const _GroupStorage = Union{ComplexF64, Array{ComplexF64}, Vector{Float64}}

using Base: WeakRef
using ..Sampler: QEnsembleFunction, QDistribution, build_discrete_samples, build_continuous_samples
using ..EnsembleSamples: AbstractEnsembleSample
using ..ParameterGroups: ParameterGroup, ParameterGroupKind,
                          ParameterGroupScalar, ParameterGroupTimeScalar, ParameterGroupTimeFunction,
                          ParameterGroupDistribution, ParameterGroupEnsembleFunction,
                          ParameterGroupPayload, WhereWhichParamGroup


# ==================================> Construction <==============================================================
"""
    ParameterValues(param_info::ParameterInfo; qspace_ref=WeakRef())

Concrete storage for the realised values of every parameter group across all
time slots and sample index combinations. The structure is intentionally thin –
it only holds the mutable data needed at run time while the immutable
descriptors remain with `ParameterGroup. 


For use cases that operate on ensemble abstract indices instead of static sample
positions, see `AbstractIndexParameters`, which wraps the same group metadata
with index-aware accessors.
"""
mutable struct ParameterValues
    where_which::WhereWhichParamGroup
    groups::Vector{ParameterGroup}
    group_values::Vector{_GroupStorage}
    group_definition_initialized::BitVector
    group_initialized::BitVector
    group_times_initialized::BitVector
    got_all_definitions::Bool
    time_group::Int
    group_dependencies::Vector{Vector{Int}}
    group_update_waves::Vector{Vector{Int}}
end
function ParameterValues(groups::Vector{ParameterGroup})
    where_which = WhereWhichParamGroup(groups)
    group_count = length(groups)

    group_values = Vector{_GroupStorage}(undef, group_count)
    group_definition_initialized = falses(group_count)
    group_initialized = falses(group_count)

    @inbounds for g in 1:group_count
        group = groups[g]
        time_count = max(group.time_count, 1)
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
    time_slot_count = (0 < time_group <= group_count) ? max(groups[time_group].time_count, 1) : 0
    group_times_initialized = time_slot_count == 0 ? BitVector() : falses(time_slot_count)
    update_order = _compute_update_order(groups, time_group)
    group_dependencies = _compute_group_dependencies(groups, time_group)
    group_update_waves = _compute_group_update_waves(group_dependencies, update_order, group_count, time_group)

    pv = ParameterValues(where_which,
                         groups,
                         group_values,
                         group_definition_initialized,
                         group_initialized,
                         group_times_initialized,
                         got_all_definitions,
                         time_group,
                         group_dependencies,
                         group_update_waves)

    #@inbounds for g in 1:group_count
    #    group = groups[g]
    #    if group.kind == ParameterGroupScalar && group.payload !== nothing
    #        _set_group!(pv, g, group.payload; allow_function=true)
    #    end
    #end
    #recompute_functions!(pv, 1)
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
    return vcat(function_groups, ensemble_groups)
end

function _compute_group_dependencies(groups::Vector{ParameterGroup}, time_group::Int)
    group_count = length(groups)
    deps = Vector{Vector{Int}}(undef, group_count)
    @inbounds for group_idx in 1:group_count
        curr_deps = unique(groups[group_idx].dependency_indices)
        if !any(==(time_group), curr_deps)
            for i in curr_deps
                if groups[i].of_t
                    push!(curr_deps, time_group)
                    break
                end
            end
        end
        sort!(curr_deps)
        deps[group_idx] = curr_deps
    end

    return deps
end


function _compute_group_update_waves(group_dependencies::Vector{Vector{Int}}, update_order::Vector{Int}, group_count::Int, time_group::Int)
    dependents = [Int[] for _ in 1:group_count]
    @inbounds for group_idx in 1:group_count
        for dep in group_dependencies[group_idx]
            push!(dependents[dep], group_idx)
        end
    end
    @inbounds for idx in 1:group_count
        if !isempty(dependents[idx])
            sort!(dependents[idx])
            unique!(sort!(dependents[idx]))
        end
    end
    waves = Vector{Vector{Int}}(undef, group_count)
    visited = BitVector()
    for start in 1:group_count
        queue = Vector{Int}()
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
        if start == time_group
            wave = copy(update_order)
        elseif !isempty(update_order)
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
    return zeros(Int, length(group.index_sizes)) # initialize with zeros if size is still unkown, at least we know the number of dimensions for the array
end


# ==================================> Display <==============================================================
function Base.show(io::IO, pv::ParameterValues)
    group_count = length(pv.groups)
    if get(io, :compact, false)
        print(io, "ParameterValues(", group_count, " groups)")
        return
    end
    labels = Vector{String}(undef, group_count)
    sizes = Vector{String}(undef, group_count)
    for idx in 1:group_count
        group = pv.groups[idx]
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
    if 0 < pv.time_group <= group_count
        order = vcat([time_group], filter(!=(time_group), order))
    end
    println(io, "ParameterValues:")
    println(io, "  ", rpad(def_hdr, 3), " ", rpad(init_hdr, 4), " ", rpad("group", name_width), "  ", rpad("size", size_width))
    for idx in order
        def_flag = pv.group_definition_initialized[idx] ? "✓" : "x"
        init_flag = pv.group_initialized[idx] ? "✓" : "x"
        println(io, "  ", rpad(def_flag, 3), " ", rpad(init_flag, 4), " ", rpad(labels[idx], name_width), "  ", rpad(sizes[idx], size_width))
    end
end

# ==================================> Helpers <==============================================================

# Mark a specific time slot as populated and update the group-level flag.
@inline function _mark_slot_initialized!(pv::ParameterValues, group_idx::Int, slot::Int)
    if group_idx == pv.time_group && !isempty(pv.group_times_initialized)
        slot_idx = clamp(slot, 1, length(pv.group_times_initialized))
        pv.group_times_initialized[slot_idx] = true
    end
    pv.group_initialized[group_idx] = true
end

# Clamp time-slot requests to the storage bounds recorded for a group.
@inline function _clamp_slot(pv::ParameterValues, group_idx::Int, slot::Int)
    group = pv.groups[group_idx]
    time_count = max(group.time_count, 1)
    return clamp(slot, 1, time_count)
end

# Ensure that a group's storage contains an evaluated value for the requested slot.
@inline function _ensure_group_slot!(pv::ParameterValues, group_idx::Int, slot::Int)::Bool
    slot_idx = _clamp_slot(pv, group_idx, slot)
    group = pv.param_groups[group_idx]
    if group_idx == pv.time_group && !isempty(pv.group_times_initialized)
        pv.group_times_initialized[slot_idx] && return true
    end
        pv.group_initialized[group_idx] && return true
    pv.group_definition_initialized[group_idx] || return false
    if group.kind == ParameterGroupTimeFunction
        _evaluate_scalar_function_group!(pv, group_idx, slot_idx)
    elseif group.kind == ParameterGroupEnsembleFunction
        _evaluate_qensemble_group!(pv, group_idx, slot_idx)
    elseif group.kind == ParameterGroupDistribution
        return false
    else
        if group_idx == pv.time_group && !isempty(pv.group_times_initialized)
            pv.group_times_initialized[slot_idx] = true
        end
        pv.group_initialized[group_idx] = true
    end
    return true
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

@inline function _compose_coords(slot_idx::Int, sample_indexes::Vector{Int})
    coords = Vector{Int}(undef, 1 + length(sample_indexes))
    coords[1] = slot_idx
    for (pos, val) in pairs(sample_indexes)
        coords[pos + 1] = val
    end
    return coords
end

@inline function _value_for_param(pv::ParameterValues, info::ParameterInfo, param_idx::Int)
    group_idx = info.param_group_by_index[param_idx]
    coords = info.param_coords[param_idx]
    group = info.param_groups[group_idx]
    time_index = group.of_t ? coords[1] - 1 : -1
    sample_vec = length(coords) > 1 ? Vector{Int}(coords[2:end]) : Int[]
    return value(pv, group_idx, time_index, sample_vec)
end

# Write a scalar parameter value into the appropriate storage buffer.
# ==================================> Group Assignment & Lookup <==============================================================
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
    if group_idx == pv.time_group && !isempty(pv.group_times_initialized)
        fill!(pv.group_times_initialized, true)
    end
    pv.group_initialized[group_idx] = true
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
    if group_idx == pv.time_group && !isempty(pv.group_times_initialized)
        fill!(pv.group_times_initialized, true)
    end
    pv.group_initialized[group_idx] = true
    return value
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
        storage_tuple = Tuple(coords)
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
set_time!(pv::ParameterValues, value::Float64; slot::Int=0) = update_t!(pv, value; slot=slot)

