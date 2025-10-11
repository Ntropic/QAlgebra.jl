using Base: WeakRef
import .CFunctions
import .CFunctions: update_t!
import .ParameterGroups

const _GroupStorageAbstract = Union{ComplexF64, Array{ComplexF64}, Vector{Float64}}

struct AbstractIndexParameters
    qspace::WeakRef
    param_info::CFunctions.ParameterInfo
    group_values::Vector{_GroupStorageAbstract}
    where_which::WhereWhichParamGroup
end

@inline function _zero_storage(group::ParameterGroup)
    time_count = max(group.time_count, 1)
    idx_sizes = group.index_sizes
    if group.kind == ParameterGroups.ParameterGroupDistribution
        factor = group.of_t ? time_count : 1
        total = isempty(idx_sizes) ? factor : factor * prod(idx_sizes)
        return zeros(Float64, total)
    elseif group.kind == ParameterGroups.ParameterGroupTimeScalar
        return fill(Float64(NaN), time_count)
    else
        dims = if group.of_t
            isempty(idx_sizes) ? (time_count,) : (time_count, idx_sizes...)
        else
            isempty(idx_sizes) ? () : (idx_sizes...,)
        end
        if isempty(dims)
            return ComplexF64(0)
        else
            return zeros(ComplexF64, dims...)
        end
    end
end

@inline function _copy_group_value(storage)
    storage isa ComplexF64 && return storage
    storage isa Array{ComplexF64} && return copy(storage)
    storage isa Vector{Float64} && return copy(storage)
    return deepcopy(storage)
end

@inline _stored_group_count(aip::AbstractIndexParameters) = length(aip.group_values)
@inline _group_stored(aip::AbstractIndexParameters, idx::Int) = 1 <= idx <= length(aip.group_values)

@inline _coords_tuple(coords::AbstractVector{<:Integer}) = Tuple(coords)

@inline function _distribution_slot(info::CFunctions.ParameterInfo, group_idx::Int, param_idx::Int)
    params = info.param_groups[group_idx].parameter_indices
    pos = findfirst(==(param_idx), params)
    pos === nothing && error("Parameter $(info.params_str[param_idx]) does not belong to group $(info.param_groups[group_idx].name).")
    return pos
end

@inline function _store_param_value!(aip::AbstractIndexParameters, param_idx::Int, value)
    info = aip.param_info
    group_idx = info.param_group_by_index[param_idx]
    _group_stored(aip, group_idx) || error("Parameter group $(info.param_groups[group_idx].name) is not stored in AbstractIndexParameters.")
    group = info.param_groups[group_idx]
    storage = aip.group_values[group_idx]
    coords = info.param_coords[param_idx]

    if storage isa ComplexF64
        isempty(coords) || length(coords) == 1 ||
            error("Expected scalar storage for parameter $(info.params_name[param_idx]).")
        aip.group_values[group_idx] = ComplexF64(value)
    elseif storage isa Array{ComplexF64}
        idxs = if !isempty(group.sample_sizes) && ndims(storage) == length(coords) - 1
            _coords_tuple(@view coords[2:end])
        else
            _coords_tuple(coords)
        end
        storage[idxs...] = ComplexF64(value)
    elseif storage isa Vector{Float64}
        if group.kind == ParameterGroups.ParameterGroupTimeScalar
            slot_idx = isempty(coords) ? 1 : coords[1]
            slot_idx = clamp(slot_idx, 1, length(storage))
            storage[slot_idx] = Float64(value)
        else
            position = _distribution_slot(info, group_idx, param_idx)
            storage[position] = Float64(value)
        end
    else
        error("Unsupported storage type $(typeof(storage)) for group $(group_idx).")
    end
    return value
end

@inline function _raw_param_value(aip::AbstractIndexParameters, param_idx::Int)
    info = aip.param_info
    group_idx = info.param_group_by_index[param_idx]
    _group_stored(aip, group_idx) || error("Parameter group $(info.param_groups[group_idx].name) is not stored in AbstractIndexParameters.")
    storage = aip.group_values[group_idx]
    group = info.param_groups[group_idx]
    coords = info.param_coords[param_idx]

    if storage isa ComplexF64
        return storage
    elseif storage isa Array{ComplexF64}
        idxs = if !isempty(group.sample_sizes) && ndims(storage) == length(coords) - 1
            _coords_tuple(@view coords[2:end])
        else
            _coords_tuple(coords)
        end
        return storage[idxs...]
    elseif storage isa Vector{Float64}
        if group.kind == ParameterGroups.ParameterGroupTimeScalar
            slot_idx = isempty(coords) ? 1 : coords[1]
            slot_idx = clamp(slot_idx, 1, length(storage))
            return storage[slot_idx]
        else
            position = _distribution_slot(info, group_idx, param_idx)
            return storage[position]
        end
    end
    error("Unsupported storage type $(typeof(storage)) for group $(group_idx).")
end

"""
    AbstractIndexParameters(qspace::QSpace)

Create a lightweight snapshot of the sample-index parameter state associated
with `qspace`. Scalar groups are copied from the live `ParameterValues`, while
time/function/distribution groups before the first ensemble-function group are
initialised with empty storage so they can be updated independently.
"""
function AbstractIndexParameters(qspace::QSpace)
    pv = qspace.sample_index_param_values
    pv.got_all_definitions ||
        error("AbstractIndexParameters requires all parameter definitions to be resolved before construction.")
    info = pv.param_info
    where_which = qspace.where_which_param_groups
    group_count = length(info.param_groups)
    stored_count = clamp(where_which.first_qensemble_group - 1, 0, group_count)
    group_values = Vector{_GroupStorageAbstract}(undef, stored_count)
    stored_scalar_groups = Set(filter(idx -> idx <= stored_count, where_which.scalar_groups))
    for idx in 1:stored_count
        group_values[idx] = idx in stored_scalar_groups ?
            _copy_group_value(pv.group_values[idx]) :
            _zero_storage(info.param_groups[idx])
    end
    return AbstractIndexParameters(WeakRef(qspace), info, group_values, where_which)
end

function Base.show(io::IO, aip::AbstractIndexParameters)
    ws = aip.where_which
    first_idx = ws.first_qensemble_group
    groups = aip.param_info.param_groups
    if get(io, :compact, false)
        print(io, "AbstractIndexParameters(", length(groups), " groups, first ensemble function index = ", first_idx, ")")
        return
    end
    println(io, "AbstractIndexParameters (first ensemble function index = ", first_idx, "):")
    labels = [group.display_signature for group in groups]
    idx_labels = string.(eachindex(groups))
    sizes = Vector{String}(undef, length(groups))
    stored_count = _stored_group_count(aip)
    for (idx, group) in enumerate(groups)
        if idx <= stored_count
            value = aip.group_values[idx]
            size_str = if value isa ComplexF64
                "1"
            elseif value isa AbstractArray
                dims = size(value)
                isempty(dims) ? "1" : join(string.(dims), "×")
            else
                string(typeof(value))
            end
        else
            size_str = "(not stored)"
        end
        sizes[idx] = size_str
    end
    name_width = isempty(labels) ? length("group") : max(length("group"), maximum(length, labels))
    size_width = isempty(sizes) ? length("size") : max(length("size"), maximum(length, sizes))
    index_width = isempty(idx_labels) ? length("index") : max(length("index"), maximum(length, idx_labels))
    order = collect(eachindex(groups))
    time_group = ws.time_group
    if 0 < time_group <= length(groups)
        order = vcat([time_group], filter(!=(time_group), order))
    end
    println(io, "  ", rpad("index", index_width), " ", rpad("group", name_width), "  ", rpad("size", size_width))
    for idx in order
        println(io, "  ", rpad(idx_labels[idx], index_width), " ", rpad(labels[idx], name_width), "  ", rpad(sizes[idx], size_width))
    end
end

@inline function _refresh_value_groups!(aip::AbstractIndexParameters, pv::CFunctions.ParameterValues, groups::Vector{Int})
    for idx in groups
        _group_stored(aip, idx) || continue
        aip.group_values[idx] = _copy_group_value(pv.group_values[idx])
    end
    return aip
end

function _evaluate_time_function_groups!(aip::AbstractIndexParameters, slot_idx::Int)
    info = aip.param_info
    limit = _stored_group_count(aip)
    for group_idx in aip.where_which.time_function_groups
        group_idx <= limit || continue
        group = info.param_groups[group_idx]
        func = group.payload
        func isa Function || continue
        for param_idx in group.parameter_indices
            coords = info.param_coords[param_idx]
            if group.of_t
                isempty(coords) && continue
                coords[1] == slot_idx || continue
            end
            refs = info.function_param_refs[param_idx]
            result = if refs === nothing
                func()
            else
                args = Vector{Any}(undef, length(refs))
                for (pos, ref) in pairs(refs)
                    args[pos] = _raw_param_value(aip, ref)
                end
                func(args...)
            end
            _store_param_value!(aip, param_idx, result)
        end
    end
    return aip
end

function update_t!(aip::AbstractIndexParameters, t::Float64; slot::Int=0)
    time_group = aip.where_which.time_group
    stored_count = _stored_group_count(aip)
    (0 < time_group <= stored_count) ||
        error("Time parameter group is not stored within this AbstractIndexParameters.")
    storage = aip.group_values[time_group]
    storage isa Vector{Float64} ||
        error("Time parameter group does not use time-scalar storage in AbstractIndexParameters.")
    slot_idx = slot + 1
    1 <= slot_idx <= length(storage) ||
        error("Slot $(slot) is out of bounds for time parameter group (length $(length(storage))).")
    storage[slot_idx] = Float64(t)
    _evaluate_time_function_groups!(aip, slot_idx)
    return aip
end

"""
    set_distribution_params!(aip::AbstractIndexParameters, assignments; slot=0)
    set_distribution_params!(aip::AbstractIndexParameters, assignment1, assignment2...; slot=0)

Populate the stored distribution parameter groups using one or more assignment
tuples `(SubSpaceIndex, values)`. The specified `slot` selects the time slice
for time-dependent distributions. Assignment values are copied into both the
backing `ParameterValues` object and the `AbstractIndexParameters` snapshot.
Returns `aip`.
"""
function set_distribution_params!(aip::AbstractIndexParameters, assignments::AbstractVector; slot::Int=0)
    return _set_distribution_params!(aip, assignments; slot=slot)
end

function set_distribution_params!(aip::AbstractIndexParameters, assignments::Tuple{SubSpaceIndex,<:AbstractVector{<:Real}}...; slot::Int=0)
    return _set_distribution_params!(aip, assignments; slot=slot)
end

function _set_distribution_params!(aip::AbstractIndexParameters, assignments; slot::Int)
    qspace = aip.qspace.value
    qspace === nothing && error("AbstractIndexParameters does not reference a live QSpace; cannot assign distribution parameters.")
    pv = qspace.sample_index_param_values
    info = pv.param_info
    sub_info = info.subspace_info
    stored_count = _stored_group_count(aip)
    touched_groups = Int[]
    for entry in assignments
        (entry isa Tuple && length(entry) == 2) ||
            error("Distribution assignments must be tuples `(SubSpaceIndex, values)`.")
        sub_idx, raw_vals_any = entry
        sub_idx isa SubSpaceIndex ||
            error("Distribution assignments must reference SubSpaceIndex entries; got $(typeof(sub_idx)).")
        raw_vals = raw_vals_any isa AbstractVector ? raw_vals_any : collect(raw_vals_any)
        ensemble_idx = sub_info.ensemble_index_by_outer_index[sub_idx.outer]
        ensemble_idx != 0 || error("SubSpaceIndex $(sub_idx) does not reference an ensemble subspace.")
        vals = Float64.(raw_vals)
        for group_idx in aip.where_which.distribution_groups
            group_idx <= stored_count ||
                error("Distribution group index $(group_idx) is not stored in AbstractIndexParameters.")
            group = info.param_groups[group_idx]
            ensemble_idx in group.index_outer_subspaces || continue
            storage = pv.group_values[group_idx]
            storage isa Vector{Float64} ||
                error("Distribution group $(group.name) does not use flat vector storage.")
            positions = Int[]
            for (pos, param_idx) in enumerate(group.parameter_indices)
                tuples = info.param_index_tuples[param_idx]
                coords = info.param_coords[param_idx]
                if group.of_t
                    coords[1] == slot + 1 || continue
                end
                for (ens, inner) in tuples
                    if ens == ensemble_idx && inner == sub_idx.inner
                        push!(positions, pos)
                        break
                    end
                end
            end
            isempty(positions) && continue
            length(vals) == length(positions) ||
                error("Expected $(length(positions)) values for group $(group.name) at $(sub_idx); got $(length(vals)).")
            for (val_idx, storage_pos) in enumerate(positions)
                storage[storage_pos] = vals[val_idx]
            end
            if group.of_t
                local_slot = slot + 1
                CFunctions._mark_slot_initialized!(pv, group_idx, local_slot)
            else
                pv.group_initialized[group_idx] = true
            end
            pv.group_definition_initialized[group_idx] = true
            _group_stored(aip, group_idx) ||
                error("Distribution group $(group.name) is not stored in AbstractIndexParameters.")
            aip.group_values[group_idx] = _copy_group_value(storage)
            push!(touched_groups, group_idx)
        end
    end
    if !isempty(touched_groups)
        sort!(touched_groups)
        unique!(touched_groups)
        for group_idx in touched_groups
            group = info.param_groups[group_idx]
            if group.of_t
                slot_idx = CFunctions._clamp_slot(pv, group_idx, slot + 1)
                CFunctions._refresh_group_dependents!(pv, group_idx; slot=slot_idx)
            else
                CFunctions._refresh_group_dependents!(pv, group_idx)
            end
        end
    end
    pv.got_all_definitions = all(pv.group_definition_initialized)
    return aip
end
