using Base: @propagate_inbounds
import ..SubSpaceIndex
using ..ParameterGroups: ParameterGroup, ParameterGroupLike, ParameterGroupScalar, ParameterGroupTimeScalar, ParameterGroupTimeFunction,
                          ParameterGroupDistribution, ParameterGroupEnsembleFunction,
                          ParameterGroupEnsembleTimeFunction

struct AbstractIndexMode <: ParameterValuesMode end

array_scaling(::Type{AbstractIndexMode}, group::ParameterGroup)::Vector{Int} =
    group.of_t ? vcat(group.time_count, copy(group.index_sizes)) : copy(group.index_sizes)

function construct_emtpy_arrays(::Type{AbstractIndexMode}, group::ParameterGroup{T}) where {T}
    dims = array_scaling(AbstractIndexMode, group)
    correct_type = parameter_group_value_type(group.kind)
    if correct_type <: Vector{Float64}
        @assert length(dims) == 1 "Vector type (as used by $(parameter_group_kind_name(group.kind)) - $(group.name)) requires exactly one dimension, got $(length(dims))"
        return fill(NaN, dims[1])
    elseif correct_type == Array{Float64}
        @assert !isempty(dims) "Parameter group $(group.name) of kind $(parameter_group_kind_name(group.kind)) requires arguments, either time or other parameters."
        return fill(NaN, Tuple(dims)...)  # Cartesian product of dims
    elseif correct_type == Float64
        return NaN
    else
        error("Unhandled value type $(correct_type) for kind $(parameter_group_kind_name(group.kind)).")
    end
end

function _compute_update_order(::Type{AbstractIndexMode}, groups::AbstractVector{ParameterGroupLike}, time_group::Int)
    function_groups = Int[]
    @inbounds for (idx, group) in enumerate(groups)
        idx == time_group && continue
        if group.kind == ParameterGroupTimeFunction
            push!(function_groups, idx)
        end
    end
    return function_groups
end

function _compute_group_update_waves(::Type{AbstractIndexMode}, group_dependencies::Vector{Vector{Int}}, update_order::Vector{Int}, group_count::Int, time_group::Int)
    return _compute_group_update_waves(SampleIndexMode, group_dependencies, update_order, group_count, time_group)
end

@propagate_inbounds function Base.getindex(pv::ParameterValues{AbstractIndexMode}, g::Int, t::Int, I::Vararg{Int})
    @boundscheck checkbounds(pv.group_values, g)
    kind = pv.groups[g].kind
    return abstract_get(Val(kind), pv, g, t, I...)
end

@inline time(pv::ParameterValues{AbstractIndexMode}, time_index::Int) = pv.group_values[pv.time_group].value[time_index+1]

function conditional_payload2update!(pv::ParameterValues{AbstractIndexMode}, g::Int)::Nothing
    group = pv.groups[g]
    kind = group.kind
    if kind == ParameterGroupDistribution || kind == ParameterGroupEnsembleFunction || kind == ParameterGroupEnsembleTimeFunction
        return nothing
    end
    if conditions_met_for_computing_values(pv, g)
        payload2values!(pv, g, 0)
        pv.group_initialized[g] = true
    end
    return nothing
end

function update_ensemble_sample_sizes!(::ParameterValues{AbstractIndexMode}, ::Int)::Nothing
    return nothing
end

function initialize_group_payload!(::Type{AbstractIndexMode}, pv::ParameterValues{AbstractIndexMode}, g::Int, payload)
    group = pv.groups[g]
    kind = group.kind
    if kind == ParameterGroupScalar
        pv.group_values[g].value = payload
        pv.group_initialized[g] = true
    elseif kind == ParameterGroupTimeFunction
        possible_times = findall(pv.group_times_initialized).-1
        if !isempty(possible_times)
            conditional_payload2update_of_time!(pv, g, possible_times)
        end
    end
    return nothing
end

@inline abstract_get(::Val{ParameterGroupScalar}, pv::ParameterValues{AbstractIndexMode}, g::Int, t::Int, I...) =
    pv.group_values[g][t, I...]

@inline abstract_get(::Val{ParameterGroupTimeScalar}, pv::ParameterValues{AbstractIndexMode}, g::Int, t::Int, I...) =
    pv.group_values[g][t, I...]

@inline abstract_get(::Val{ParameterGroupTimeFunction}, pv::ParameterValues{AbstractIndexMode}, g::Int, t::Int, I...) =
    pv.group_values[g][t, I...]

@inline abstract_get(::Val{ParameterGroupDistribution}, pv::ParameterValues{AbstractIndexMode}, g::Int, t::Int, I...) =
    pv.group_values[g][t, I...]

function abstract_get(::Val{ParameterGroupEnsembleFunction}, pv::ParameterValues{AbstractIndexMode}, g::Int, t::Int, I...)
    payload = pv.groups[g].payload
    arg_ids = payload.argument_group_indices
    args = ntuple(j -> pv[arg_ids[j], t, I[j]], length(arg_ids))
    return payload.func(args...)
end

function abstract_get(::Val{ParameterGroupEnsembleTimeFunction}, pv::ParameterValues{AbstractIndexMode}, g::Int, t::Int, I...)
    payload = pv.groups[g].payload
    arg_ids = payload.argument_group_indices
    time_val = time(pv, t)
    args = ntuple(j -> pv[arg_ids[j+1], t, I[j]], length(arg_ids) - 1)
    return payload.func(time_val, args...)
end

function resolve_ensemble_values!(pv::ParameterValues{AbstractIndexMode}, idx::SubSpaceIndex, values::Vector{Float64})
    outer = idx.outer
    maps = pv.ensemble_distribution_groups
    (1 <= outer <= length(maps)) || throw(ArgumentError("Invalid subspace index $(outer) for resolve_ensemble_values!"))
    groups = maps[outer]
    isempty(groups) && return

    inner = idx.inner
    @assert length(values) == length(groups) "Expected $(length(groups)) values for ensemble $(outer), got $(length(values))."
    for (local_idx, group_idx) in pairs(groups)
        storage = pv.group_values[group_idx].value::Vector{Float64}
        @inbounds storage[inner] = values[local_idx]
        pv.group_initialized[group_idx] = true
    end

    return nothing
end

function refresh_ensemble_values!(pv::ParameterValues{AbstractIndexMode}, idx::SubSpaceIndex)::Nothing
    outer = idx.outer
    maps = pv.ensemble_distribution_groups
    (1 <= outer <= length(maps)) || throw(ArgumentError("Invalid subspace index $(outer) for refresh_ensemble_values!"))
    groups = maps[outer]
    isempty(groups) && return nothing

    dependents = pv.group_update_waves
    for group_idx in groups
        for dep in dependents[group_idx]
            storage = pv.group_values[dep].value
            if storage isa Vector{Float64} || storage isa Array{Float64}
                fill!(storage, NaN)
            else
                pv.group_values[dep].value = NaN
            end
            pv.group_initialized[dep] = false
        end
    end

    return nothing
end
