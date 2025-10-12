const _GroupStorage = Union{Float64, Array{Float64}, Vector{Float64}}

using Base: WeakRef
using ..Sampler: QEnsembleFunction, QDistribution, build_discrete_samples, build_continuous_samples
using ..EnsembleSamples: AbstractEnsembleSample
using ..ParameterGroups: ParameterGroup, ParameterGroupLike, ParameterGroupKind,
                          ParameterGroupScalar, ParameterGroupTimeScalar, ParameterGroupTimeFunction, ParameterGroupStorageUnion
                          ParameterGroupDistribution, ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction,
                          WhereWhichParamGroup, parameter_group_value_type, parameter_group_kind_name, parameter_group_storage_target_type
using Base: @propagate_inbounds, checkbounds

include("ParameterValuesOps/ParameterValues_setget.jl") # Fast paths for accessing values. 

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
    groups::Vector{ParameterGroupLike}
    group_values::Vector{GroupValsAny}
    group_definition_initialized::BitVector
    group_initialized::BitVector
    group_times_initialized::BitVector
    got_all_definitions::Bool
    time_group::Int
    how_many_times::Int
    group_dependencies::Vector{Vector{Int}}
    group_update_waves::Vector{Vector{Int}}
end
function ParameterValues(groups::AbstractVector{ParameterGroupLike})
    where_which = WhereWhichParamGroup(groups)
    group_count = length(groups)

    group_values = Vector{GroupValsAny}(undef, group_count)
    group_definition_initialized = falses(group_count)
    group_initialized = falses(group_count)

    @inbounds for (g, group) in enumerate(groups)
        val = construct_emtpy_arrays(group)
        group_values[g] = GroupVals(val, group.of_t) 
        if !isnothing(group.payload) || group.is_time_group
            group_definition_initialized[g] = true
        end
    end
    got_all_definitions = all(group_definition_initialized)

    time_group = where_which.time_group
    time_slot_count = (0 < time_group <= group_count) ? max(groups[time_group].time_count, 1) : 0
    group_times_initialized = time_slot_count == 0 ? BitVector() : falses(time_slot_count)
    how_many_times = groups[time_group].time_count

    update_order = _compute_update_order(groups, time_group)
    group_dependencies = _compute_group_dependencies(groups, time_group)
    group_update_waves = _compute_group_update_waves(group_dependencies, update_order, group_count, time_group)

    pv = ParameterValues(where_which, groups, group_values, group_definition_initialized, group_initialized,
                         group_times_initialized, got_all_definitions, time_group, how_many_times, group_dependencies, group_update_waves)
    return pv
end

function array_scaling(group::ParameterGroup)::Vector{Int}
    index_sizes = isempty(group.index_sizes) ? zeros(Int, length(group.sample_sizes)) : copy(group.index_sizes)
    return group.of_t ? vcat(group.time_count, index_sizes) : index_sizes
end 
function construct_emtpy_arrays(group::ParameterGroup{T}) where {T}
    array_dims = array_scaling(group)
    correct_type = parameter_group_value_type(group.kind)
    if correct_type <: Vector{Float64}
        @assert length(array_dims) == 1 "Vector type (as used by $(parameter_group_kind_name(group.kind)) - $(group.name)) requires exactly one dimension, got $(length(array_dims))"
        if group.kind == ParameterGroupTimeScalar
            return fill(NaN, array_dims[1])
        else
            return Vector{Float64}(undef,array_dims[1])
        end
    elseif correct_type == Array{Float64}
        @assert length(array_dims) > 0 "Parameter group $(group.name) of kind $(parameter_group_kind_name(group.kind)) requires arguments, either time or other parameters."
        return Array{Float64}(undef, array_dims...)
    elseif correct_type == Float64
        return NaN # Float64
    else
        error("Unhandled value type $(correct_type) for kind $(parameter_group_kind_name(group.kind)).")
    end
end

function _compute_update_order(groups::AbstractVector{ParameterGroupLike}, time_group::Int)
    function_groups = Int[]
    ensemble_groups = Int[]
    @inbounds for (idx, group) in enumerate(groups)
        idx == time_group && continue
        if group.kind == ParameterGroupTimeFunction
            push!(function_groups, idx)
        elseif group.kind == ParameterGroupEnsembleFunction || group.kind == ParameterGroupEnsembleTimeFunction
            push!(ensemble_groups, idx)
        end
    end
    return vcat(function_groups, ensemble_groups)
end

function _compute_group_dependencies(groups::AbstractVector{ParameterGroupLike}, time_group::Int)
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

# ==================================> Display <==============================================================
function Base.show(io::IO, pv::ParameterValues)
    group_count = length(pv.groups)
    time_group = pv.time_group
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
        sizes[idx] = storage isa Float64 ? "1" :
             (storage isa Vector{Float64} ? string(length(storage)) :
             (storage isa Array{Float64} ? join(string.(size(storage)), "×") : ""))

    end
    name_width = isempty(labels) ? length("group") : max(length("group"), maximum(length, labels))
    size_width = isempty(sizes) ? length("size") : max(length("size"), maximum(length, sizes))
    def_hdr, init_hdr = "def", "init"
    order = collect(1:group_count)
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
end

# ==================================> Get and Set <==============================================================

@propagate_inbounds function Base.getindex(pv::ParameterValues, g::Int, t::Int, I...)
    @boundscheck checkbounds(pv.group_values, g)
    V = pv.group_values[g]       # GroupVals{S,OT}
    @inbounds return V[t, I...]  # dispatches to the correct GroupVals method
end
@propagate_inbounds function Base.setindex!(pv::ParameterValues, val, g::Int, time_index::Int, I...)
    @boundscheck checkbounds(pv.group_values, g)
    V = pv.group_values[g]
    @inbounds V[time_index, I...] = val
    return val
end
@inbounds time(pv::ParameterValues, time_index::Int) = pv.group_values[pv.time_group].value[time_index]

# ==================================> Compute Values <===========================================================

# Conditions for computing values are met - needs a time_index and time_index free variant. 
# for a time dependent groups need to check which time_indexes are defined and give donitions for them. 
# use got_all_definitions (as the outer most check, since it gets rid of all checks except for time checks) 
# and the individual group_initialized, for a condition to be met, the values group_definition_initialized needs to be initialized. 
@inline function conditions_met_for_computing_values(pv::ParameterValues, group_index::Int)::Bool
    pv.got_all_definitions && return true
    if pv.group_definition_initialized[group_index]
        return all(@inbounds pv.group_initialized[i] for i in pv.group_dependencies[group_index])
    else   
        false
    end
end

# Compute the values for a group using its payload. 
@inline function payload2values!(pv::ParameterValues, g::Int, time_index::Int)::Nothing
    kind = pv.groups[g].kind :: ParameterGroupKind
    payload2values!(kind, pv, g, time_index) 
end
payload2values!(::Val{ParameterGroupTimeScalar}, pv::ParameterValues, g::Int, time_index::Int) = error("Values for time groups are set via update_t! ")
function payload2values!(::Val{ParameterGroupScalar}, pv::ParameterValues, g::Int, time_index::Int)::Nothing
    pv.group_values[g].value = pv.groups[g].payload
    return nothing
end 
function payload2values!(::Val{ParameterGroupTimeFunction}, pv::ParameterValues, g::Int, time_index::Int)::Nothing
    pv.group_values[g].value[time_index+1] = pv.groups[g].payload(time(pv, time_index))
    return nothing
end 
function payload2values!(::Val{ParameterGroupEnsembleFunction}, pv::ParameterValues, g::Int, time_index::Int)::Nothing
    payload = pv.groups[g].payload 
    curr_vals::Array{Float64} = pv.group_values[g].value
    curr_groups_vals::Vector{Vector{Float64}} = [pv.group_values[s].value for s in payload.argument_group_indices]
    # each argument is a Vector! we construct all index combination 
    for indexes in eachindex(curr_vals)
        curr_vals[indexes] = payload.func([v(i) for (v,i) in zip(curr_groups_vals, indexes)]...)
    end
    return nothing
end
function payload2values!(::Val{ParameterGroupEnsembleFunction}, pv::ParameterValues, g::Int, time_index::Int)::Nothing
    payload = pv.groups[g].payload
    output  = pv.group_values[g].value              # Array{Float64}
    arg_ids = payload.argument_group_indices        # indices of input groups
    arg_vecs = [pv.group_values[i].value for i in arg_ids]  # all argument vectors

    @inbounds for idx in CartesianIndices(output)
        # pick one element from each argument vector
        args = ntuple(j -> arg_vecs[j][idx.I[j]], length(arg_vecs))
        output[idx] = payload.func(args...)
    end
    return nothing
end
function payload2values!(::Val{ParameterGroupEnsembleTimeFunction}, pv::ParameterValues, g::Int, time_index::Int)::Nothing
    payload = pv.groups[g].payload
    output  = pv.group_values[g].value              # Array{Float64}
    arg_ids = payload.argument_group_indices
    arg_vecs = [pv.group_values[i].value for i in arg_ids]
    t = time_index + 1                              # external 0-based → 1-based storage
    sample_axes = ntuple(d -> axes(output, d + 1), ndims(output) - 1)
    @inbounds for idx in CartesianIndices(sample_axes)
        args = ntuple(j -> arg_vecs[j][idx.I[j]], length(arg_vecs))
        output[t, idx] = payload.func(time_index, args...)
    end
    return nothing
end

function conditional_t_update!(pv::ParameterValues, g::Int, time_index::Int=0)::Nothing
    if conditions_met_for_computing_values(pv, g)
        payload2values!(pv, g, time_index)
        pv.group_initialized[g] = true 
    end
    return nothing 
end

function conditional_payload2update(pv::ParameterValues, g::Int) 

end


function conditional_payload2update_of_time(pv::ParameterValues, g::Int, time_index::Int) 

end

# Add payload - store the new payload in groups - if conditions for update are met compute values - then start the update wave
# 2 variants: 
#   update_t!(param_values::ParameterValues, t_val::Real; time_index::Int=0) for updating a time values (time_indexes start at 0, but are stored at index 1 ) 
#   resolve_param!(param_values::ParameterValues, payload::Union{Number, Function, QDistribution})
#       - determines target type via parameter_group_storage_target_type => asserts for input 
#       - for each type we need to handle this separately => Function to QEnsembleFunction is already performed in QSpace_parameters.jl 
#       - after setting payload: trigger conditional update wave
function update_t!(pv::ParameterValues, value::Float64; slot::Int=0)::Nothing
    @assert slot >= 0 && slot < pv.how_many_times
    time_group = pv.time_group
    time_index = slot+1
    pv.group_values[time_group].value[time_index] = value 
    pv.group_times_initialized[time_index] = true 
    # update wave 
    for update_index in pv.group_update_waves[time_group] 
        conditional_t_update!(pv, update_index, slot)
    end
    return Nothing
end

function resolve_param!(pv::ParameterValues, g::Int, payload::Union{Number, Function, QDistribution, QEnsembleFunction})::Nothing
    group = pv.groups[g]
    kind = group.kind
    correct_type = parameter_group_value_type(kind)
    @assert typeof(payload) <: correct_type "The parameter group $(group.name) is a $(parameter_group_kind_name(kind)) and expects are payloads $(correct_type), got $(typeof(payload))."
    @assert kind != ParameterGroupTimeScalar "Don't set the time via resolve_param, use update_t!"
    @assert kind ∉ [ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction] "ParameterGroupEnsembleFunctions must be set during ParameterDefinitions."
    group.payload = payload 
    update_wave = pv.group_update_waves[g]
    possible_times = findall(pv.group_times_initialized)
    if group.of_t
        for time_index in possible_times
            conditional_payload2update_of_time!(pv, g, time_index)
        end
    else
        conditional_payload2update(pv, g)
    end
    for up in update_wave
        if pv.groups[up].of_time
            for time_index in possible_times
                payload2values!(pv, up, time_index)
            end
        else
            payload2values!(pv, up, 0)
        end
    end
    return Nothing
end