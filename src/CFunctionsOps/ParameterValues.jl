struct ParameterValues
    param_info::ParameterInfo
    group_values::Vector{Array{Any}}
    param_cache::Vector{Any}
    group_has_function::BitVector
    function_dirty::BitVector
    dependencies::Vector{Vector{Int}}
    source_groups::Vector{Vector{Int}}
    params_by_group::Vector{Vector{Int}}
    param_coords::Vector{Vector{Int}}
    function_param_refs::Vector{Union{Nothing,Vector{Int}}}
    group_name_map::Dict{Symbol,Int}
    param_name_map::Dict{Symbol,Int}
    time_param_lookup::Dict{Int,Int}
end

function ParameterValues(param_info::ParameterInfo)
    group_count = length(param_info.group_distributions)
    param_count = length(param_info.params_name)

    group_values = Vector{Array{Any}}(undef, group_count)
    for g in 1:group_count
        dims = (param_info.group_time_counts[g],)
        index_sizes = param_info.group_index_sizes[g]
        if !isempty(index_sizes)
            dims = (param_info.group_time_counts[g], index_sizes...)
        end
        arr = Array{Any}(undef, dims...)
        fill!(arr, nothing)
        group_values[g] = arr
    end

    param_cache = fill(nothing, param_count)

    group_has_function = BitVector(map(!isnothing, param_info.group_functions))
    function_dirty = copy(group_has_function)

    dependencies = [Int[] for _ in 1:group_count]
    source_groups = [Int[] for _ in 1:group_count]

    function_param_refs = param_info.function_param_refs
    params_by_group = param_info.params_by_group
    param_coords = param_info.param_coords

    for g in 1:group_count
        if !group_has_function[g]
            continue
        end
        for idx in params_by_group[g]
            refs = function_param_refs[idx]
            refs === nothing && continue
            for ref in refs
                src = param_info.param_group_by_index[ref]
                if src != g && !(g in dependencies[src])
                    push!(dependencies[src], g)
                end
                if !(src in source_groups[g])
                    push!(source_groups[g], src)
                end
            end
        end
    end

    group_name_map = Dict{Symbol,Int}(Symbol(param_info.outer_labels_symbols[g]) => g for g in 1:group_count)

    param_name_map = Dict{Symbol,Int}()
    for idx in 1:param_count
        sym = Symbol(param_info.params_str[idx])
        param_name_map[sym] = idx
        if param_info.param_is_t[idx]
            coord = param_coords[idx]
            t_idx = coord[1] - 1
            param_name_map[Symbol("t$(t_idx)")] = idx
        end
    end

    time_param_lookup = Dict{Int,Int}()
    for idx in 1:param_count
        if param_info.param_is_t[idx]
            t_idx = param_coords[idx][1] - 1
            time_param_lookup[t_idx] = idx
        end
    end

    ParameterValues(param_info, group_values, param_cache, group_has_function, function_dirty,
                    dependencies, source_groups, params_by_group, param_coords, function_param_refs,
                    group_name_map, param_name_map, time_param_lookup)
end

function _mark_dependents_dirty!(pv::ParameterValues, group_idx::Int)
    for dep in pv.dependencies[group_idx]
        pv.function_dirty[dep] = true
    end
end

function _coords_tuple(coords::Vector{Int})
    return Tuple(coords)
end

function _set_param_value!(pv::ParameterValues, param_idx::Int, value; mark_dependents::Bool=true)
    group_idx = pv.param_info.param_group_by_index[param_idx]
    pv.group_has_function[group_idx] && error("Cannot set value for function-defined parameter group $(pv.param_info.params_name[param_idx]).")
    coords = pv.param_coords[param_idx]
    pv.group_values[group_idx][_coords_tuple(coords)...] = value
    pv.param_cache[param_idx] = value
    if mark_dependents
        _mark_dependents_dirty!(pv, group_idx)
    end
end

function set_param!(pv::ParameterValues, param_idx::Int, value)
    _set_param_value!(pv, param_idx, value)
end

function _set_group!(pv::ParameterValues, group_idx::Int, value)
    pv.group_has_function[group_idx] && error("Cannot set values for function-defined parameter group $(pv.param_info.params_name[group_idx]).")
    arr = pv.group_values[group_idx]
    dims = size(arr)
    if isa(value, Number)
        fill!(arr, value)
    elseif isa(value, AbstractArray)
        size(value) == dims || error("Value shape $(size(value)) does not match expected $(dims) for group $(pv.param_info.params_str[ pv.params_by_group[group_idx][1] ]).")
        arr .= value
    else
        error("Unsupported value type $(typeof(value)) for group assignment.")
    end
    for idx in pv.params_by_group[group_idx]
        coords = pv.param_coords[idx]
        pv.param_cache[idx] = arr[_coords_tuple(coords)...]
    end
    pv.function_dirty[group_idx] = false
    _mark_dependents_dirty!(pv, group_idx)
end

function set_param!(pv::ParameterValues, name::Symbol, value)
    if haskey(pv.group_name_map, name)
        _set_group!(pv, pv.group_name_map[name], value)
        return
    end
    idx = get_parameter_index(pv, name)
    set_param!(pv, idx, value)
end

set_param!(pv::ParameterValues, name::String, value) = set_param!(pv, Symbol(name), value)

function set_time!(pv::ParameterValues, value)
    idx = get(pv.time_param_lookup, 0, nothing)
    idx === nothing && error("No time parameter t0 registered in ParameterValues.")
    _set_param_value!(pv, idx, value)
    for dep in pv.dependencies[pv.param_info.param_group_by_index[idx]]
        pv.function_dirty[dep] = true
    end
end

function get_parameter_index(pv::ParameterValues, name::Symbol)
    idx = get(pv.param_name_map, name, nothing)
    idx === nothing && error("No parameter named $(name) registered in ParameterValues.")
    return idx
end

get_parameter_index(pv::ParameterValues, name::String) = get_parameter_index(pv, Symbol(name))

function _get_param_value(pv::ParameterValues, param_idx::Int)
    val = pv.param_cache[param_idx]
    val === nothing && error("Parameter $(pv.param_info.params_str[param_idx]) is unset.")
    return val
end

function _compute_function_group!(pv::ParameterValues, group_idx::Int, stack::Vector{Int})
    pv.group_has_function[group_idx] || return
    pv.function_dirty[group_idx] || return
    if group_idx in stack
        error("Circular dependency detected among ensemble functions involving group $(Symbol(pv.param_info.outer_labels_symbols[group_idx])).")
    end
    push!(stack, group_idx)
    for src in pv.source_groups[group_idx]
        _compute_function_group!(pv, src, stack)
    end
    func = pv.param_info.group_functions[group_idx]
    refs_vector = pv.param_info.function_param_refs
    for idx in pv.params_by_group[group_idx]
        refs = refs_vector[idx]
        refs === nothing && continue
        args = map(refs) do ref
            _get_param_value(pv, ref)
        end
        value = func.func(args...)
        # update without marking dependents (they will be handled by caller)
        coords = pv.param_coords[idx]
        pv.group_values[group_idx][_coords_tuple(coords)...] = value
        pv.param_cache[idx] = value
    end
    pv.function_dirty[group_idx] = false
    pop!(stack)
end

function ensure_functions!(pv::ParameterValues)
    stack = Int[]
    for g in eachindex(pv.group_has_function)
        if pv.group_has_function[g] && pv.function_dirty[g]
            _compute_function_group!(pv, g, stack)
        end
    end
end

function param_value(pv::ParameterValues, param_idx::Int)
    ensure_functions!(pv)
    return _get_param_value(pv, param_idx)
end

function param_value(pv::ParameterValues, name::Symbol)
    idx = get_parameter_index(pv, name)
    return param_value(pv, idx)
end

param_value(pv::ParameterValues, name::String) = param_value(pv, Symbol(name))
