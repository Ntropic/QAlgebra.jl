struct ParameterDicts
    group_name_to_index::Dict{Symbol,Int}
    param_name_to_indices::Dict{Symbol,Vector{Int}}
    time_slot_to_param::Dict{Int,Int}
end

@inline function _register_param_key!(dict::Dict{Symbol,Vector{Int}}, key::Symbol, idx::Int)
    key_str = String(key)
    isempty(key_str) && return nothing
    entry = get!(dict, key, Int[])
    in(idx, entry) || push!(entry, idx)
    return nothing
end

@inline function _normalize_param_placeholder(name::String)
    stripped = strip(name)
    isempty(stripped) && return nothing
    cleaned = replace(stripped, '(' => '_', ')' => "", '{' => '_', '}' => "", ',' => "_", ' ' => "")
    cleaned = strip(cleaned, '_')
    isempty(cleaned) && return nothing
    return Symbol(cleaned)
end

@inline function _numeric_param_key(base::String, coords::Vector{Int}, of_time::Bool)
    core = strip(base)
    isempty(core) && return nothing
    parts = String[core]
    time_idx = coords[1] - 1
    if of_time || time_idx != 0
        push!(parts, "t$(time_idx)")
    end
    for idx in coords[2:end]
        push!(parts, string(idx))
    end
    length(parts) == 1 && return nothing
    return Symbol(join(parts, "_"))
end

function build_parameter_dicts(info::ParameterInfo)::ParameterDicts
    group_name_to_index = Dict{Symbol,Int}()
    for (idx, sym) in enumerate(info.outer_labels_symbols)
        group_name_to_index[sym] = idx
    end

    param_name_to_indices = Dict{Symbol,Vector{Int}}()
    time_slot_to_param = Dict{Int,Int}()

    for (idx, param) in enumerate(info.params)
        coords = param.coords
        group_idx = param.group_index
        base_symbol = String(info.outer_labels_symbols[group_idx])
        sym_str = Symbol(param.param_str)
        sym_name = Symbol(param.param_name)
        _register_param_key!(param_name_to_indices, sym_str, idx)
        _register_param_key!(param_name_to_indices, sym_name, idx)
        placeholder = _normalize_param_placeholder(param.param_name)
        placeholder !== nothing && _register_param_key!(param_name_to_indices, placeholder, idx)
        numeric_key = _numeric_param_key(base_symbol, coords, param.param_of_t)
        numeric_key !== nothing && _register_param_key!(param_name_to_indices, numeric_key, idx)
        if param.is_t
            t_idx = coords[1] - 1
            _register_param_key!(param_name_to_indices, Symbol("t$(t_idx)"), idx)
            time_slot_to_param[t_idx] = idx
        end
    end

    return ParameterDicts(group_name_to_index, param_name_to_indices, time_slot_to_param)
end
