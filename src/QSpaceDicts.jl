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

struct SubSpaceDicts
    by_outer::Dict{Symbol,Int}
    by_inner::Dict{Symbol,Tuple{Int,Int}}
end

struct AbstractOperatorDicts
    by_name::Dict{Symbol,Int}
end

function build_subspace_dicts(subspaces::Vector{SubSpace})::SubSpaceDicts
    outer_map = Dict{Symbol,Int}()
    inner_map = Dict{Symbol,Tuple{Int,Int}}()
    for (idx, ss) in enumerate(subspaces)
        if haskey(outer_map, ss.key_symbol)
            error("Duplicate outer subspace key $(ss.key_symbol) detected while building QSpace.")
        end
        outer_map[ss.key_symbol] = idx
        for (inner_idx, sym) in enumerate(ss.keys_symbols)
            if haskey(inner_map, sym)
                error("Duplicate inner subspace key $(sym) detected while building QSpace.")
            end
            inner_map[sym] = (idx, inner_idx)
        end
    end
    return SubSpaceDicts(outer_map, inner_map)
end

function build_operator_dicts(operatortypes::Vector{OperatorType})::AbstractOperatorDicts
    map = Dict{Symbol,Int}()
    for (idx, optype) in enumerate(operatortypes)
        if haskey(map, optype.name_sym)
            error("Duplicate operator type symbol $(optype.name_sym) detected while building QSpace.")
        end
        map[optype.name_sym] = idx
    end
    return AbstractOperatorDicts(map)
end
