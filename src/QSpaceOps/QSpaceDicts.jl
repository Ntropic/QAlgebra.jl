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
    for (idx, sym) in enumerate(info.params_symbols)
        group_name_to_index[sym] = idx
    end
    return ParameterDicts(group_name_to_index, Dict{Symbol,Vector{Int}}(), Dict{Int,Int}())
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
        for (inner_idx, sym) in enumerate(subspace_symbols(ss))
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
