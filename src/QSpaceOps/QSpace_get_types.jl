
export get_parameter_group, get_subspace, get_subspace_index, get_ensemble, get_operator_type
using ..ParameterGroups: ParameterGroupLike
using ..StringUtils: var_unsubstitution, symbol2formatted

# ==================> HELPERS <==========================================
@inline function _normalize_parameter_lookup(name::Union{Symbol,String})::Tuple{Symbol, String, String}
    raw = String(name)
    normalized = var_unsubstitution(raw)
    sym = isempty(normalized) ? Symbol(raw) : Symbol(normalized)
    return sym, raw, normalized
end

@inline function _group_label(group::ParameterGroupLike)
    formatted = group.param_str
    plain = group.param_raw
    return "$(formatted) ($(plain))"
end
function _parameter_group_options(qspace::QSpace)
    [_group_label(group) for group in qspace.param_info.param_groups]
end

function _match_parameter_group_strings(qspace::QSpace, raw::String, normalized::String)
    matches = Int[]
    for (idx, group) in enumerate(qspace.param_info.param_groups)
        signature_plain = default_group_signature(group, qspace.subspace_info)
        if raw == signature_plain || raw == group.param_raw || raw == string(group.param_symbol)
            push!(matches, idx)
        elseif !isempty(normalized) && normalized == var_unsubstitution(signature_plain)
            push!(matches, idx)
        end
    end
    return unique(matches)
end
function _resolve_subspace_location(qspace::QSpace, name::Symbol)
    dicts = qspace.subspace_dicts
    if haskey(dicts.by_outer, name)
        return qspace.subspaces[dicts.by_outer[name]], nothing
    elseif haskey(dicts.by_inner, name)
        idx, inner_idx = dicts.by_inner[name]
        return qspace.subspaces[idx], inner_idx
    else
        error("No subspace associated with key $(name).")
    end
end

# ==================> Getters <==========================================
function get_parameter_group(qspace::QSpace, name::Union{Symbol,String})
    sym, raw, normalized = _normalize_parameter_lookup(name)
    dict = qspace.parameter_dicts.group_name_to_index
    idx = get(dict, sym, nothing)
    if idx === nothing && sym != Symbol(raw)
        idx = get(dict, Symbol(raw), nothing)
    end
    if idx === nothing
        matches = _match_parameter_group_strings(qspace, raw, normalized)
        if isempty(matches)
            options = join(_parameter_group_options(qspace), ", ")
            norm_hint = (!isempty(normalized) && normalized != raw) ? " (normalized: \"$(normalized)\")" : ""
            error("No parameter group matching \"$(raw)\"$(norm_hint) registered in QSpace. Available groups: $(options).")
        elseif length(matches) > 1
            names = [_group_label(qspace.param_info.param_groups[m]) for m in matches]
            error("Parameter group name \"$(raw)\" is ambiguous. Matches: $(join(names, ", ")).")
        else
            idx = matches[1]
        end
    end
    return idx
end

function get_subspace(qspace::QSpace, name::Symbol)
    subspace, _ = _resolve_subspace_location(qspace, name)
    return subspace
end
get_subspace(qspace::QSpace, name::String) = get_subspace(qspace, Symbol(name))

function get_subspace_index(qspace::QSpace, name::Symbol)
    _, inner_idx = _resolve_subspace_location(qspace, name)
    inner_idx === nothing && error("Key $(name) refers to an outer subspace; no inner index to return.")
    return inner_idx
end
get_subspace_index(qspace::QSpace, name::String) = get_subspace_index(qspace, Symbol(name))

function get_ensemble(qspace::QSpace, name::Symbol)
    subspace = get_subspace(qspace, name)
    ens = subspace.ensemble
    ens === nothing && error("Subspace $(subspace.key) is not an ensemble.")
    return ens
end
get_ensemble(qspace::QSpace, name::String) = get_ensemble(qspace, Symbol(name))

function get_operator_type(qspace::QSpace, name::Symbol)
    idx = get(qspace.operator_dicts.by_name, name, nothing)
    idx === nothing && error("No operator type named $(name) registered in QSpace.")
    return qspace.operatortypes[idx]
end
get_operator_type(qspace::QSpace, name::String) = get_operator_type(qspace, Symbol(name))
