
export get_parameter_group, get_subspace, get_subspace_index, get_ensemble, get_operator_type

# ==================> HELPERS <==========================================
@inline function _normalize_parameter_lookup(name::Union{Symbol,String})::Tuple{Symbol, String, String}
    normalized = var_unsubstitution(name)
    sym = isempty(normalized) ? Symbol(raw) : Symbol(normalized)
    return sym, String(name), normalized
end

@inline function _group_label(group::ParameterGroup)
    formatted, _ = symbol2formatted(String(group.name))
    plain = String(group.name)
    return "$(formatted) ($(plain))"
end
function _parameter_group_options(qspace::QSpace)
    [_group_label(group) for group in qspace.param_info.param_groups]
end

function _parameter_options(qspace::QSpace)
    unique(_parameter_group_options(qspace))
end
function _match_parameter_group_strings(qspace::QSpace, raw::String, normalized::String)
    matches = Int[]
    for (idx, group) in enumerate(qspace.param_info.param_groups)
        if raw == group.display_signature || raw == string(group.name)
            push!(matches, idx)
        elseif !isempty(normalized) && normalized == var_unsubstitution(group.display_signature)
            push!(matches, idx)
        end
    end
    return unique(matches)
end
function _match_parameter_strings(qspace::QSpace, raw::String, normalized::String)
    matches = Int[]
    for (idx, param) in enumerate(qspace.params)
        if raw == param.param_str || raw == param.param_name || raw == param.param_latex || raw == param.param_name_no_t
            push!(matches, idx)
            continue
        end
        if !isempty(normalized)
            norm_param = var_unsubstitution(param.param_str)
            if normalized == norm_param || normalized == var_unsubstitution(param.param_name) || normalized == var_unsubstitution(param.param_latex)
                push!(matches, idx)
            end
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
