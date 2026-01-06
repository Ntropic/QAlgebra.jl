using Combinatorics
using SparseArrays
using ..CFunctions: ParameterInfo, ParameterValues
import ..CFunctions: param2string
using ..StringUtils: symbol2formatted, str2sub, t_suffix, normalize_underscore_indices, format_normalized_indices, split_index, var_unsubstitution, reverse_var_substitution
using ..SparsePermutationTools: SparsePermutation, denseperm
using ..Sampler: QDistribution, QEnsembleFunction, AbstractEnsembleSample
using ..ParameterGroups: ParameterGroupKind, ParameterGroup, ParameterGroupLike, AbstractSubSpace,
                          ParameterGroupScalar, ParameterGroupTimeScalar, ParameterGroupTimeFunction,
                          ParameterGroupDistribution, ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction,
                          parameter_group_input_type, parameter_group_storage_type, parameter_group_value_type, parameter_group_kind_name, 
                          WhereWhichParamGroup
using ..QIndexes: AbstractIndex, TimeIndex

@inline function _infer_group_kind(name::String, of_t::Bool, indices::Vector{String}, function_args::Vector{String})::ParameterGroupKind
    has_indices = !isempty(indices)
    has_args = !isempty(function_args)
    if name == "t"
        return ParameterGroupTimeScalar
    elseif has_indices
        if !has_args
            return ParameterGroupDistribution
        end
        return ("t" in function_args || of_t) ? ParameterGroupEnsembleTimeFunction : ParameterGroupEnsembleFunction
    elseif has_args || (of_t && name != "t")
        return ParameterGroupTimeFunction
    else
        return ParameterGroupScalar
    end
end

function _build_group_template(name::String, of_t::Bool, indices::Vector{String},
                               function_args::Vector{String}, kind::ParameterGroupKind)
    base_str, base_latex = symbol2formatted(name)
    param_symbol = Symbol(name)
    param_raw = String(name)
    return ParameterGroup(param_symbol, param_raw, base_str, base_latex, kind, of_t, copy(indices), copy(function_args),
                          String[], Int[], Int[], Int[], Int[], Tuple{Symbol,Symbol}[], Tuple{String,String}[], Int[], Int[], name == "t", nothing)
end

function _assign_group_definition!(group::ParameterGroup, payload)
    payload === nothing && return group
    kind = group.kind
    allowed = parameter_group_input_type(kind)
    payload isa allowed ||
        error("Parameter group $(group.param_symbol) expects a $(allowed) payload, got $(typeof(payload)).")
    if kind == ParameterGroupTimeScalar
        values = fill(Float64(NaN), 1)
        values[1] = Float64(payload)
        group.payload = values
    elseif kind == ParameterGroupEnsembleFunction || kind == ParameterGroupEnsembleTimeFunction
        if payload isa QEnsembleFunction
            expected = Symbol.(group.function_args)
            payload.argument_symbols == expected ||
                error("Ensemble function for $(group.param_symbol) expects arguments $(expected), got $(payload.argument_symbols).")
            group.payload = payload
        else
            group.payload = _build_qensemble_function(group.param_raw, group.indices, group.function_args, payload)
        end
    else
        group.payload = payload
    end
    storage_type = parameter_group_storage_type(kind)
    group.payload === nothing || group.payload isa storage_type ||
        error("Converted payload for $(group.param_symbol) does not match expected storage type $(storage_type) (got $(typeof(group.payload))).")
    return group
end

function _extract_group_payload(var::Union{AbstractString,Symbol})
    return String(var), nothing
end
function _extract_group_payload(var::Pair{T,V}) where {T<:Union{AbstractString,Symbol},V}
    return String(first(var)), last(var)
end
function _extract_group_payload(var::Tuple{T,V}) where {T<:Union{AbstractString,Symbol},V}
    return String(var[1]), var[2]
end
"""
    ParameterDefinitions(params...)

Construct a parameter definition list from symbols or strings. Use `(t)` to mark
time-dependent groups. Ensemble parameters must be provided together with either a
`QDistribution` or a `QEnsembleFunction`, supplied as `(definition, payload)` tuples
or `definition => payload`. Constructing a parameter that references ensemble indices
without one of these payloads throws an error. Distributions/functions are stored
once per parameter group and exposed via `qspace.param_info.param_groups`
after construction.

For example, providing a multi-ensemble function can be written as:

```
ParameterDefinitions( "alpha" => 2.0, 
                      "beta(t)" => t->t^2, 
                      "delta_i" => QUniform(0,1,10),
                      "eta_i",
                      "gamma_{i1,i2}(t, delta_i1, delta_i2)" => (t, gi, gj)->t*gi+gj)
```
"""
struct ParameterDefinitions
    groups::Vector{ParameterGroupLike}
    function ParameterDefinitions(params...)
        groups = ParameterGroupLike[]
        for var in params
            label, payload = _extract_group_payload(var)
            label_clean = reverse_var_substitution(label)
            pre, brace_elements = brace_separate(label_clean)
            name, indices = normalize_underscore_indices(pre)
            of_t = (name == "t") || ("t" in brace_elements)
            function_args = String.(brace_elements)
            index_tokens = String.(indices)
            inferred_kind = _infer_group_kind(name, of_t, index_tokens, function_args)
            group = _build_group_template(name, of_t, index_tokens, function_args, inferred_kind)
            payload === nothing || _assign_group_definition!(group, payload)
            push!(groups, group)
        end

        _finalize_group_dependencies!(groups)

        sort!(groups, by = g -> ((g.kind == ParameterGroupTimeScalar && g.param_symbol == :t) ? -1 : Int(g.kind), g.param_raw))

        return new(groups)
    end
end

function _finalize_group_dependencies!(groups::Vector{ParameterGroupLike})
    for group in groups
        deps = String[]
        if group.kind in (ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction, ParameterGroupTimeFunction)
            if group.kind == ParameterGroupTimeFunction && isempty(group.indices)
                invalid_args = String[]
                for arg in group.function_args
                    arg == "t" && continue
                    push!(invalid_args, arg)
                end
                if !isempty(invalid_args)
                    base_sig = group.param_raw
                    if !isempty(group.indices)
                        base_sig *= "_{" * join(group.indices, ",") * "}"
                    end
                    if group.of_t
                        base_sig *= "(t)"
                    end
                    message_missing = join(invalid_args, ", ")
                    error("Time function $(base_sig) without ensemble indices may only reference `t`; remove argument(s) $(message_missing) or convert the group to an ensemble function.")
                end
            end
            for arg in group.function_args
                if arg == "t"
                    push!(deps, "t")
                    continue
                end
                arg_name, _ = normalize_underscore_indices(arg)
                arg_group = _find_group_by_name(groups, arg_name)
                arg_group !== nothing ||
                    error("Function for parameter group $(group.param_symbol) references unknown group \"$arg_name\".")
                push!(deps, arg_name)
            end
        end
        if group.kind == ParameterGroupEnsembleFunction || group.kind == ParameterGroupEnsembleTimeFunction
            for arg in group.function_args
                arg == "t" && continue
                arg_name, _ = normalize_underscore_indices(arg)
                target = _find_group_by_name(groups, arg_name)
                target === nothing &&
                    error("Ensemble function for parameter group $(group.param_symbol) references unknown group \"$arg_name\".")
                if target.kind != ParameterGroupDistribution
                    error("Ensemble function for parameter group $(group.param_symbol) must reference distribution arguments; \"$arg_name\" is not distribution-backed.")
                end
            end
        end
        group.dependency_names = deps
    end
end

function _find_group_by_name(groups::Vector{ParameterGroupLike}, name::String)
    for group in groups
        if group.param_raw == name
            return group
        end
    end
    return nothing
end

@inline function _find_group_index(groups::Vector{ParameterGroupLike}, name::Union{String,Symbol})
    target = name isa Symbol ? String(name) : name
    @inbounds for (idx, group) in enumerate(groups)
        group.param_raw == target && return idx
    end
    return nothing
end

function _ensure_signature_alignment!(group::ParameterGroupLike, indices, function_args::Vector{String})
    parsed_indices = [String(x) for x in indices]
    parsed_args = [String(x) for x in function_args]
    parsed_indices == group.indices ||
        error("Signature for parameter group $(group.param_symbol) expects indices $(group.indices), got $(parsed_indices).")
    parsed_args == group.function_args ||
        error("Signature for parameter group $(group.param_symbol) expects arguments $(group.function_args), got $(parsed_args).")
    return nothing
end

function set_parameter_group_definition!(defs::ParameterDefinitions, signature::Union{AbstractString,Symbol}, payload)
    label = String(signature)
    label_clean = reverse_var_substitution(label)
    pre, brace_elements = brace_separate(label_clean)
    name, indices = normalize_underscore_indices(pre)
    group = _find_group_by_name(defs.groups, name)
    group === nothing && error("Unknown parameter group \"$name\" in ParameterDefinitions.")
    _ensure_signature_alignment!(group, indices, Vector{String}(brace_elements))
    _assign_group_definition!(group, payload)
    return group
end

set_parameter_group_definition!(defs::ParameterDefinitions, pair::Pair) = set_parameter_group_definition!(defs, pair[1], pair[2])

function AbstractIndex2string(subspace_info::SubSpaceInfo, index::AbstractIndex; do_latex::Bool=false, as_index::Bool=true)::String
    ensemble_labels = subspace_info.ensemble_labels[index.ensemble]
    index_symbol = ensemble_labels[1+index.summation] 
    subindex = index.slot == 0 ? "" : string(index.slot)
    if do_latex
        if as_index  
            return "_{"*index_symbol*"_{"*subindex*"}}"
        else
            return index_symbol*"_{"*subindex*"}"
        end
    else
        if as_index
            return str2sub(index_symbol*subindex)
        else
            return index_symbol * str2sub(subindex)
        end
    end
end



@inline function _format_group_argument(arg::String)::String
    arg == "t" && return "t"
    base, idxs = normalize_underscore_indices(arg)
    base_str = symbol2formatted(base)[1]
    return base_str * format_normalized_indices(idxs; do_latex=false)
end

function _format_group_signature(group::ParameterGroupLike)::String
    base_str = symbol2formatted(group.param_raw)[1]
    if !isempty(group.indices)
        base_str *= format_normalized_indices(group.indices; do_latex=false)
    end
    args = String[]
    if !isempty(group.function_args)
        for arg in group.function_args
            push!(args, _format_group_argument(arg))
        end
    elseif group.of_t && group.param_symbol != :t
        push!(args, "t")
    end
    if !isempty(args)
        base_str *= "(" * join(args, ",") * ")"
    end
    return base_str
end

function Base.show(io::IO, param_def::ParameterDefinitions)
    group_count = length(param_def.groups)
    var_str_vec = Vector{String}(undef, group_count)
    @inbounds for (idx, group) in enumerate(param_def.groups)
        var_str_vec[idx] = _format_group_signature(group)
    end
    println(io, "ParameterDefinitions: [" * join(var_str_vec, ", ") * "]")
end

function _build_qensemble_function(name::String, indices::Vector{String}, brace_elements::Vector{String}, f::Function)
    isempty(brace_elements) && error("Parameter group \"$name\" requires a parentheses list specifying argument order when providing an ensemble function, e.g. \"$name(t, alpha)\" => (t, alpha) -> ...")
    return QEnsembleFunction(name, indices, brace_elements, f)
end

function contiguous_blocks(v::Vector{Int})
    seen = Set{Int}()
    last = nothing
    for x in v
        if x != last
            if x in seen
                return false
            end
            push!(seen, x)
            last = x
        end
    end
    return true
end
function find_blocks(v::Vector{Int})
    blocks = Int[]
    lengths = Int[]
    i = 1
    while i ≤ length(v)
        val = v[i]
        len = 1
        i += 1
        while i ≤ length(v) && v[i] == val
            len += 1
            i += 1
        end
        push!(blocks, val)
        push!(lengths, len)
    end
    return blocks, lengths
end

function unformatted_var_name(param_name, indices::Vector{String})
    index_str = join(indices, ",")
    if length(indices) > 0 
        return param_name * "_{" * index_str * "}"
    else 
        return param_name * "_" * index_str 
    end
end

function build_subspace_index_maps(parameters::Vector{Any}, ensemble_index_maps, s_info::SubSpaceInfo)
    N = length(parameters)
    nsub = length(s_info.subspaces)
    result = Vector{Array{SparsePermutation,2}}(undef, nsub)

    ensemble_set = Set(s_info.where_ensembles)

    for s in 1:nsub
        if s in ensemble_set
            m = length(s_info.inner_labels_symbols[s])  # ensemble size for subspace s
            pairmat = Array{SparsePermutation}(undef, m, m)   # (inner_to, inner_from)
            for inner_to in 1:m, inner_from in 1:m
                idxs = Int[]
                diffs = Int[]
                @inbounds for (idx, param) in enumerate(parameters)
                    dest = idx
                    if param.indexed_param
                        has_s = false
                        inners = Vector{Int}(undef, length(param.param_indices))
                        for a in eachindex(param.param_indices)
                            ix = param.param_indices[a]
                            v = ix.inner
                            if ix.outer == s && ix.inner == inner_from
                                v = inner_to
                                has_s = true
                            end
                            inners[a] = v
                        end
                        if has_s
                            inners_sorted = sort(inners)
                            dest = ensemble_index_maps[param.group_index][param.t_index+1][inners_sorted...]
                        end
                    end
                    if dest != idx
                        push!(idxs, idx)
                        push!(diffs, dest - idx)
                    end
                end
                pairmat[inner_to, inner_from] = SparsePermutation(N, sparsevec(idxs, diffs, N))
            end
            result[s] = pairmat
        else
            # non-ensemble subspace -> empty (0×0) matrix
            result[s] = Array{SparsePermutation}(undef, 0, 0)
        end
    end
    for s in 1:nsub
        for mat in result[s]
            if any(==(0), denseperm(mat))
                error("build_subspace_index_maps: unmapped index (0) detected in subspace $s")
            end
        end
    end
    return result
end

function ensure_time_group!(groups::Vector{ParameterGroupLike}, used_symbols::Set{Symbol})
    if all(group.param_symbol != :t for group in groups) && !(:t in used_symbols)
        push!(groups, _build_group_template("t", true, String[], String[], ParameterGroupTimeScalar))
    end
end

function validate_group_payload(group::ParameterGroupLike)
    storage_type = parameter_group_storage_type(group.kind)
    payload = group.payload
    payload isa storage_type ||
        error("Parameter group $(group.param_symbol) payload does not match expected storage type $(storage_type) (got $(typeof(payload))).")
end

function resolve_ensemble_metadata!(group::ParameterGroupLike, subspace_info::SubSpaceInfo)
    indexes = group.indices
    if isempty(indexes)
        group.ensemble_indices = Int[]
        group.unique_ensemble_indices = Int[]
        group.subspace_indices = Int[]
        group.index_symbol_pairs = Tuple{Symbol,Symbol}[]
        group.index_string_pairs = Tuple{String,String}[]
        group.index_sizes = Int[]
        group.sample_sizes = Int[]
        return
    end

    outer_indices = Int[]
    for token in indexes
        outer = _find_ensemble_outer(subspace_info, token)
        outer === nothing && error("Unknown ensemble index \"$token\" referenced in parameter group $(group.param_symbol).")
        push!(outer_indices, outer)
    end

    group.subspace_indices = outer_indices
    group.ensemble_indices = subspace_info.ensemble_index_by_subspace_index[outer_indices]
    group.unique_ensemble_indices = unique(group.ensemble_indices)
    group.index_symbol_pairs = _build_index_symbol_pairs(outer_indices, subspace_info.subspaces)
    group.index_string_pairs = _build_index_string_pairs(outer_indices, subspace_info.subspaces)
    group.index_sizes = fill(1, length(indexes))
    group.sample_sizes = Int[]
end

function _find_ensemble_outer(subspace_info::SubSpaceInfo, token::String)
    base, _ = split_index(token)
    token_sym = Symbol(base)
    for (outer, sub) in enumerate(subspace_info.subspaces)
        if sub.is_ensemble_ss
            if sub.key == base || sub.key_symbol == token_sym
                return outer
            end
            if has_summation(sub)
                if secondary_label(sub) == base || secondary_symbol(sub) == token_sym
                    return outer
                end
            end
        end
    end
    return nothing
end

function _build_index_symbol_pairs(outer_indices::Vector{Int}, subspaces::Vector{SubSpace})
    count = length(outer_indices)
    pairs = Vector{Tuple{Symbol,Symbol}}(undef, count)
    for idx in 1:count
        sub = subspaces[outer_indices[idx]]
        pairs[idx] = (sub.key_symbol, secondary_symbol(sub))
    end
    return pairs
end

function _build_index_string_pairs(outer_indices::Vector{Int}, subspaces::Vector{SubSpace})
    count = length(outer_indices)
    pairs = Vector{Tuple{String,String}}(undef, count)
    for idx in 1:count
        sub = subspaces[outer_indices[idx]]
        pairs[idx] = (sub.key, secondary_label(sub))
    end
    return pairs
end

function register_group_with_subspaces!(group::ParameterGroupLike, subspaces::Vector{SubSpace})
    for outer_ind in unique(group.subspace_indices)
        outer_ind > 0 || continue
        ensemble_cfg = subspaces[outer_ind].ensemble
        if ensemble_cfg !== nothing && !(group.param_symbol in ensemble_cfg.param_groups)
            push!(ensemble_cfg.param_groups, group.param_symbol)
        end
    end
end

function collect_parameter_symbols(parameters::Vector{Any})
    return [param.param_symbol for param in parameters]
end

function collect_outer_label_strings(params_symbols::Vector{Symbol})
    params_raw = String.(params_symbols)
    params_str = Vector{String}(undef, length(params_symbols))
    params_latex = Vector{String}(undef, length(params_symbols))
    for (i, sym) in enumerate(params_symbols)
        base_str, base_latex = symbol2formatted(String(sym))
        params_str[i] = base_str
        params_latex[i] = base_latex
    end
    return params_raw, params_str, params_latex
end


function finalize_group_dependencies!(groups::Vector{ParameterGroupLike})
    for (g, group) in enumerate(groups)
        deps = Int[]
        for dep_name in group.dependency_names
            dep_idx = _find_group_index(groups, dep_name)
            dep_idx === nothing && error("Parameter group $(group.param_raw) depends on undefined group \"$dep_name\".")
            push!(deps, dep_idx::Int)
        end
        group.dependency_indices = deps
    end
    max_outer_idx = 0
    for group in groups
        if group.kind == ParameterGroupDistribution
            for outer_idx in group.unique_ensemble_indices
                outer_idx > 0 || continue
                max_outer_idx = max(max_outer_idx, outer_idx)
            end
        end
    end
    max_outer_idx == 0 && return
    dist_memberships = [Int[] for _ in 1:max_outer_idx]
    for (idx, group) in enumerate(groups)
        if group.kind == ParameterGroupDistribution
            for outer_idx in group.unique_ensemble_indices
                outer_idx > 0 || continue
                push!(dist_memberships[outer_idx], idx)
            end
        end
    end
    for members in dist_memberships
        length(members) <= 1 && continue
        for idx in members
            deps = groups[idx].dependency_indices
            for other in members
                other == idx && continue
                push!(deps, other)
            end
            unique!(deps)
            sort!(deps)
        end
    end
end

function resolve_function_arguments!(groups::Vector{ParameterGroupLike})
    time_group_idx = _find_group_index(groups, "t")
    for group in groups
        kind = group.kind
        payload = group.payload
        args = group.function_args
        indexes = group.indices
        if isempty(args) || !(payload isa QEnsembleFunction) ||
           !(kind in (ParameterGroupTimeFunction, ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction))
            continue
        end

        index_positions = Dict{String,Vector{Int}}()
        for (pos, name) in enumerate(indexes)
            push!(get!(index_positions, name, Int[]), pos)
        end
        required_tokens = Set(indexes)
        used_tokens = Set{String}()

        for (arg_pos, arg_str) in enumerate(args)
            if arg_str == "t"
                time_group_idx === nothing &&
                    error("Argument \"t\" in ensemble function $(group.param_symbol) requires a time parameter group.")
                payload.argument_group_indices[arg_pos] = time_group_idx
                continue
            end

            base, tokens = normalize_underscore_indices(arg_str)
            target_idx = _find_group_index(groups, base)
            target_idx === nothing &&
                error("Ensemble function $(group.param_symbol) references unknown group \"$base\".")
            target_group = groups[target_idx]
            if kind in (ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction) &&
               target_group.kind != ParameterGroupDistribution
                error("Ensemble function $(group.param_symbol) expects distribution arguments; \"$base\" is not distribution-backed.")
            end

            payload.argument_group_indices[arg_pos] = target_idx
            self_positions = payload.argument_self_index_positions[arg_pos]
            if isempty(tokens)
                isempty(self_positions) || error("Argument $arg_str should not specify indices for scalar group $base.")
            else
                length(tokens) == length(self_positions) ||
                    error("Argument $arg_str must specify $(length(self_positions)) index tokens for group $base.")
                for (tok, self_pos) in zip(tokens, self_positions)
                    positions = get(index_positions, tok) do
                        error("Index token $tok referenced in argument $arg_str is undefined for group $(group.param_symbol).")
                    end
                    self_pos in positions ||
                        error("Argument $arg_str references index position $self_pos which does not match token $tok in $(group.param_symbol).")
                    push!(used_tokens, tok)
                end
            end
        end

        if !isempty(required_tokens)
            missing = setdiff(required_tokens, used_tokens)
            if !isempty(missing)
                signature = group.param_raw
                if !isempty(indexes)
                    signature *= "_{" * join(indexes, ",") * "}"
                end
                if group.of_t
                    signature *= "(t)"
                end
                missing_list = join(collect(missing), ", ")
                error("Ensemble function for $(signature) does not reference index(es) $(missing_list).")
            end
        end
    end
    return groups
end

function normalize_ensemble_group_lists!(subspaces::Vector{SubSpace}, param_groups::Vector{ParameterGroupLike})
    group_count = length(param_groups)
    for subspace in subspaces
        ens = subspace.ensemble
        ens === nothing && continue
        dist_idxs = Int[]
        func_idxs = Int[]
        ordered_syms = Symbol[]
        seen_syms = Set{Symbol}()
        for group_sym in ens.param_groups
            if !(group_sym in seen_syms)
                push!(ordered_syms, group_sym)
                push!(seen_syms, group_sym)
            end
            group_idx = _find_group_index(param_groups, group_sym)
            group_idx === nothing && continue
            group_idx = group_idx::Int
            group = param_groups[group_idx]
            if group.kind == ParameterGroupDistribution
                group_idx in dist_idxs || push!(dist_idxs, group_idx)
            elseif group.kind == ParameterGroupEnsembleFunction || group.kind == ParameterGroupEnsembleTimeFunction
                group_idx in func_idxs || push!(func_idxs, group_idx)
            end
        end
        ens.param_groups = ordered_syms
        ens.distribution_group_indices = dist_idxs
        ens.ensemble_function_group_indices = func_idxs
    end
end

# Return parameter mapping vector for switching from t_index2 to t_index1.
map_by_tindex(t_index1::Int, t_index2::Int, pinfo::ParameterInfo) = denseperm(pinfo.t_index_transform[t_index1+1, t_index2+1])

function default_group_signature(group::ParameterGroupLike,
                                 subspace_info::SubSpaceInfo;
                                 do_latex::Bool=false)::String
    if isempty(group.subspace_indices)
        indices = AbstractIndex[]
    else
        indices = Vector{AbstractIndex}(undef, length(group.subspace_indices))
        @inbounds for (pos, outer) in enumerate(group.subspace_indices)
            ens = subspace_info.ensemble_index_by_subspace_index[outer]
            slot = 1
            indices[pos] = AbstractIndex(outer, ens, slot, false)
        end
    end
    time_idx = TimeIndex(0)
    return param2string(group, indices, time_idx, subspace_info; do_latex=do_latex)
end

function _synchronize_time_payload!(group::ParameterGroupLike)
    if group.kind != ParameterGroupTimeScalar
        return 0
    end
    desired_slots = 1
    existing = group.payload
    if existing === nothing
        group.payload = fill(Float64(NaN), desired_slots)
    elseif existing isa Vector{Float64}
        if length(existing) < desired_slots
            group.payload = vcat(existing, fill(Float64(NaN), desired_slots - length(existing)))
        end
    else
        converted = fill(Float64(NaN), desired_slots)
        converted[1] = Float64(existing)
        group.payload = converted
    end
    return length(group.payload)
end

function _initial_max_time_index(groups::Vector{ParameterGroupLike})::Int
    for group in groups
        group.is_time_group || continue
        payload = group.payload
        if payload isa Vector{Float64}
            return max(length(payload), 1) - 1
        elseif payload isa Array{Float64}
            return max(size(payload, 1), 1) - 1
        else
            return 0
        end
    end
    return -1
end

function _build_param_info(groups::Vector{ParameterGroupLike}, subspace_info::SubSpaceInfo)
    params_symbols = [group.param_symbol for group in groups]
    params_raw = [group.param_raw for group in groups]
    params_str = [group.param_str for group in groups]
    params_latex = [group.param_latex for group in groups]

    return ParameterInfo(
        params_symbols,
        params_raw,
        params_str,
        params_latex,
        subspace_info,
        groups,
    )
end

function ParameterDefinitions2Parameters(vd::ParameterDefinitions,
                                         subspace_info::SubSpaceInfo,
                                         subspaces::Vector{SubSpace},
                                         used_symbols::Set{Symbol})
    groups = deepcopy(vd.groups)
    ensure_time_group!(groups, used_symbols)
    sort!(groups, by = g -> (Int(g.kind), g.param_raw))

    for (idx, group) in enumerate(groups)
        _synchronize_time_payload!(group)
        resolve_ensemble_metadata!(group, subspace_info)
        register_group_with_subspaces!(group, subspaces)
    end

    finalize_group_dependencies!(groups)
    resolve_function_arguments!(groups)
    normalize_ensemble_group_lists!(subspaces, groups)

    param_info = _build_param_info(groups, subspace_info)

    initial_max_t = _initial_max_time_index(groups)
    sample_index_param_values = ParameterValues(param_info; max_t_ind=initial_max_t)
    parameter_dicts = build_parameter_dicts(param_info)

    return param_info, sample_index_param_values, parameter_dicts
end
