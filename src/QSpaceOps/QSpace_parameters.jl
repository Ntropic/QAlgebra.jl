using Combinatorics
using SparseArrays
using ..CFunctions: ParameterInfo, ParameterIndexes, ParameterDicts, ParameterValues, build_parameter_dicts
using ..StringUtils: symbol2formatted, str2sub, var_unsubstitution, var_unsubstitution, reverse_var_substitution
using ..SparsePermutationTools: SparsePermutation, denseperm
using ..Sampler: QDistribution, QEnsembleFunction, AbstractEnsembleSample
using ..ParameterGroups: ParameterGroupKind, ParameterGroup, ParameterGroupLike, AbstractSubSpace,
                          ParameterGroupScalar, ParameterGroupTimeScalar, ParameterGroupTimeFunction,
                          ParameterGroupDistribution, ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction,
                          parameter_group_input_type, parameter_group_storage_type, parameter_group_value_type, parameter_group_kind_name, 
                          WhereWhichParamGroup

""" 
    Parameter(param_name::String, param_of_t::Bool, var_of_ensemble::Bool, var_ensemble_index::Int=0;
              var_suffix::String="")

Container describing a single parameter instance in the `QSpace`. Ensemble distributions
are tracked per parameter group and materialise in `EnsembleSamples` stored on the ensembles.
"""
mutable struct Parameter
    param_symbol::Symbol
    param_name::String
    param_str::String
    param_name_no_t::String
    param_str_no_t::String
    param_latex::String
    index_comb_symbol::Vector{Symbol}
    param_of_t::Bool
    is_t::Bool
    t_index::Int 
    group_index::Int 
    indexed_param::Bool
    param_indexes::Vector{SubSpaceIndex}
    coords::Vector{Int}
    acts_on::Vector{BitVector}
    ensemble_tuples::Vector{Tuple{Int,Int}}
    function_refs::Union{Nothing,Vector{Int}}
    indexed_slot::Int
end

@inline function _format_argument_display(arg::String)
    arg == "t" && return "t"
    base, idxs = underscore_separate(arg)
    disp, _ = symbol2formatted(base, idxs)
    return disp
end
function _group_display_signature(name::String, indexes::Vector{String}, function_args::Vector{String})
    base_str, _ = symbol2formatted(name, indexes)
    isempty(function_args) && return base_str
    arg_strs = [_format_argument_display(arg) for arg in function_args]
    return base_str * "(" * join(arg_strs, ",") * ")"
end

@inline function _infer_group_kind(name::String, of_t::Bool, indexes::Vector{String}, function_args::Vector{String})::ParameterGroupKind
    has_indexes = !isempty(indexes)
    has_args = !isempty(function_args)
    if name == "t"
        return ParameterGroupTimeScalar
    elseif has_indexes
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

function _build_group_template(name::String, of_t::Bool, indexes::Vector{String},
                               function_args::Vector{String}, kind::ParameterGroupKind, display_signature::String)
    return ParameterGroup(Symbol(name), display_signature, kind, of_t, copy(indexes), copy(function_args), 
                            String[], Int[], Int[], Int[], BitVector(), Int[], AbstractSubSpace[],
                            of_t ? 0 : 1, Int[], Int[], name == "t", nothing)
end

function _assign_group_definition!(group::ParameterGroupLike, payload)
    payload === nothing && return group
    kind = group.kind
    allowed = parameter_group_input_type(kind)
    payload isa allowed ||
        error("Parameter group $(group.name) expects a $(allowed) payload, got $(typeof(payload)).")
    if kind == ParameterGroupTimeScalar
        slots = max(group.time_count, 1)
        values = fill(Float64(NaN), slots)
        values[1] = Float64(payload)
        group.payload = values
    elseif kind == ParameterGroupEnsembleFunction || kind == ParameterGroupEnsembleTimeFunction
        if payload isa QEnsembleFunction
            expected = Symbol.(group.function_args)
            payload.argument_symbols == expected ||
                error("Ensemble function for $(group.name) expects arguments $(expected), got $(payload.argument_symbols).")
            group.payload = payload
        else
            group.payload = _build_qensemble_function(String(group.name), group.indexes, group.function_args, payload)
        end
    else
        group.payload = payload
    end
    storage_type = parameter_group_storage_type(kind)
    group.payload === nothing || group.payload isa storage_type ||
        error("Converted payload for $(group.name) does not match expected storage type $(storage_type) (got $(typeof(group.payload))).")
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
or `definition => payload`. Constructing a parameter that references ensemble indexes
without one of these payloads throws an error. Distributions/functions are stored
once per parameter group and exposed via `qspace.param_info.param_groups`
after construction.

For example, providing a multi-ensemble function can be written as:

```
ParameterDefinitions(
    "gamma_{i,j}(t, alpha, beta)" => ((t, alpha, beta) -> t + alpha + beta),
)
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
            name, indexes = underscore_separate(pre)
            of_t = (name == "t") || ("t" in brace_elements)
            function_args = String.(brace_elements)
            index_tokens = String.(indexes)
            inferred_kind = _infer_group_kind(name, of_t, index_tokens, function_args)
            display_signature = _group_display_signature(name, index_tokens, function_args)
            group = _build_group_template(name, of_t, index_tokens, function_args, inferred_kind, display_signature)
            payload === nothing || _assign_group_definition!(group, payload)
            push!(groups, group)
        end

        _finalize_group_dependencies!(groups)

        sort!(groups, by = g -> ((g.kind == ParameterGroupTimeScalar && g.name == :t) ? -1 : Int(g.kind), String(g.name)))

        return new(groups)
    end
end

function Base.getproperty(defs::ParameterDefinitions, sym::Symbol)
    if sym === :var_param
        return getfield(defs, :groups)
    end
    return getfield(defs, sym)
end

function Base.propertynames(::ParameterDefinitions, private::Bool=false)
    names = (:groups,)
    return private ? (names..., :var_param) : names
end

function _finalize_group_dependencies!(groups::Vector{ParameterGroupLike})
    for group in groups
        deps = String[]
        if group.kind in (ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction, ParameterGroupTimeFunction)
            for arg in group.function_args
                if arg == "t"
                    push!(deps, "t")
                    continue
                end
                arg_name, _ = underscore_separate(arg)
                arg_group = _find_group_by_name(groups, arg_name)
                arg_group !== nothing ||
                    error("Function for parameter group $(group.name) references unknown group \"$arg_name\".")
                push!(deps, arg_name)
            end
        end
        if group.kind == ParameterGroupEnsembleFunction || group.kind == ParameterGroupEnsembleTimeFunction
            for arg in group.function_args
                arg == "t" && continue
                arg_name, _ = underscore_separate(arg)
                target = _find_group_by_name(groups, arg_name)
                target === nothing &&
                    error("Ensemble function for parameter group $(group.name) references unknown group \"$arg_name\".")
                if target.kind != ParameterGroupDistribution
                    error("Ensemble function for parameter group $(group.name) must reference distribution arguments; \"$arg_name\" is not distribution-backed.")
                end
            end
        end
        group.dependency_names = deps
    end
end

function _find_group_by_name(groups::Vector{ParameterGroupLike}, name::String)
    for group in groups
        if String(group.name) == name
            return group
        end
    end
    return nothing
end

@inline function _find_group_index(groups::Vector{ParameterGroupLike}, name::Union{String,Symbol})
    target = name isa Symbol ? String(name) : name
    @inbounds for (idx, group) in enumerate(groups)
        String(group.name) == target && return idx
    end
    return nothing
end

function _ensure_signature_alignment!(group::ParameterGroupLike,
                                      indexes,
                                      function_args::Vector{String})
    parsed_indexes = [String(x) for x in indexes]
    parsed_args = [String(x) for x in function_args]
    parsed_indexes == group.indexes ||
        error("Signature for parameter group $(group.name) expects indexes $(group.indexes), got $(parsed_indexes).")
    parsed_args == group.function_args ||
        error("Signature for parameter group $(group.name) expects arguments $(group.function_args), got $(parsed_args).")
    return nothing
end

function set_parameter_group_definition!(defs::ParameterDefinitions, signature::Union{AbstractString,Symbol}, payload)
    label = String(signature)
    label_clean = reverse_var_substitution(label)
    pre, brace_elements = brace_separate(label_clean)
    name, indexes = underscore_separate(pre)
    group = _find_group_by_name(defs.groups, name)
    group === nothing && error("Unknown parameter group \"$name\" in ParameterDefinitions.")
    _ensure_signature_alignment!(group, indexes, Vector{String}(brace_elements))
    _assign_group_definition!(group, payload)
    return group
end

set_parameter_group_definition!(defs::ParameterDefinitions, pair::Pair) =
    set_parameter_group_definition!(defs, pair[1], pair[2])


function Base.show(io::IO, param_def::ParameterDefinitions)
    var_str_vec = []
    for group in param_def.groups
        param_str = symbol2formatted(String(group.name))[1]
        if group.of_t
            param_str *= "(t)"
        end
        param_str *= str2sub(join(group.indexes, ","))
        push!(var_str_vec, param_str) 
    end 
    println(io, "ParameterDefinitions: [" * join(var_str_vec, ", ") * "]")
end

function _build_qensemble_function(name::String, indexes::Vector{String}, brace_elements::Vector{String}, f::Function)
    isempty(brace_elements) && error("Parameter \"$name\" requires a parentheses list specifying argument order when providing an ensemble function, e.g. \"$name(t, alpha)\" => (t, alpha) -> ...")
    return QEnsembleFunction(name, indexes, brace_elements, f)
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

function unformatted_var_name(param_name, indexes::Vector{String})
    index_str = join(indexes, ",")
    if length(indexes) > 0 
        return param_name * "_{" * index_str * "}"
    else 
        return param_name * "_" * index_str 
    end
end

function build_subspace_index_maps(parameters::Vector{Parameter}, ensemble_index_maps, s_info::SubSpaceInfo)
    N = length(parameters)
    nsub = length(s_info.outer_labels_symbols)
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
                        inners = Vector{Int}(undef, length(param.param_indexes))
                        for a in eachindex(param.param_indexes)
                            ix = param.param_indexes[a]
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
    if all(group.name != :t for group in groups) && !(:t in used_symbols)
        display_signature = _group_display_signature("t", String[], String[])
        push!(groups, _build_group_template("t", true, String[], String[], ParameterGroupTimeScalar, display_signature))
    end
end

function validate_group_payload(group::ParameterGroupLike)
    storage_type = parameter_group_storage_type(group.kind)
    payload = group.payload
    payload isa storage_type ||
        error("Parameter group $(group.name) payload does not match expected storage type $(storage_type) (got $(typeof(payload))).")
end

function resolve_ensemble_metadata!(group::ParameterGroupLike, subspace_info::SubSpaceInfo)
    index_syms = Symbol.(group.indexes)
    if isempty(index_syms)
        group.ensemble_outer_indices = Int[]
        group.index_outer_subspaces = Int[]
        group.sample_sizes = Int[]
        presence = BitVector(undef, length(subspace_info.where_ensembles))
        fill!(presence, false)
        group.ensemble_presence = presence
        return Int[], Int[]
    end

    outer_subsystem_inds = -ones(Int, length(index_syms))
    for (i, index_sym) in enumerate(index_syms)
        subsystem_ind = nothing
        for (outer_ind, outer_labels_symbol) in enumerate(subspace_info.outer_labels_symbols)
            if index_sym in subspace_info.inner_labels_symbols[outer_ind] || index_sym == outer_labels_symbol
                subsystem_ind = outer_ind
                if !(subsystem_ind in subspace_info.where_ensembles)
                    error("The index $index_sym is a subsystem index, but not an Ensemble subsystem index.")
                end
            end
        end
        if subsystem_ind === nothing
            ensemble_indexes = subspace_info.inner_labels[subspace_info.where_ensembles]
            error("Index $index_sym is not affiliated with a subsystem. The defined ensemble subsystems have the indexes $ensemble_indexes.")
        end
        outer_subsystem_inds[i] = subsystem_ind
    end

    if !contiguous_blocks(outer_subsystem_inds)
        error("Ensemble indexes must be contiguous. Indexes belonging to the same ensemble must be grouped.")
    end

    unique_outers = unique(filter(x -> x > 0, outer_subsystem_inds))
    group.ensemble_outer_indices = copy(unique_outers)
    group.index_outer_subspaces = copy(outer_subsystem_inds)
    group.sample_sizes = Int[]
    presence = BitVector(undef, length(subspace_info.where_ensembles))
    fill!(presence, false)
    for (pos, outer_idx) in enumerate(subspace_info.where_ensembles)
        presence[pos] = outer_idx in unique_outers
    end
    group.ensemble_presence = presence
    return outer_subsystem_inds, unique_outers
end

function register_group_with_subspaces!(group::ParameterGroupLike, subspaces::Vector{SubSpace})
    empty!(group.ensemble_subspaces)
    for outer_ind in group.ensemble_outer_indices
        outer_ind > 0 || continue
        subspace = subspaces[outer_ind]
        if !(subspace in group.ensemble_subspaces)
            push!(group.ensemble_subspaces, subspace)
        end
        ensemble_cfg = subspaces[outer_ind].ensemble
        if ensemble_cfg !== nothing && !(group.name in ensemble_cfg.param_groups)
            push!(ensemble_cfg.param_groups, group.name)
        end
    end
end

function build_indexed_parameters!(parameters::Vector{Parameter},
                                   group::ParameterGroupLike, group_index::Int,
                                   outer_subsystem_inds::Vector{Int}, subspace_info::SubSpaceInfo,
                                   t_indices::Vector{Int})::Vector{Array{Int}}
    param_name = String(group.name)
    blocks, block_lengths = find_blocks(outer_subsystem_inds)
    inner_label_symbols = subspace_info.inner_labels_symbols[blocks]
    ensemble_lengths = [length(inner_label_symbols[i]) for i in 1:length(inner_label_symbols) for _ in 1:block_lengths[i]]
    block_combinations = [collect(Combinatorics.with_replacement_combinations(1:ensemble_len, block_len))
                          for (ensemble_len, block_len) in zip(ensemble_lengths, block_lengths)]

    index_map_vec::Vector{Array{Int}} = [zeros(Int, ensemble_lengths...) for _ in t_indices]
    for comb_comb in Iterators.product(block_combinations...)
        inner_subspace_inds = vcat(comb_comb...)
        symbol_comb = vcat([inner_labels[comb] for (inner_labels, comb) in zip(inner_label_symbols, comb_comb)]...)
        str_comb = string.(symbol_comb)
        var_name_str, var_name_latex = symbol2formatted(param_name, str_comb)
        curr_var_name = unformatted_var_name(param_name, str_comb)
        subspace_indexes = [SubSpaceIndex(outer, inner, subspace_info)
                            for (outer, inner) in zip(outer_subsystem_inds, inner_subspace_inds)]
        for t_ind in t_indices
            t_suff = group.of_t ? "(" * t_suffix(t_ind) * ")" : ""
            t_suff_latex = group.of_t ? "(" * t_suffix(t_ind, do_latex=true) * ")" : ""
            push!(parameters, Parameter(group.name, curr_var_name * t_suff, var_name_str * t_suff,
                                        curr_var_name, var_name_str, var_name_latex * t_suff_latex, symbol_comb,
                                        group.of_t, false, t_ind, group_index, true, subspace_indexes,
                                        Int[], Vector{BitVector}(), Tuple{Int,Int}[], nothing, 0))
            index_map_vec[t_ind + 1][inner_subspace_inds...] = length(parameters)
        end
    end
    return index_map_vec
end

function build_scalar_parameters!(parameters::Vector{Parameter},
                                  group::ParameterGroupLike, group_index::Int,
                                  t_indices::Vector{Int})::Vector{Int}
    param_name = String(group.name)
    var_name_str, var_name_latex = symbol2formatted(param_name)
    t_index_map_vec = Vector{Int}(undef, length(t_indices))
    for t_ind in t_indices
        is_t = (group.name == :t)
        t_suff = ""
        t_suff_latex = ""
        if group.of_t && !is_t
            t_suff = "(" * t_suffix(t_ind) * ")"
            t_suff_latex = "(" * t_suffix(t_ind, do_latex=true) * ")"
        elseif is_t && t_ind > 0
            t_suff = str2sub(string(t_ind))
            t_suff_latex = "_$t_ind"
        end
        push!(parameters, Parameter(group.name, param_name * t_suff, var_name_str * t_suff,
                                    param_name, var_name_str, var_name_latex * t_suff_latex, Symbol[],
                                    group.of_t, is_t, t_ind, group_index, false, SubSpaceIndex[],
                                    Int[], Vector{BitVector}(), Tuple{Int,Int}[], nothing, 0))
        t_index_map_vec[t_ind + 1] = length(parameters)
    end
    return t_index_map_vec
end

function parameter_group2params!(parameters::Vector{Parameter},
                                 group::ParameterGroupLike, group_index::Int,
                                 subspace_info::SubSpaceInfo, subspaces::Vector{SubSpace},
                                 max_t_ind::Int)
    validate_group_payload(group)

    t_indices = group.of_t ? collect(0:max_t_ind) : [0]
    outer_subsystem_inds, unique_outers = resolve_ensemble_metadata!(group, subspace_info)

    if !isempty(unique_outers)
        if length(unique_outers) > 1 && !(group.kind in (ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction))
            name = String(group.name)
            error("Parameter \"$name\" spans multiple ensemble subspaces and therefore must be declared as an ensemble function (e.g. \"$name(t, ...)\" => args -> ...).")
        end
        if group.kind == ParameterGroupDistribution && group.payload isa QEnsembleFunction
            name = String(group.name)
            error("Parameter \"$name\" expects a QDistribution but received an ensemble function definition.")
        elseif (group.kind == ParameterGroupEnsembleFunction || group.kind == ParameterGroupEnsembleTimeFunction) && group.payload isa QDistribution
            name = String(group.name)
            error("Parameter \"$name\" expects an ensemble function but received a QDistribution.")
        end
    end

    register_group_with_subspaces!(group, subspaces)

    if !isempty(group.indexes)
        ensemble_maps = build_indexed_parameters!(parameters, group, group_index,
                                                  outer_subsystem_inds, subspace_info, t_indices)
        return ensemble_maps, Vector{Int}()
    else
        t_maps = build_scalar_parameters!(parameters, group, group_index, t_indices)
        ensemble_placeholder = [zeros(Int, 0) for _ in t_indices]
        return ensemble_placeholder, t_maps
    end
end

function build_t_index_transform(parameters::Vector{Parameter},
                                 ensemble_index_maps,
                                 t_index_maps,
                                 max_t_ind::Int)
    N = length(parameters)
    T = max_t_ind + 1
    transform = Array{SparsePermutation}(undef, T, T)
    for t_to in 0:max_t_ind
        for t_from in 0:max_t_ind
            idxs = Int[]
            diffs = Int[]
            @inbounds for (idx, param) in enumerate(parameters)
                dest = idx
                if param.param_of_t && param.t_index == t_from
                    g = param.group_index
                    if param.indexed_param
                        inners = [ix.inner for ix in param.param_indexes]
                        dest = ensemble_index_maps[g][t_to + 1][inners...]
                    else
                        dest = t_index_maps[g][t_to + 1]
                    end
                end
                if dest != idx
                    push!(idxs, idx)
                    push!(diffs, dest - idx)
                end
            end
            transform[t_to + 1, t_from + 1] = SparsePermutation(N, sparsevec(idxs, diffs, N))
        end
    end
    return transform
end

function collect_parameter_symbols(parameters::Vector{Parameter})
    return [param.param_symbol for param in parameters]
end

function collect_outer_label_strings(outer_labels_symbols::Vector{Symbol})
    outer_labels = String.(outer_labels_symbols)
    outer_labels_str = Vector{String}(undef, length(outer_labels_symbols))
    outer_labels_latex = Vector{String}(undef, length(outer_labels_symbols))
    for (i, sym) in enumerate(outer_labels_symbols)
        base_str, base_latex = symbol2formatted(String(sym))
        outer_labels_str[i] = base_str
        outer_labels_latex[i] = base_latex
    end
    return outer_labels, outer_labels_str, outer_labels_latex
end

function decorate_parameters!(parameters::Vector{Parameter}, subspace_info::SubSpaceInfo)
    ensemble_sizes = subspace_info.how_many_by_ensemble
    slot_counter = 0
    for param in parameters
        coord_length = 1 + length(param.param_indexes)
        coords = Vector{Int}(undef, coord_length)
        coords[1] = param.param_of_t ? param.t_index + 1 : 1
        acts = [falses(n) for n in ensemble_sizes]
        tuples = Tuple{Int,Int}[]
        for (inner_pos, sub_idx) in enumerate(param.param_indexes)
            coords[inner_pos + 1] = sub_idx.inner
            ensemble = subspace_info.ensemble_index_by_outer_index[sub_idx.outer]
            if ensemble != 0
                acts[ensemble][sub_idx.inner] = true
                push!(tuples, (ensemble, sub_idx.inner))
            end
        end
        param.coords = coords
        param.acts_on = acts
        param.ensemble_tuples = tuples
        param.function_refs = nothing
        if param.indexed_param
            slot_counter += 1
            param.indexed_slot = slot_counter
        else
            param.indexed_slot = 0
        end
    end
end

function ParameterIndexes(subspace_info::SubSpaceInfo, parameters::Vector{Parameter}, indexes_by_t_index::Vector{Vector{Int}})::ParameterIndexes
    labels::Vector{String} = []
    label_map = Dict{Tuple{Int,Int},Int}()
    global_idx = 1
    for outer_idx in subspace_info.where_ensembles
        inner_labels = subspace_info.inner_labels[outer_idx]
        for (inner_pos, label) in enumerate(inner_labels)
            push!(labels, label)
            label_map[(outer_idx, inner_pos)] = global_idx
            global_idx += 1
        end
    end
    label_parameter_indexes = [Int[] for _ in 1:length(labels)]
    for (param_idx, param) in enumerate(parameters)
        for sub_idx in param.param_indexes
            ensemble = subspace_info.ensemble_index_by_outer_index[sub_idx.outer]
            ensemble == 0 && continue
            label_idx = label_map[(sub_idx.outer, sub_idx.inner)]
            push!(label_parameter_indexes[label_idx], param_idx)
        end
    end
    t_labels = String[]
    t_labels_latex = String[]
    for (i, _) in enumerate(indexes_by_t_index)
        push!(t_labels, t_suffix(i-1, do_latex=false))
        push!(t_labels_latex, t_suffix(i-1, do_latex=true))
    end
    return ParameterIndexes(labels, t_labels, t_labels_latex, label_parameter_indexes, indexes_by_t_index)
end

function assign_params_to_groups!(groups::Vector{ParameterGroupLike}, parameters::Vector{Parameter})
    group_count = length(groups)
    params_by_group = [Int[] for _ in 1:group_count]
    size_buffers = [isempty(group.indexes) ? Int[] : zeros(Int, length(group.indexes)) for group in groups]

    for (idx, param) in enumerate(parameters)
        g = param.group_index
        push!(params_by_group[g], idx)
        coords = param.coords
        t_coord = coords[1]
        groups[g].time_count = max(groups[g].time_count, t_coord)
        buffer = size_buffers[g]
        if !isempty(buffer)
            for dim in 1:length(buffer)
                buffer[dim] = max(buffer[dim], coords[dim + 1])
            end
        end
    end

    for g in 1:group_count
        groups[g].parameter_indices = copy(params_by_group[g])
        if isempty(groups[g].indexes)
            groups[g].index_sizes = Int[]
        else
            groups[g].index_sizes = copy(size_buffers[g])
        end
    end

    return params_by_group
end

function build_time_param_lookup!(groups::Vector{ParameterGroupLike}, parameters::Vector{Parameter})
    lookup = Dict{Int,Int}()
    for (idx, param) in enumerate(parameters)
        if param.is_t
            lookup[param.t_index] = idx
            groups[param.group_index].is_time_group = true
        end
    end
    return lookup
end

function resolve_function_arguments!(groups::Vector{ParameterGroupLike},
                                     parameters::Vector{Parameter},
                                     params_by_group::Vector{Vector{Int}},
                                     time_param_lookup::Dict{Int,Int})
    group_count = length(groups)

    for g in 1:group_count
        group = groups[g]
        kind = group.kind
        payload = group.payload
        args = group.function_args
        indexes = group.indexes
        if !(kind in (ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction, ParameterGroupTimeFunction))
            continue
        end

        for param_idx in params_by_group[g]
            param = parameters[param_idx]
            refs = Vector{Int}(undef, length(args))
            main_index_map = Dict{String,Int}()
            required_index_tokens = Set(indexes)
            used_index_tokens = Set{String}()
            for (name, sub_idx) in zip(indexes, param.param_indexes)
                main_index_map[name] = sub_idx.inner
            end
            for (arg_pos, arg_str) in enumerate(args)
                if arg_str == "t"
                    t_key = param.coords[1] - 1
                    ref_param_idx = get(time_param_lookup, t_key) do
                        error("No time parameter found for t$(t_key) when evaluating ensemble function for $(param.param_name).")
                    end
                    refs[arg_pos] = ref_param_idx
                    if payload isa QEnsembleFunction
                        payload.argument_group_indices[arg_pos] = parameters[ref_param_idx].group_index
                    end
                    continue
                end
                arg_name, arg_tokens = underscore_separate(arg_str)
                target_group_idx = _find_group_index(groups, arg_name)
                target_group_idx === nothing &&
                    error("Unknown parameter group $arg_name referenced in ensemble function for $(param.param_name).")
                target_group_idx = target_group_idx::Int
                target_def = groups[target_group_idx]
                if kind == ParameterGroupEnsembleFunction || kind == ParameterGroupEnsembleTimeFunction
                    groups[target_group_idx].kind == ParameterGroupDistribution ||
                        error("Ensemble function for $(param.param_name) must reference distribution arguments; \"$(String(target_def.name))\" is not distribution-backed.")
                end
                if isempty(target_def.indexes)
                    !isempty(arg_tokens) && error("Argument $arg_str should not specify indexes for scalar parameter group $(String(target_def.name)).")
                else
                    length(arg_tokens) == length(target_def.indexes) || error("Argument $arg_str must specify $(length(target_def.indexes)) index tokens for parameter group $(String(target_def.name)).")
                end
                target_coords = Vector{Int}(undef, 1 + length(target_def.indexes))
                target_coords[1] = target_def.of_t ? param.coords[1] : 1
                for (tok_idx, tok) in enumerate(arg_tokens)
                    inner_val = get(main_index_map, tok) do
                        error("Index token $tok referenced in ensemble function for $(param.param_name) is undefined.")
                    end
                    push!(used_index_tokens, tok)
                    target_coords[tok_idx + 1] = inner_val
                end
                ref_idx = nothing
                for candidate_idx in params_by_group[target_group_idx]
                    if parameters[candidate_idx].coords == target_coords
                        ref_idx = candidate_idx
                        break
                    end
                end
                ref_idx === nothing && error("Unable to locate parameter for $arg_str in group $(String(target_def.name)) when evaluating $(param.param_name).")
                refs[arg_pos] = ref_idx
                if payload isa QEnsembleFunction
                    payload.argument_group_indices[arg_pos] = target_group_idx
                    if !isempty(target_def.indexes)
                        example_idx = first(params_by_group[target_group_idx])
                        example_param = parameters[example_idx]
                        self_positions = payload.argument_self_index_positions[arg_pos]
                        length(self_positions) == length(target_def.indexes) ||
                            error("Argument $arg_str provides $(length(self_positions)) index tokens but parameter group $(String(target_def.name)) declares $(length(target_def.indexes)).")
                        for (tok_idx, self_pos) in enumerate(self_positions)
                            parent_index = param.param_indexes[self_pos]
                            target_index = example_param.param_indexes[tok_idx]
                            parent_index.outer == target_index.outer ||
                                error("Argument $arg_str references index belonging to ensemble $(parent_index.outer), but parameter group $(target_def.name) expects ensemble $(target_index.outer).")
                        end
                    end
                end
            end
            if !isempty(required_index_tokens) && !(required_index_tokens ⊆ used_index_tokens)
                _missing = setdiff(required_index_tokens, used_index_tokens)
                signature = String(group.name)
                if !isempty(indexes)
                    signature *= "_{" * join(indexes, ",") * "}"
                end
                if group.of_t
                    signature *= "(t)"
                end
                message_missing = join(collect(_missing), ", ")
                error("Ensemble function for $(signature) does not reference index(es) $(message_missing). Each underscore index in $(signature) must appear in the function arguments.")
            end
            param.function_refs = refs
        end
    end
end

function finalize_group_dependencies!(groups::Vector{ParameterGroupLike})
    for (g, group) in enumerate(groups)
        deps = Int[]
        for dep_name in group.dependency_names
            dep_idx = _find_group_index(groups, dep_name)
            dep_idx === nothing && error("Parameter group $(String(group.name)) depends on undefined group \"$dep_name\".")
            push!(deps, dep_idx::Int)
        end
        group.dependency_indices = deps
        if !group.of_t
            group.time_count = max(group.time_count, 1)
        end
    end
    max_outer_idx = 0
    for group in groups
        if group.kind == ParameterGroupDistribution
            for outer_idx in group.ensemble_outer_indices
                outer_idx > 0 || continue
                max_outer_idx = max(max_outer_idx, outer_idx)
            end
        end
    end
    max_outer_idx == 0 && return
    dist_memberships = [Int[] for _ in 1:max_outer_idx]
    for (idx, group) in enumerate(groups)
        if group.kind == ParameterGroupDistribution
            for outer_idx in group.ensemble_outer_indices
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

function ParameterDefinitions2Parameters(vd::ParameterDefinitions, subspace_info::SubSpaceInfo,
                                         subspaces::Vector{SubSpace}, used_symbols::Set{Symbol},
                                         max_t_ind::Int)::Tuple{Vector{Parameter}, ParameterInfo, ParameterValues, ParameterDicts}
    groups = deepcopy(vd.groups)
    ensure_time_group!(groups, used_symbols)

    group_count = length(groups)

    outer_labels_symbols::Vector{Symbol} = Symbol[]
    parameters::Vector{Parameter} = Parameter[]
    ensemble_index_maps = Vector{Vector{Array{Int}}}()
    t_index_maps = Vector{Vector{Int}}()

    for (group_index, group) in enumerate(groups)
        param_name_sym = group.name
        push!(outer_labels_symbols, param_name_sym)
        if param_name_sym in used_symbols && param_name_sym != :t
            error("Variable name $param_name_sym is already used in the system! Please choose a distinct name.")
        end

        ensemble_maps, t_maps = parameter_group2params!(parameters, group, group_index,
                                                        subspace_info, subspaces, max_t_ind)
        push!(ensemble_index_maps, ensemble_maps)
        push!(t_index_maps, t_maps)

        push!(used_symbols, param_name_sym)
    end

    t_index_transform = build_t_index_transform(parameters, ensemble_index_maps, t_index_maps, max_t_ind)
    subspace_index_maps = build_subspace_index_maps(parameters, ensemble_index_maps, subspace_info)

    inner_labels_symbols_flat = collect_parameter_symbols(parameters)
    outer_labels, outer_labels_str, outer_labels_latex = collect_outer_label_strings(outer_labels_symbols)
    indexes_of_t = Int[]
    max_t = -1
    for (idx, param) in enumerate(parameters)
        if param.param_of_t
            push!(indexes_of_t, idx)
            max_t = max(max_t, param.t_index)
        end
    end
    indexes_by_t_index = max_t < 0 ? Vector{Vector{Int}}() : [Int[] for _ in 0:max_t]
    if max_t >= 0
        for idx in indexes_of_t
            param = parameters[idx]
            push!(indexes_by_t_index[param.t_index + 1], idx)
        end
    end
    ensemble_sizes = subspace_info.how_many_by_ensemble
    decorate_parameters!(parameters, subspace_info)
    param_indexes = ParameterIndexes(subspace_info, parameters, indexes_by_t_index)
    params_by_group = assign_params_to_groups!(groups, parameters)
    time_param_lookup = build_time_param_lookup!(groups, parameters)
    resolve_function_arguments!(groups, parameters, params_by_group, time_param_lookup)
    finalize_group_dependencies!(groups)
    param_groups = groups
    normalize_ensemble_group_lists!(subspaces, param_groups)

    param_info = ParameterInfo(outer_labels_symbols, inner_labels_symbols_flat, outer_labels, outer_labels_str, outer_labels_latex,
        subspace_index_maps, t_index_transform, indexes_by_t_index, indexes_of_t, ensemble_sizes,
        subspace_info, param_indexes, param_groups, Any[parameters...])

    param_dicts = build_parameter_dicts(param_info)
    sample_index_param_values = ParameterValues(param_groups)

    return parameters, param_info, sample_index_param_values, param_dicts
end

#Return parameter mapping vector for switching ensemble inner index within a subspace.
function map_by_subspace(i_to::SubSpaceIndex, i_from::SubSpaceIndex, pinfo::ParameterInfo)::Vector{Int}
    @assert i_to.outer == i_from.outer "Subspace mapping requires the same outer subspace."
    M = pinfo.subspace_index_maps[i_to.outer]
    if size(M,1) == 0  # non-ensemble subspace -> identity map
        return collect(1:length(pinfo.parameters))
    else
        return denseperm(M[i_to.inner, i_from.inner])
    end
end

# Return parameter mapping vector for switching from t_index2 to t_index1.
map_by_tindex(t_index1::Int, t_index2::Int, pinfo::ParameterInfo) = denseperm(pinfo.t_index_transform[t_index1+1, t_index2+1])
