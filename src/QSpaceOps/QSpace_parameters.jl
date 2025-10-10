using Combinatorics
using SparseArrays
using ..CFunctions: ParameterInfo, ParameterIndexes, ParameterDicts, ParameterValues, build_parameter_dicts
using ..StringUtils: symbol2formatted, str2sub, unformat_symbol, var_unsubstitution, reverse_var_substitution
using ..SparsePermutationTools: SparsePermutation, denseperm
using ..Sampler: QDistribution, QEnsembleFunction, AbstractEnsembleSample
using ..ParameterGroups: ParameterGroupKind, ParameterGroup,
                          ParameterGroupScalar, ParameterGroupTimeScalar, ParameterGroupTimeFunction,
                          ParameterGroupDistribution, ParameterGroupEnsembleFunction, WhereWhichParamGroup

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
end

mutable struct ParameterGroupDefinition
    name::String
    of_t::Bool
    indexes::Vector{String}
    function_args::Vector{String}
    payload::Any
    kind::ParameterGroupKind
    dependencies::Vector{String}
    display_signature::String
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

@inline function _infer_group_kind(name::String, of_t::Bool, indexes, function_args::Vector{String})::ParameterGroupKind
    has_indexes = !isempty(indexes)
    has_args = !isempty(function_args)
    if name == "t"
        return ParameterGroupTimeScalar
    elseif has_indexes
        return has_args ? ParameterGroupEnsembleFunction : ParameterGroupDistribution
    elseif has_args || (of_t && name != "t")
        return ParameterGroupTimeFunction
    else
        return ParameterGroupScalar
    end
end

function _assign_group_definition!(group::ParameterGroupDefinition, payload)
    payload === nothing && return group
    kind = group.kind
    if kind in (ParameterGroupScalar, ParameterGroupTimeScalar)
        payload isa Function && error("Parameter group $(group.name) expects a literal value, not a function.")
        group.payload = payload
    elseif kind == ParameterGroupTimeFunction
        payload isa Function ||
            error("Parameter group $(group.name) expects a Function definition.")
        group.payload = payload
    elseif kind == ParameterGroupDistribution
        payload isa QDistribution ||
            error("Parameter group $(group.name) expects a QDistribution, got $(typeof(payload)).")
        group.payload = payload
    elseif kind == ParameterGroupEnsembleFunction
        if payload isa QEnsembleFunction
            expected = Symbol.(group.function_args)
            payload.argument_symbols == expected ||
                error("Ensemble function for $(group.name) expects arguments $(expected), got $(payload.argument_symbols).")
            group.payload = QEnsembleFunction(group.name, group.indexes, group.function_args, payload.func)
        elseif payload isa Function
            group.payload = _build_qensemble_function(group.name, group.indexes, group.function_args, payload)
        else
            error("Parameter group $(group.name) expects a QEnsembleFunction or plain Function, got $(typeof(payload)).")
        end
    else
        error("Unsupported parameter group kind $(kind).")
    end
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
    var_param::Vector{ParameterGroupDefinition}
    function ParameterDefinitions(params...)
    var_param = ParameterGroupDefinition[]
    for var in params
        label, payload = _extract_group_payload(var)
        label_clean = reverse_var_substitution(label)
        pre, brace_elements = brace_separate(label_clean)
            name, indexes = underscore_separate(pre)
            of_t = (name == "t") || ("t" in brace_elements)
            function_args = [String(arg) for arg in brace_elements]
            index_tokens = [String(ix) for ix in indexes]
            inferred_kind = _infer_group_kind(name, of_t, index_tokens, function_args)
            display_signature = _group_display_signature(name, index_tokens, function_args)
            group_def = ParameterGroupDefinition(name, of_t, index_tokens, function_args, nothing, inferred_kind, String[], display_signature)
            payload === nothing || _assign_group_definition!(group_def, payload)
            push!(var_param, group_def)
    end

    _finalize_group_dependencies!(var_param)

        sort!(var_param, by = g -> ((g.kind == ParameterGroupTimeScalar && g.name == "t") ? -1 : Int(g.kind), g.name))

    return new(var_param)
end
end

function _finalize_group_dependencies!(group_defs::Vector{ParameterGroupDefinition})
    for group in group_defs
        deps = String[]
        if group.kind in (ParameterGroupEnsembleFunction, ParameterGroupTimeFunction)
            for arg in group.function_args
                if arg == "t"
                    push!(deps, "t")
                    continue
                end
                arg_name, _ = underscore_separate(arg)
                arg_group = _find_group_by_name(group_defs, arg_name)
                arg_group !== nothing ||
                    error("Function for parameter group $(group.name) references unknown group \"$arg_name\".")
                push!(deps, arg_name)
            end
        end
        if group.kind == ParameterGroupEnsembleFunction
            for arg in group.function_args
                arg == "t" && continue
                arg_name, _ = underscore_separate(arg)
                target = _find_group_by_name(group_defs, arg_name)
                target === nothing &&
                    error("Ensemble function for parameter group $(group.name) references unknown group \"$arg_name\".")
                if target.kind != ParameterGroupDistribution
                    error("Ensemble function for parameter group $(group.name) must reference distribution arguments; \"$arg_name\" is not distribution-backed.")
                end
            end
        end
        group.dependencies = deps
    end
end

function _find_group_by_name(groups::Vector{ParameterGroupDefinition}, name::String)
    for group in groups
        if group.name == name
            return group
        end
    end
    return nothing
end

@inline function _find_group_index(groups::Vector{ParameterGroupDefinition}, name::Union{String,Symbol})
    target = name isa Symbol ? String(name) : name
    @inbounds for (idx, group) in enumerate(groups)
        group.name == target && return idx
    end
    return nothing
end

function _ensure_signature_alignment!(group::ParameterGroupDefinition,
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
    group = _find_group_by_name(defs.var_param, name)
    group === nothing && error("Unknown parameter group \"$name\" in ParameterDefinitions.")
    _ensure_signature_alignment!(group, indexes, Vector{String}(brace_elements))
    _assign_group_definition!(group, payload)
    return group
end

set_parameter_group_definition!(defs::ParameterDefinitions, pair::Pair) =
    set_parameter_group_definition!(defs, pair[1], pair[2])


function Base.show(io::IO, param_def::ParameterDefinitions)
    var_str_vec = []
    for group_def in param_def.var_param
        param_str = symbol2formatted(group_def.name)[1]
        if group_def.of_t
            param_str *= "(t)"
        end
        param_str *= str2sub(join(group_def.indexes, ","))
        push!(var_str_vec, param_str) 
    end 
    println(io, "ParameterDefinitions: [" * join(var_str_vec, ", ") * "]")
end

function _build_qensemble_function(name::String, indexes::Vector{String}, brace_elements::Vector{String}, f::Function)
    isempty(brace_elements) && error("Parameter \"$name\" requires a parentheses list specifying argument order when providing an ensemble function, e.g. \"$name(t, alpha)\" => (t, alpha) -> ...")
    return QEnsembleFunction(name, indexes, brace_elements, f)
end

function ParameterIndexes(subspace_info::SubSpaceInfo, indexed_parameter_indexes::Vector{Int}, where_acting_by_parameter::Vector{Vector{BitVector}}, indexes_by_t_index::Vector{Vector{Int}})::ParameterIndexes
    labels::Vector{String} = []
    t_labels::Vector{String} = []
    t_labels_latex::Vector{String} = []
    for w in subspace_info.where_ensembles
        inner = subspace_info.inner_labels[w]
        append!(labels, inner) 
    end
    label_parameter_indexes::Vector{Vector{Int}} = [Int[] for _ in 1:length(labels)]
    for (param_idx, position) in pairs(indexed_parameter_indexes)
        position == 0 && continue
        param_acting = where_acting_by_parameter[position]
        flattened_acting = vcat(param_acting...)
        inds = findall(flattened_acting)
        for ind in inds 
            append!(label_parameter_indexes[ind], param_idx) 
        end
    end
    for (i, t_indexes) in enumerate(indexes_by_t_index)
        push!(t_labels, t_suffix(i-1, do_latex=false))
        push!(t_labels_latex, t_suffix(i-1, do_latex=true))
    end
    return ParameterIndexes(labels, t_labels, t_labels_latex, label_parameter_indexes, indexes_by_t_index)
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

function ParameterDefinitions2Parameters(vd::ParameterDefinitions, subspace_info::SubSpaceInfo,
                                         subspaces::Vector{SubSpace}, used_symbols::Set{Symbol},
                                         max_t_ind::Int)::Tuple{Vector{Parameter}, ParameterInfo, ParameterValues, ParameterDicts}
    # --- start from a local copy and auto-add t if not present ---
    var_param = copy(vd.var_param)
    if all(group.name != "t" for group in var_param) && !(:t in used_symbols)
        push!(var_param, ParameterGroupDefinition("t", true, String[], String[], nothing,
                                                  ParameterGroupTimeScalar, String[], _group_display_signature("t", String[], String[])))
    end

    group_defs = var_param
    group_count = length(group_defs)

    group_names = Vector{Symbol}(undef, group_count)
    group_display_signatures = Vector{String}(undef, group_count)
    group_kinds = Vector{ParameterGroupKind}(undef, group_count)
    group_payloads = Vector{Any}(undef, group_count)
    group_of_t = BitVector(undef, group_count)
    group_indexes = Vector{Vector{String}}(undef, group_count)
    group_function_args = Vector{Vector{String}}(undef, group_count)
    group_dependency_names = Vector{Vector{String}}(undef, group_count)
    group_outer_indices = Vector{Vector{Int}}(undef, group_count)
    group_index_outer_subspaces = Vector{Vector{Int}}(undef, group_count)
    group_ensemble_presence = Vector{BitVector}(undef, group_count)
    group_sample_sizes = Vector{Vector{Int}}(undef, group_count)

    outer_labels_symbols::Vector{Symbol} = Symbol[]
    parameters::Vector{Parameter} = Parameter[]
    ensemble_index_maps::Vector{Vector{Array{Int}}} = Vector{Vector{Array{Int}}}()
    t_index_maps::Vector{Vector{Int}} = Vector{Vector{Int}}()
    param_of_indexes::BitVector = Bool[]

    for (group_index, group_def) in enumerate(group_defs)
        param_name = group_def.name
        of_t = group_def.of_t
        index_strs = group_def.indexes
        payload = group_def.payload
        kind = group_def.kind

        group_names[group_index] = Symbol(param_name)
        group_display_signatures[group_index] = group_def.display_signature
        group_kinds[group_index] = kind
        group_payloads[group_index] = payload
        group_of_t[group_index] = of_t
        group_indexes[group_index] = copy(index_strs)
        group_function_args[group_index] = copy(group_def.function_args)
        group_dependency_names[group_index] = copy(group_def.dependencies)

        if kind in (ParameterGroupScalar, ParameterGroupTimeScalar) && payload isa Function
            error("Parameter group $(param_name) expects a literal value, not a function.")
        elseif kind == ParameterGroupTimeFunction && !(payload === nothing || payload isa Function)
            error("Parameter group $(param_name) expects a Function definition.")
        elseif kind == ParameterGroupDistribution && !(payload === nothing || payload isa QDistribution)
            error("Parameter group $(param_name) expects a QDistribution, got $(typeof(payload)).")
        elseif kind == ParameterGroupEnsembleFunction && !(payload === nothing || payload isa QEnsembleFunction)
            error("Parameter group $(param_name) expects a QEnsembleFunction, got $(typeof(payload)).")
        end

        var_name_sym::Symbol = Symbol(param_name)
        push!(outer_labels_symbols, var_name_sym)
        if var_name_sym in used_symbols && var_name_sym != :t
            error("Variable name $var_name_sym is already used in the system! Please choose a distinct name.")
        end

        t_vals = of_t ? collect(0:max_t_ind) : [0]

        # resolve ensemble indexes
        index_str_syms = Symbol.(index_strs)
        outer_subsystem_inds = -ones(Int, length(index_str_syms))
        for (i, index_sym) in enumerate(index_str_syms)
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

        belongs_to_ensemble = !isempty(index_str_syms)
        unique_outers = belongs_to_ensemble ? unique(filter(x->x>0, outer_subsystem_inds)) : Int[]
        group_outer_indices[group_index] = unique_outers
        group_index_outer_subspaces[group_index] = copy(outer_subsystem_inds)
        group_sample_sizes[group_index] = Int[]
        presence = BitVector(undef, length(subspace_info.where_ensembles))
        fill!(presence, false)
        for (pos, outer_idx) in enumerate(subspace_info.where_ensembles)
            presence[pos] = outer_idx in unique_outers
        end
        group_ensemble_presence[group_index] = presence

        if belongs_to_ensemble
            if length(unique_outers) > 1 && kind != ParameterGroupEnsembleFunction
                error("Parameter \"$param_name\" spans multiple ensemble subspaces and therefore must be declared as an ensemble function (e.g. \"$param_name(t, ...)\" => args -> ...).")
            end
            if kind == ParameterGroupDistribution && payload isa QEnsembleFunction
                error("Parameter \"$param_name\" expects a QDistribution but received an ensemble function definition.")
            elseif kind == ParameterGroupEnsembleFunction && payload isa QDistribution
                error("Parameter \"$param_name\" expects an ensemble function but received a QDistribution.")
            end
        end

        if !isempty(unique_outers)
            param_sym = Symbol(param_name)
            for outer_ind in unique_outers
                ensemble_cfg = subspaces[outer_ind].ensemble
                if ensemble_cfg !== nothing && !(param_sym in ensemble_cfg.param_groups)
                    push!(ensemble_cfg.param_groups, param_sym)
                end
            end
        end

        blocks, block_lengths = find_blocks(outer_subsystem_inds)
        inner_label_symbols = subspace_info.inner_labels_symbols[blocks]
        ensemble_lengths = [length(inner_label_symbols[i]) for i in 1:length(inner_label_symbols) for _ in 1:block_lengths[i]]
        block_combinations = [collect(Combinatorics.with_replacement_combinations(1:ensemble_len, block_len))
                              for (ensemble_len, block_len) in zip(ensemble_lengths, block_lengths)]

        if !isempty(block_combinations)
            # -------- indexed group: real ensemble map, placeholder t-map --------
            index_map_vec::Vector{Array{Int}} = [zeros(Int, ensemble_lengths...) for _ in t_vals]

            for comb_comb in Iterators.product(block_combinations...)
                inner_subspace_inds = vcat(comb_comb...)
                symbol_comb = vcat([inner_labels[comb] for (inner_labels, comb) in zip(inner_label_symbols, comb_comb)]...)
                str_comb = string.(symbol_comb)
                var_name_str, var_name_latex = symbol2formatted(param_name, str_comb)
                curr_var_name = unformatted_var_name(param_name, str_comb)
                subspace_indexes::Vector{SubSpaceIndex} = [SubSpaceIndex(outer, inner, subspace_info)
                                                    for (outer, inner) in zip(outer_subsystem_inds, inner_subspace_inds)]
                for t_ind in t_vals
                    t_suff        = of_t ? "(" * t_suffix(t_ind) * ")" : ""
                    t_suff_latex  = of_t ? "(" * t_suffix(t_ind, do_latex=true) * ")" : ""
                    push!(parameters, Parameter(Symbol(param_name), curr_var_name*t_suff, var_name_str*t_suff,
                                                curr_var_name, var_name_str, var_name_latex*t_suff_latex, symbol_comb,
                                                of_t, false, t_ind, group_index, true, subspace_indexes))
                    index_map_vec[t_ind+1][inner_subspace_inds...] = length(parameters)
                    push!(param_of_indexes, true)
                end
            end
            push!(ensemble_index_maps, index_map_vec)         # real
            push!(t_index_maps, Vector{Int}())                # placeholder to keep g-alignment

        else
            # -------- non-indexed group: real t-map, placeholder ensemble map --------
            var_name_str, var_name_latex = symbol2formatted(param_name)
            t_index_map_vec::Vector{Int} = Vector{Int}(undef, length(t_vals))

            for t_ind in t_vals
                is_t = (var_name_sym==:t)
                t_suff = ""
                t_suff_latex = ""

                if of_t && !is_t
                    # e.g. beta(t)
                    t_suff       = "(" * t_suffix(t_ind) * ")"
                    t_suff_latex = "(" * t_suffix(t_ind, do_latex=true) * ")"
                elseif is_t && t_ind > 0
                    # the actual variable t_1, t_2, ...
                    t_suff       = str2sub(string(t_ind))
                    t_suff_latex = "_$t_ind"
                end
                push!(parameters, Parameter(Symbol(param_name), param_name*t_suff, var_name_str*t_suff,
                                            param_name, var_name_str, var_name_latex*t_suff_latex, Symbol[],
                                            of_t, is_t, t_ind, group_index, false, SubSpaceIndex[]))
                t_index_map_vec[t_ind+1] = length(parameters)
                push!(param_of_indexes, false)
            end

            push!(t_index_maps, t_index_map_vec)                              # real
            push!(ensemble_index_maps, [zeros(Int, 0) for _ in t_vals])       # placeholder to keep g-alignment
        end

        push!(used_symbols, Symbol(param_name))
    end

    # ---- build t_index_transform (t_to, t_from) -> mapping vector over full parameter space ----
    N = length(parameters)
    T = max_t_ind + 1
    t_index_transform = Array{SparsePermutation}(undef, T, T)
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
                        dest = ensemble_index_maps[g][t_to+1][inners...]
                    else
                        dest = t_index_maps[g][t_to+1]
                    end
                end
                if dest != idx
                    push!(idxs, idx)
                    push!(diffs, dest - idx)
                end
            end
            t_index_transform[t_to+1, t_from+1] = SparsePermutation(N, sparsevec(idxs, diffs, N))
        end
    end

    subspace_index_maps = build_subspace_index_maps(parameters, ensemble_index_maps, subspace_info)

    inner_labels_symbols_flat = [param.param_symbol for param in parameters]
    param_names  = [param.param_name  for param in parameters]
    param_strs   = [param.param_str   for param in parameters]
    param_latex  = [param.param_latex for param in parameters]

    param_of_t::BitVector = [param.param_of_t for param in parameters]
    param_is_t::BitVector = [param.is_t       for param in parameters]

    outer_labels = String.(outer_labels_symbols)
    outer_labels_str = Vector{String}(undef, length(outer_labels_symbols))
    outer_labels_latex = Vector{String}(undef, length(outer_labels_symbols))
    for (i, sym) in enumerate(outer_labels_symbols)
        base_str, base_latex = symbol2formatted(String(sym))
        outer_labels_str[i] = base_str
        outer_labels_latex[i] = base_latex
    end
    param_group_by_index = zeros(Int, length(parameters))
    t_index_by_index = zeros(Int, length(parameters))
    indexes_of_t = Int[]

    for (i, param) in enumerate(parameters)
        param_group_by_index[i] = param.group_index
        t_index_by_index[i] = param.t_index - !param.param_of_t
        if param.param_of_t
            push!(indexes_of_t, i)
        end
    end
    indexes_by_t_index = [findall(==(t_ind), t_index_by_index) for t_ind in 0:maximum(t_index_by_index)]

    indexed_parameter_indexes = zeros(Int, length(parameters))
    where_acting_by_parameter = Vector{Vector{BitVector}}()
    ensemble_sizes = subspace_info.how_many_by_ensemble

    curr_ind = 1
    for (i, param) in enumerate(parameters)
        if param.indexed_param
            indexed_parameter_indexes[i] = curr_ind
            curr_ind += 1
            curr_bools = [falses(n) for n in ensemble_sizes]
            for curr_ind in param.param_indexes
                outer = curr_ind.outer
                inner = curr_ind.inner
                outer_ind = subspace_info.ensemble_index_by_outer_index[outer]
                curr_bools[outer_ind][inner] = true
            end
            push!(where_acting_by_parameter, curr_bools)
        end
    end
    params_acting_by_index = [[falses(length(parameters)) for _ in 1:n] for n in ensemble_sizes]
    for (param_idx, storage_idx) in pairs(indexed_parameter_indexes)
        storage_idx == 0 && continue
        param_acts = where_acting_by_parameter[storage_idx]
        for (ensemble_idx, bits) in enumerate(param_acts)
            for inner_idx in findall(bits)
                params_acting_by_index[ensemble_idx][inner_idx][param_idx] = true
            end
        end
    end

    param_indexes = ParameterIndexes(subspace_info, indexed_parameter_indexes, where_acting_by_parameter, indexes_by_t_index)

    param_index_tuples = Vector{Vector{Tuple{Int,Int}}}(undef, length(parameters))
    for (idx, param) in enumerate(parameters)
        if param.indexed_param
            tuples = Vector{Tuple{Int,Int}}(undef, length(param.param_indexes))
            @inbounds for (inner_pos, sub_idx) in enumerate(param.param_indexes)
                ensemble = subspace_info.ensemble_index_by_outer_index[sub_idx.outer]
                ensemble != 0 || error("Parameter index does not belong to an ensemble subspace.")
                tuples[inner_pos] = (ensemble, sub_idx.inner)
            end
            param_index_tuples[idx] = tuples
        else
            param_index_tuples[idx] = Tuple{Int,Int}[]
        end
    end

    params_by_group = [Int[] for _ in 1:group_count]
    group_time_counts = fill(1, group_count)
    param_coords = Vector{Vector{Int}}(undef, length(parameters))
    for (idx, param) in enumerate(parameters)
        g = param.group_index
        push!(params_by_group[g], idx)
        t_coord = param.param_of_t ? param.t_index + 1 : 1
        group_time_counts[g] = max(group_time_counts[g], t_coord)
        coords = Vector{Int}(undef, 1 + length(param.param_indexes))
        coords[1] = t_coord
        for (k, sub_idx) in enumerate(param.param_indexes)
            coords[k+1] = sub_idx.inner
        end
        param_coords[idx] = coords
    end

    group_index_sizes = [Int[] for _ in 1:group_count]
    for g in 1:group_count
        idx_names = group_indexes[g]
        if isempty(idx_names)
            group_index_sizes[g] = Int[]
        else
            size_vec = zeros(Int, length(idx_names))
            for idx in params_by_group[g]
                coords = param_coords[idx]
                for dim in 1:length(idx_names)
                    size_vec[dim] = max(size_vec[dim], coords[dim+1])
                end
            end
            group_index_sizes[g] = size_vec
        end
    end

    group_is_t = falses(group_count)
    time_param_lookup = Dict{Int,Int}()
    for (idx, param) in enumerate(parameters)
        if param.is_t
            time_param_lookup[param.t_index] = idx
            group_is_t[param.group_index] = true
        end
    end

    function_param_refs = Vector{Union{Nothing,Vector{Int}}}(undef, length(parameters))
    fill!(function_param_refs, nothing)
    for g in 1:group_count
        kind = group_kinds[g]
        payload = group_payloads[g]
        args = group_function_args[g]
        if !(kind in (ParameterGroupEnsembleFunction, ParameterGroupTimeFunction))
            continue
        end

        for idx in params_by_group[g]
            param = parameters[idx]
            refs = Vector{Int}(undef, length(args))
            main_index_map = Dict{String,Int}()
            required_index_tokens = Set(group_indexes[g])
            used_index_tokens = Set{String}()
            for (name, sub_idx) in zip(group_indexes[g], param.param_indexes)
                main_index_map[name] = sub_idx.inner
            end
            for (arg_pos, arg_str) in enumerate(args)
                if arg_str == "t"
                    t_key = param_coords[idx][1] - 1
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
                target_group_idx = _find_group_index(group_defs, arg_name)
                target_group_idx === nothing &&
                    error("Unknown parameter group $arg_name referenced in ensemble function for $(param.param_name).")
                target_group_idx = target_group_idx::Int
                target_def = group_defs[target_group_idx]
                if kind == ParameterGroupEnsembleFunction
                    group_kinds[target_group_idx] == ParameterGroupDistribution ||
                        error("Ensemble function for $(param.param_name) must reference distribution arguments; \"$(target_def.name)\" is not distribution-backed.")
                end
                if isempty(target_def.indexes)
                    !isempty(arg_tokens) && error("Argument $arg_str should not specify indexes for scalar parameter group $(target_def.name).")
                else
                    length(arg_tokens) == length(target_def.indexes) || error("Argument $arg_str must specify $(length(target_def.indexes)) index tokens for parameter group $(target_def.name).")
                end
                target_coords = Vector{Int}(undef, 1 + length(target_def.indexes))
                target_coords[1] = target_def.of_t ? param_coords[idx][1] : 1
                for (tok_idx, tok) in enumerate(arg_tokens)
                    inner_val = get(main_index_map, tok) do
                        error("Index token $tok referenced in ensemble function for $(param.param_name) is undefined.")
                    end
                    push!(used_index_tokens, tok)
                    target_coords[tok_idx+1] = inner_val
                end
                ref_idx = findfirst(j -> param_coords[j] == target_coords, params_by_group[target_group_idx])
                ref_idx === nothing && error("Unable to locate parameter for $arg_str in group $(target_def.name) when evaluating $(param.param_name).")
                refs[arg_pos] = params_by_group[target_group_idx][ref_idx]
                if payload isa QEnsembleFunction
                    payload.argument_group_indices[arg_pos] = target_group_idx
                    if !isempty(target_def.indexes)
                        example_idx = first(params_by_group[target_group_idx])
                        example_param = parameters[example_idx]
                        self_positions = payload.argument_self_index_positions[arg_pos]
                        length(self_positions) == length(target_def.indexes) ||
                            error("Argument $arg_str provides $(length(self_positions)) index tokens but parameter group $(target_def.name) declares $(length(target_def.indexes)).")
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
                missing = setdiff(required_index_tokens, used_index_tokens)
                signature = group_defs[g].name
                if !isempty(group_indexes[g])
                    signature *= "_{" * join(group_indexes[g], ",") * "}"
                end
                if group_of_t[g]
                    signature *= "(t)"
                end
                message_missing = join(collect(missing), ", ")
                error("Ensemble function for $(signature) does not reference index(es) $(message_missing). Each underscore index in $(signature) must appear in the function arguments.")
            end
            function_param_refs[idx] = refs
        end
    end

    # compute dependency indices using resolved group map
    name_to_index = Dict{String,Int}()
    for (idx, def) in enumerate(group_defs)
        name_to_index[def.name] = idx
    end
    group_dependency_indices = Vector{Vector{Int}}(undef, group_count)
    for g in 1:group_count
        deps = Vector{Int}()
        for dep_name in group_dependency_names[g]
            dep_idx = get(name_to_index, dep_name, nothing)
            dep_idx === nothing && error("Parameter group $(group_defs[g].name) depends on undefined group \"$dep_name\".")
            push!(deps, dep_idx)
        end
        group_dependency_indices[g] = deps
    end

    param_groups = Vector{ParameterGroup}(undef, group_count)
    for g in 1:group_count
        param_groups[g] = ParameterGroup(
            group_names[g],
            group_display_signatures[g],
            group_kinds[g],
            group_of_t[g],
            group_indexes[g],
            group_function_args[g],
            group_dependency_names[g],
            group_dependency_indices[g],
            group_outer_indices[g],
            group_index_outer_subspaces[g],
            group_ensemble_presence[g],
            params_by_group[g],
            group_time_counts[g],
            group_index_sizes[g],
            group_sample_sizes[g],
            group_is_t[g],
            group_payloads[g],
        )
    end

    symbol_to_index = Dict{Symbol,Int}(group_names[i] => i for i in 1:group_count)
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
            group_idx = get(symbol_to_index, group_sym, nothing)
            group_idx === nothing && continue
            group = param_groups[group_idx]
            if group.kind == ParameterGroupDistribution
                group_idx in dist_idxs || push!(dist_idxs, group_idx)
            elseif group.kind == ParameterGroupEnsembleFunction
                group_idx in func_idxs || push!(func_idxs, group_idx)
            end
        end
        ens.param_groups = ordered_syms
        ens.distribution_group_indices = dist_idxs
        ens.ensemble_function_group_indices = func_idxs
    end

    param_info = CFunctions.ParameterInfo(outer_labels_symbols, inner_labels_symbols_flat, outer_labels, outer_labels_str, outer_labels_latex, param_names,
        param_strs, param_latex, param_of_indexes, param_group_by_index,
        t_index_by_index, indexed_parameter_indexes,
        where_acting_by_parameter, params_acting_by_index, param_index_tuples, subspace_index_maps, t_index_transform, indexes_by_t_index,
        indexes_of_t, ensemble_sizes, param_of_t, param_is_t, function_param_refs,
        param_coords, subspace_info, param_indexes, param_groups)

    param_dicts = build_parameter_dicts(param_info)
    sample_index_param_values = ParameterValues(param_info)

    return parameters, param_info, sample_index_param_values, param_dicts
end

#Return parameter mapping vector for switching ensemble inner index within a subspace.
function map_by_subspace(i_to::SubSpaceIndex, i_from::SubSpaceIndex, pinfo::ParameterInfo)::Vector{Int}
    @assert i_to.outer == i_from.outer "Subspace mapping requires the same outer subspace."
    M = pinfo.subspace_index_maps[i_to.outer]
    if size(M,1) == 0  # non-ensemble subspace -> identity map
        return collect(1:length(pinfo.param_group_by_index))
    else
        return denseperm(M[i_to.inner, i_from.inner])
    end
end

# Return parameter mapping vector for switching from t_index2 to t_index1.
map_by_tindex(t_index1::Int, t_index2::Int, pinfo::ParameterInfo) = denseperm(pinfo.t_index_transform[t_index1+1, t_index2+1])
# from t_index2 to t_index1 
