using Combinatorics
using SparseArrays
using ..CFunctions: ParameterInfo, ParameterIndexes, ParameterDicts, ParameterValues
using ..StringUtils: symbol2formatted, str2sub
using ..SparsePermutationTools: SparsePermutation, denseperm
using ..QDistributions: QDistribution, QEnsembleFunction

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
    while occursin("__", cleaned)
        cleaned = replace(cleaned, "__" => "_")
    end
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

""" 
    Parameter(param_name::String, param_of_t::Bool, var_of_ensemble::Bool, var_ensemble_index::Int=0;
              var_suffix::String="")

Container describing a single parameter instance in the `QSpace`. Ensemble distributions
are tracked per parameter group and available through `ParameterValues.ensemble_group_distributions`.
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

struct ParameterGroupDefinition
    name::String
    of_t::Bool
    indexes::Vector{String}
    distribution::Union{Nothing,QDistribution}
    ensemble_function::Union{Nothing,QEnsembleFunction}
    scalar_function::Union{Nothing,Function}
    function_args::Vector{String}
end

"""
    ParameterDefinitions(params...)

Construct a parameter definition list from symbols or strings. Use `(t)` to mark
time-dependent groups. Ensemble parameters must be provided together with either a
`QDistribution` or a `QEnsembleFunction`, supplied as `(definition, payload)` tuples
or `definition => payload`. Constructing a parameter that references ensemble indexes
without one of these payloads throws an error. Distributions/functions are stored
once per parameter group and exposed via `qspace.param_values.ensemble_group_distributions`
and `qspace.param_values.ensemble_group_functions` after construction.

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
            pre, brace_elements = brace_separate(label)
            name, indexes = underscore_separate(pre)
            of_t = "t" in brace_elements

            dist, qfun, scalar_fun = _coerce_payload(name, brace_elements, payload, !isempty(indexes))

            if qfun === nothing
                extra_args = filter(x -> x != "t", brace_elements)
                if !isempty(extra_args)
                    error("Parameter \"$label\" lists arguments $(extra_args) but no ensemble function was provided. Supply one via \"$label\" => (args -> ...).")
                end
            else
                if !of_t && (:t in qfun.argument_symbols)
                    of_t = true
                end
            end

            if dist !== nothing && qfun !== nothing
                error("Parameter \"$name\" received both a QDistribution and a QEnsembleFunction. Provide only one.")
            end

            function_args = copy(brace_elements)
            push!(var_param, ParameterGroupDefinition(name, of_t, indexes, dist, qfun, scalar_fun, function_args))
        end
        return new(var_param)
    end
end
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

_to_param_string(name::String) = name
_to_param_string(name::Symbol) = String(name)
_to_param_string(name) = error("Unsupported parameter label type $(typeof(name)). Expected String or Symbol.")

function _extract_group_payload(var)
    if var isa Pair
        lhs, rhs = var
        return _to_param_string(lhs), rhs
    elseif var isa Tuple && length(var) == 2
        return _to_param_string(var[1]), var[2]
    else
        return _to_param_string(var), nothing
    end
end

function _build_qensemble_function(name::String, brace_elements::Vector{String}, f::Function)
    isempty(brace_elements) && error("Parameter \"$name\" requires a parentheses list specifying argument order when providing an ensemble function, e.g. \"$name(t, alpha)\" => (t, alpha) -> ...")
    arg_symbols = Symbol.(brace_elements)
    return QEnsembleFunction(name, arg_symbols, f)
end

function _ensure_group_capacity!(vec::BitVector, idx::Int)
    if length(vec) < idx
        new_len = max(idx, max(length(vec) * 2, INITIAL_PARAMETER_GROUP_MASK_SIZE))
        resize!(vec, new_len)
    end
end

function _mark_parameter_group_acting!(subspace::SubSpace, group_index::Int)
    mask = subspace.parameter_group_acting
    _ensure_group_capacity!(mask, group_index)
    mask[group_index] = true
end

function _set_parameter_group_distribution!(subspace::SubSpace, group_index::Int, is_distribution::Bool)
    mask = subspace.parameter_group_distribution
    _ensure_group_capacity!(mask, group_index)
    mask[group_index] = is_distribution
end

function _coerce_payload(name::String, brace_elements::Vector{String}, payload, has_indexes::Bool)
    dist = nothing
    qfun = nothing
    scalar_fun = nothing
    if payload === nothing
        return dist, qfun, scalar_fun
    elseif payload isa QDistribution
        dist = payload
    elseif payload isa Function
        if has_indexes
            qfun = _build_qensemble_function(name, brace_elements, payload)
        else
            scalar_fun = payload
        end
    elseif payload isa QEnsembleFunction
        qfun = payload
    else
        error("Unsupported payload type $(typeof(payload)) for parameter \"$name\".")
    end
    return dist, qfun, scalar_fun
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

function ParameterInfo(parameters::Vector{Parameter}, outer_labels_symbols::Vector{Symbol},
                       group_defs::Vector{ParameterGroupDefinition},
                       ensemble_group_functions::Vector{Union{Nothing,QEnsembleFunction}},
                       scalar_group_functions::Vector{Union{Nothing,Function}},
                       param_of_indexes::BitVector,
                       ss_ensemble_indexes_by_group::Vector{Vector{Int}},
                       ss_ensemble_present_by_group::Vector{BitVector},
                       subspace_index_maps::Vector{Array{SparsePermutation,2}},
                       t_index_transform::Array{SparsePermutation,2},
                       subspace_info::SubSpaceInfo)
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
            curr_bools = [falses( n) for n in ensemble_sizes]
            for curr_ind in param.param_indexes 
                outer = curr_ind.outer 
                inner = curr_ind.inner 
                outer_ind = subspace_info.ensemble_index_by_outer_index[outer]
                curr_bools[outer_ind][inner] = true  
            end
            push!(where_acting_by_parameter, curr_bools)
        end
    end
    params_acting_by_index = [ [ falses(length(parameters)) for _ in 1:n ] for n in ensemble_sizes ]
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

    group_count = length(group_defs)
    group_of_t = BitVector(gd.of_t for gd in group_defs)
    group_is_t = falses(group_count)
    group_name_to_index = Dict{Symbol,Int}(Symbol(def.name) => idx for (idx, def) in enumerate(group_defs))
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
        idx_names = group_defs[g].indexes
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
        qfun = ensemble_group_functions[g]
        scalar_fun = scalar_group_functions[g]
        args = group_defs[g].function_args
        qfun === nothing && scalar_fun === nothing && continue

        for idx in params_by_group[g]
            param = parameters[idx]
            refs = Vector{Int}(undef, length(args))
            main_index_map = Dict{String,Int}()
            required_index_tokens = Set(group_defs[g].indexes)
            used_index_tokens = Set{String}()
            for (name, sub_idx) in zip(group_defs[g].indexes, param.param_indexes)
                main_index_map[name] = sub_idx.inner
            end
            for (arg_pos, arg_str) in enumerate(args)
                if arg_str == "t"
                    t_key = param_coords[idx][1] - 1
                    refs[arg_pos] = get(time_param_lookup, t_key) do
                        error("No time parameter found for t$(t_key) when evaluating ensemble function for $(param.param_name).")
                    end
                    continue
                end
                arg_name, arg_tokens = underscore_separate(arg_str)
                target_group_idx = get(group_name_to_index, Symbol(arg_name)) do
                    error("Unknown parameter group $arg_name referenced in ensemble function for $(param.param_name).")
                end
                target_def = group_defs[target_group_idx]
                if qfun !== nothing
                    ensemble_group_functions[target_group_idx] === nothing || error("Ensemble function for $(param.param_name) cannot depend on function-defined group $(target_def.name).")
                end
                if isempty(target_def.indexes)
                    !isempty(arg_tokens) && error("Argument $arg_str should not specify indexes for scalar parameter group $(target_def.name).")
                else
                    length(arg_tokens) == length(target_def.indexes) || error("Argument $arg_str must specify indexes $(target_def.indexes) for parameter group $(target_def.name).")
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
            end
            if !isempty(required_index_tokens) && !(required_index_tokens ⊆ used_index_tokens)
                missing = setdiff(required_index_tokens, used_index_tokens)
                signature = group_defs[g].name
                if !isempty(group_defs[g].indexes)
                    signature *= "_{" * join(group_defs[g].indexes, ",") * "}"
                end
                if group_defs[g].of_t
                    signature *= "(t)"
                end
                message_missing = join(collect(missing), ", ")
                error("Ensemble function for $(signature) does not reference index(es) $(message_missing). Each underscore index in $(signature) must appear in the function arguments.")
            end
            function_param_refs[idx] = refs
        end
    end

    group_name_to_index = Dict{Symbol,Int}(Symbol(outer_labels_symbols[g]) => g for g in 1:group_count)
    param_name_to_indices = Dict{Symbol,Vector{Int}}()
    for idx in 1:length(parameters)
        coords = param_coords[idx]
        group_idx = param_group_by_index[idx]
        base_symbol = String(outer_labels_symbols[group_idx])
        sym_str = Symbol(param_strs[idx])
        sym_name = Symbol(param_names[idx])
        _register_param_key!(param_name_to_indices, sym_str, idx)
        _register_param_key!(param_name_to_indices, sym_name, idx)
        placeholder = _normalize_param_placeholder(param_names[idx])
        placeholder !== nothing && _register_param_key!(param_name_to_indices, placeholder, idx)
        numeric_key = _numeric_param_key(base_symbol, coords, param_of_t[idx])
        numeric_key !== nothing && _register_param_key!(param_name_to_indices, numeric_key, idx)
        if param_is_t[idx]
            t_idx = coords[1] - 1
            _register_param_key!(param_name_to_indices, Symbol("t$(t_idx)"), idx)
        end
    end

    time_slot_to_param = Dict(time_param_lookup)

    param_dicts = ParameterDicts(group_name_to_index, param_name_to_indices, time_slot_to_param)

    return CFunctions.ParameterInfo(outer_labels_symbols, inner_labels_symbols_flat, outer_labels, outer_labels_str, outer_labels_latex, param_names,
        param_strs, param_latex, param_of_indexes, param_group_by_index,
        t_index_by_index, ss_ensemble_indexes_by_group, ss_ensemble_present_by_group, indexed_parameter_indexes,
        where_acting_by_parameter, params_acting_by_index, param_index_tuples, subspace_index_maps, t_index_transform, indexes_by_t_index,
        indexes_of_t, ensemble_sizes, param_of_t, param_is_t, function_param_refs,
        group_time_counts, group_index_sizes, param_coords, params_by_group, group_of_t, group_is_t,
        subspace_info, param_indexes, param_dicts)
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
                                         max_t_ind::Int)::Tuple{Vector{Parameter}, ParameterInfo, ParameterValues}
    # --- start from a local copy and auto-add t if not present ---
    var_param = copy(vd.var_param)
    if all(group.name != "t" for group in var_param) && !(:t in used_symbols)
        push!(var_param, ParameterGroupDefinition("t", true, String[], nothing, nothing, nothing, String[]))
    end

    ensemble_group_distributions = Vector{Union{Nothing,QDistribution}}(undef, length(var_param))
    ensemble_group_functions = Vector{Union{Nothing,QEnsembleFunction}}(undef, length(var_param))
    scalar_group_functions = Vector{Union{Nothing,Function}}(undef, length(var_param))
    fill!(ensemble_group_distributions, nothing)
    fill!(ensemble_group_functions, nothing)
    fill!(scalar_group_functions, nothing)
    group_defs = var_param
    group_name_to_index = Dict{Symbol,Int}(Symbol(def.name) => i for (i, def) in enumerate(group_defs))
    outer_labels_symbols::Vector{Symbol} = Symbol[]
    parameters::Vector{Parameter} = Parameter[]
    ensemble_index_maps::Vector{Vector{Array{Int}}} = Vector{Vector{Array{Int}}}()
    t_index_maps::Vector{Vector{Int}} = Vector{Vector{Int}}()
    ss_ensemble_indexes_by_group::Vector{Vector{Int}} = Vector{Vector{Int}}()
    ss_ensemble_present_by_group::Vector{BitVector} = Vector{BitVector}()
    param_of_indexes::BitVector = Bool[]

    # ---- build variables + index maps per group ----
    for (group_index, group_def) in enumerate(var_param)
        param_name = group_def.name
        of_t = group_def.of_t
        index_strs = group_def.indexes
        dist = group_def.distribution
        qfun = group_def.ensemble_function
        scalar_fun = group_def.scalar_function
        ensemble_group_distributions[group_index] = dist
        ensemble_group_functions[group_index] = qfun
        scalar_group_functions[group_index] = scalar_fun
        var_name_sym::Symbol = Symbol(param_name)
        push!(outer_labels_symbols, var_name_sym)
        if var_name_sym in used_symbols
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

        if belongs_to_ensemble
            if dist === nothing && qfun === nothing
                error("Parameter \"$param_name\" references ensemble indexes and must be provided with a QDistribution or QEnsembleFunction.")
            end
            if length(unique_outers) > 1
                qfun !== nothing || error("Parameter \"$param_name\" spans multiple ensemble subspaces and therefore requires a QEnsembleFunction, e.g. \"$param_name(t, ...)\" => (args -> ...).")
                dist === nothing || error("Parameter \"$param_name\" spans multiple ensemble subspaces; QDistribution is not supported in this case.")
            end
        end

        if !isempty(unique_outers)
            param_sym = Symbol(param_name)
            for outer_ind in unique_outers
                ensemble_cfg = subspaces[outer_ind].ensemble
                if ensemble_cfg !== nothing && !(param_sym in ensemble_cfg.parameter_groups)
                    push!(ensemble_cfg.parameter_groups, param_sym)
                end
                _mark_parameter_group_acting!(subspaces[outer_ind], group_index)
                _set_parameter_group_distribution!(subspaces[outer_ind], group_index, dist !== nothing)
            end
        end

        blocks, block_lengths = find_blocks(outer_subsystem_inds)
        inner_label_symbols = subspace_info.inner_labels_symbols[blocks]
        ensemble_lengths = [length(inner_label_symbols[i]) for i in 1:length(inner_label_symbols) for _ in 1:block_lengths[i]]
        block_combinations = [collect(Combinatorics.with_replacement_combinations(1:ensemble_len, block_len))
                              for (ensemble_len, block_len) in zip(ensemble_lengths, block_lengths)]

        push!(ss_ensemble_indexes_by_group, unique(outer_subsystem_inds))
        ensemble_bool_vec::BitVector = in.(subspace_info.where_ensembles, Ref(unique(outer_subsystem_inds)))
        push!(ss_ensemble_present_by_group, ensemble_bool_vec)

        if !isempty(block_combinations)
            # -------- indexed group: real ensemble map, placeholder t-map --------
            index_map_vec::Vector{Array{Int}} = [zeros(Int, ensemble_lengths...) for _ in t_vals]

            for comb_comb in Iterators.product(block_combinations...)
                inner_subspace_inds = vcat(comb_comb...)
                symbol_comb = vcat([inner_labels[comb] for (inner_labels, comb) in zip(inner_label_symbols, comb_comb)]...)
                str_comb = string.(symbol_comb)
                var_name_str, var_name_latex = symbol2formatted(param_name, str_comb)
                curr_var_name = unformatted_var_name(param_name, str_comb)
                param_indexes::Vector{SubSpaceIndex} = [SubSpaceIndex(outer, inner, subspace_info)
                                                    for (outer, inner) in zip(outer_subsystem_inds, inner_subspace_inds)]
                for t_ind in t_vals
                    t_suff        = of_t ? "(" * t_suffix(t_ind) * ")" : ""
                    t_suff_latex  = of_t ? "(" * t_suffix(t_ind, do_latex=true) * ")" : ""
                    push!(parameters, Parameter(Symbol(param_name), curr_var_name*t_suff, var_name_str*t_suff,
                                                curr_var_name, var_name_str, var_name_latex*t_suff_latex, symbol_comb,
                                                of_t, false, t_ind, group_index, true, param_indexes))
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

    # ---- build subspace_index_maps (always Matrix{Vector{Int}}, 0×0 for non-ensembles) ----
    subspace_index_maps = build_subspace_index_maps(parameters, ensemble_index_maps, subspace_info)
    #println(ensemble_index_maps)
    #println(t_index_maps)
    # ---- finalize ParameterInfo ----
    var_info = ParameterInfo(parameters, outer_labels_symbols, group_defs,
                             ensemble_group_functions, scalar_group_functions,
                             param_of_indexes, ss_ensemble_indexes_by_group, ss_ensemble_present_by_group,
                             subspace_index_maps, t_index_transform, subspace_info)

    param_values = ParameterValues(var_info;
        ensemble_group_distributions=ensemble_group_distributions,
        ensemble_group_functions=ensemble_group_functions,
        group_functions=scalar_group_functions)

    return parameters, var_info, param_values
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
