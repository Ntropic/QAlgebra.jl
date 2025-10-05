using Combinatorics
using SparseArrays
using ..CFunctions: ParameterInfo, ParameterIndexes
using ..SparsePermutationTools: SparsePermutation, denseperm

""" 
    Parameter(param_name::String, param_of_t::Bool, var_of_ensemble::Bool, var_ensemble_index::Int=0; param_values::Union{Nothing,Number,Vector{Number},Function}=nothing, var_suffix::String="")

Parameter is a struct that represents a variable in the state space, and information of how to access and print it.
"""
mutable struct Parameter
    param_symbol::Symbol
    param_name::String
    param_str::String
    param_name_no_t::String
    param_str_no_t::String
    param_latex::String
    index_comb_symbol::Vector{Symbol}
    param_values::Union{Nothing,Number,Vector{Number},Function}
    param_of_t::Bool
    is_t::Bool
    t_index::Int 
    group_index::Int 
    indexed_param::Bool
    param_indexes::Vector{SubSpaceIndex}
end

"""
    ParameterDefinitions(params...)

A helper to construct Parameter Info and Parameter Vector of all Parameters. 
Supports arbitrarily many String or Symbol inputs which define parameters. Use "(t)" at the end of the name to  specify that it is time dependent. 
"""
struct ParameterDefinitions
    var_param::Vector{Tuple{String, Bool, Vector{String}}}
    function ParameterDefinitions(params...)
        var_param::Vector{Tuple{String, Bool, Vector{String}}} = []
        for var in params
            pre, brace_elements = brace_separate(var) 
            name, indexes = underscore_separate(pre)
            of_t = false
            if "t" in brace_elements 
                of_t = true
            end
            if length(brace_elements) > 1
                error("Currently only supports functions of t. Please send us your suggestions for what else you would want supported.")
            end
            push!(var_param, (name, of_t, indexes))
        end
        return new(var_param)
    end
end
function Base.show(io::IO, param_def::ParameterDefinitions)
    var_str_vec = []
    for (p, do_t, elem) in param_def.var_param
        param_str = symbol2formatted(p)[1]
        if do_t
            param_str *= "(t)"
        end
        param_str *= str2sub(join(elem, ","))
        push!(var_str_vec, param_str) 
    end 
    println(io, "ParameterDefinitions: [" * join(var_str_vec, ", ") * "]")
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

function ParameterInfo(parameters::Vector{Parameter}, outer_labels_symbols::Vector{Symbol}, param_of_indexes::BitVector,
                       ss_ensemble_indexes_by_group::Vector{Vector{Int}}, ss_ensemble_present_by_group::Vector{BitVector},
                       subspace_index_maps::Vector{Array{SparsePermutation,2}}, t_index_transform::Array{SparsePermutation,2}, subspace_info::SubSpaceInfo)
    inner_labels_symbols_flat = [param.param_symbol for param in parameters]
    param_names  = [param.param_name  for param in parameters]
    param_strs   = [param.param_str   for param in parameters]
    param_latex  = [param.param_latex for param in parameters]

    param_of_t::BitVector = [param.param_of_t for param in parameters]
    param_is_t::BitVector = [param.is_t       for param in parameters]
    param_values = [param.param_values for param in parameters]

    outer_labels = String.(outer_labels_symbols)
    outer_group_by_index = zeros(Int, length(parameters))
    t_index_by_index = zeros(Int, length(parameters))
    indexes_of_t = Int[]

    for (i, param) in enumerate(parameters)
        outer_group_by_index[i] = param.group_index
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
    acting_parameters_by_index = [ [ falses(length(parameters)) for _ in 1:n ] for n in ensemble_sizes ]
    for (param_idx, storage_idx) in pairs(indexed_parameter_indexes)
        storage_idx == 0 && continue
        param_acts = where_acting_by_parameter[storage_idx]
        for (ensemble_idx, bits) in enumerate(param_acts)
            for inner_idx in findall(bits)
                acting_parameters_by_index[ensemble_idx][inner_idx][param_idx] = true
            end
        end
    end

    param_indexes = ParameterIndexes(subspace_info, indexed_parameter_indexes, where_acting_by_parameter, indexes_by_t_index)

    parameter_index_tuples = Vector{Vector{Tuple{Int,Int}}}(undef, length(parameters))
    for (idx, param) in enumerate(parameters)
        if param.indexed_param
            tuples = Vector{Tuple{Int,Int}}(undef, length(param.param_indexes))
            @inbounds for (inner_pos, sub_idx) in enumerate(param.param_indexes)
                ensemble = subspace_info.ensemble_index_by_outer_index[sub_idx.outer]
                ensemble != 0 || error("Parameter index does not belong to an ensemble subspace.")
                tuples[inner_pos] = (ensemble, sub_idx.inner)
            end
            parameter_index_tuples[idx] = tuples
        else
            parameter_index_tuples[idx] = Tuple{Int,Int}[]
        end
    end

    return ParameterInfo(outer_labels_symbols, inner_labels_symbols_flat, outer_labels, param_names,
        param_strs, param_latex, param_of_indexes, outer_group_by_index,
        t_index_by_index, ss_ensemble_indexes_by_group, ss_ensemble_present_by_group, indexed_parameter_indexes,
        where_acting_by_parameter, acting_parameters_by_index, parameter_index_tuples, subspace_index_maps, t_index_transform, indexes_by_t_index,
        indexes_of_t, ensemble_sizes, param_of_t, param_is_t, param_values, subspace_info, param_indexes)
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
                                         max_t_ind::Int)::Tuple{Vector{Parameter}, ParameterInfo}
    # --- start from a local copy and auto-add t if not present ---
    var_param = copy(vd.var_param)
    if all(name != "t" for (name, _of_t, _idxs) in var_param) && !(:t in used_symbols)
        push!(var_param, ("t", true, String[]))
    end

    outer_labels_symbols::Vector{Symbol} = Symbol[]
    parameters::Vector{Parameter} = Parameter[]
    ensemble_index_maps::Vector{Vector{Array{Int}}} = Vector{Vector{Array{Int}}}()
    t_index_maps::Vector{Vector{Int}} = Vector{Vector{Int}}()
    ss_ensemble_indexes_by_group::Vector{Vector{Int}} = Vector{Vector{Int}}()
    ss_ensemble_present_by_group::Vector{BitVector} = Vector{BitVector}()
    param_of_indexes::BitVector = Bool[]

    # ---- build variables + index maps per group ----
    for (group_index, (param_name, of_t, index_strs)) in enumerate(var_param)
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

        if !isempty(index_str_syms)
            unique_outers = unique(filter(x->x>0, outer_subsystem_inds))
            if !isempty(unique_outers)
                param_sym = Symbol(param_name)
                for outer_ind in unique_outers
                    ensemble_cfg = subspaces[outer_ind].ensemble
                    if ensemble_cfg !== nothing && !(param_sym in ensemble_cfg.parameter_groups)
                        push!(ensemble_cfg.parameter_groups, param_sym)
                    end
                end
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
                                                nothing, of_t, false, t_ind, group_index, true, param_indexes))
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
                                            nothing, of_t, is_t, t_ind, group_index, false, SubSpaceIndex[]))
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
    var_info = ParameterInfo(parameters, outer_labels_symbols, param_of_indexes,
                             ss_ensemble_indexes_by_group, ss_ensemble_present_by_group,
                             subspace_index_maps, t_index_transform, subspace_info)

    return parameters, var_info
end

#Return parameter mapping vector for switching ensemble inner index within a subspace.
function map_by_subspace(i_to::SubSpaceIndex, i_from::SubSpaceIndex, pinfo::ParameterInfo)::Vector{Int}
    @assert i_to.outer == i_from.outer "Subspace mapping requires the same outer subspace."
    M = pinfo.subspace_index_maps[i_to.outer]
    if size(M,1) == 0  # non-ensemble subspace -> identity map
        return collect(1:length(pinfo.outer_group_by_index))
    else
        return denseperm(M[i_to.inner, i_from.inner])
    end
end

# Return parameter mapping vector for switching from t_index2 to t_index1.
map_by_tindex(t_index1::Int, t_index2::Int, pinfo::ParameterInfo) = denseperm(pinfo.t_index_transform[t_index1+1, t_index2+1])
# from t_index2 to t_index1 
