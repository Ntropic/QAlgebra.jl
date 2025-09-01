using Combinatorics

""" 
    Parameter(var_name::String, var_of_t::Bool, var_of_ensemble::Bool, var_ensemble_index::Int=0; var_val::Union{Nothing,Number,Vector{Number},Function}=nothing, var_suffix::String="")

Parameter is a struct that represents a variable in the state space, and information of how to access and print it.
"""
mutable struct Parameter
    var_symbol::Symbol
    var_name::String
    var_str::String
    var_name_no_t::String
    var_str_no_t::String
    var_latex::String
    index_comb_symbol::Vector{Symbol}
    var_val::Union{Nothing,Number,Vector{Number},Function}
    var_of_t::Bool
    is_t::Bool
    t_index::Int 
    group_index::Int 
    indexed_var::Bool
    var_indexes::Vector{SubSpaceIndex}
end

"""
    ParameterDefinitions(vars...)

A helper to construct Parameter Info and Parameter Vector of all Parameters. 
Supports arbitrarily many String or Symbol inputs which define parameters. Use "(t)" at the end of the name to  specify that it is time dependent. 
"""
struct ParameterDefinitions
    var_param::Vector{Tuple{String, Bool, Vector{String}}}
    function ParameterDefinitions(vars...)
        var_param::Vector{Tuple{String, Bool, Vector{String}}} = []
        for var in vars
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
        var_str = symbol2formatted(p)[1]
        if do_t
            var_str *= "(t)"
        end
        var_str *= str2sub(join(elem, ","))
        push!(var_str_vec, var_str) 
    end 
    println(io, "ParameterDefinitions: [" * join(var_str_vec, ", ") * "]")
end

struct ParameterInfo 
    outer_labels_symbols::Vector{Symbol}
    inner_labels_symbols_flat::Vector{Symbol}
    
    outer_labels::Vector{String}
    inner_labels_flat::Vector{String}

    expanded_is_indexed::Vector{Bool}   
    outer_group_by_index::Vector{Int}   
    t_index_by_index::Vector{Int}       # -1 for parameters that aren't of t. 
    ss_ensemble_indexes_by_group::Vector{Vector{Int}}    # which ss ensembles are used for indexing in each group. 
    ss_ensemble_present_by_group::Vector{Vector{Bool}}   # which ss ensembles are present in each group.

    indexed_parameter_indexes::Vector{Int}                  # which parameters have indexes?
    where_acting_by_parameter::Vector{Vector{Vector{Bool}}}  # for each variable, where are they acting. 

    # Maps indexes for index transformation, once for switching subsystem indexes and once for time indexes
    subspace_index_maps::Vector{Array{Vector{Int},2}}
    t_index_transform::Array{Vector{Int},2}
    indexes_by_t_index::Vector{Vector{Int}}   # for each t_index which indexes have it? 
    indexes_of_t::Vector{Int}

    function ParameterInfo(parameters::Vector{Parameter},
                           outer_labels_symbols::Vector{Symbol},
                           expanded_is_indexed::Vector{Bool},
                           ss_ensemble_indexes_by_group::Vector{Vector{Int}},
                           ss_ensemble_present_by_group::Vector{Vector{Bool}},
                           subspace_index_maps::Vector{Array{Vector{Int},2}},
                           t_index_transform::Array{Vector{Int},2}, 
                           subspace_info::SubSpaceInfo)
        inner_labels_symbols_flat::Vector{Symbol} = Symbol[param.var_symbol for param in parameters]
        inner_labels_flat::Vector{String} = String[param.var_name for param in parameters]
        outer_labels::Vector{String} = String.(outer_labels_symbols)
        outer_group_by_index::Vector{Int} = zeros(Int, length(parameters))
        t_index_by_index::Vector{Int} = zeros(Int, length(parameters))
        indexes_of_t::Vector{Int} = []
        for (i, param) in enumerate(parameters)
            outer_group_by_index[i] = param.group_index
            t_index_by_index[i] = param.t_index - !param.var_of_t
            if param.var_of_t 
                push!(indexes_of_t, i)
            end
        end
        indexes_by_t_index::Vector{Vector{Int}} = [findall(==(t_ind), t_index_by_index) for t_ind in 0:maximum(t_index_by_index)]

        indexed_parameter_indexes::Vector{Int} = []
        where_acting_by_parameter::Vector{Vector{Vector{Bool}}} = []
        ensemble_sizes = subspace_info.how_many_by_ensemble
        
        for (i, param) in enumerate(parameters)
            if param.indexed_var 
                push!(indexed_parameter_indexes, i)
                curr_bools::Vector{Vector{Bool}} = [zeros(Bool, n) for n in ensemble_sizes]
                for curr_ind in param.var_indexes 
                    outer = curr_ind.outer 
                    inner = curr_ind.inner 
                    outer_ind = subspace_info.ensemble_index_by_outer_index[outer]
                    curr_bools[outer_ind][inner] = true  
                end
                push!(where_acting_by_parameter, curr_bools) # param.var_indexes)
            end
        end
        return new(outer_labels_symbols, inner_labels_symbols_flat, outer_labels, inner_labels_flat,
                   expanded_is_indexed, outer_group_by_index, t_index_by_index,
                   ss_ensemble_indexes_by_group, ss_ensemble_present_by_group, 
                   indexed_parameter_indexes, where_acting_by_parameter,
                   subspace_index_maps, t_index_transform, indexes_by_t_index, indexes_of_t)
    end
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

function unformatted_var_name(var_name, indexes::Vector{String})
    index_str = join(indexes, ",")
    if length(indexes) > 0 
        return var_name * "_{" * index_str * "}"
    else 
        return var_name * "_" * index_str 
    end
end

function build_subspace_index_maps(parameters::Vector{Parameter}, ensemble_index_maps, s_info::SubSpaceInfo)
    N = length(parameters)
    nsub = length(s_info.outer_labels_symbols)
    result = Vector{Array{Vector{Int},2}}(undef, nsub)

    ensemble_set = Set(s_info.where_ensembles)

    for s in 1:nsub
        if s in ensemble_set
            m = length(s_info.inner_labels_symbols[s])  # ensemble size for subspace s
            pairmat = Array{Vector{Int}}(undef, m, m)   # (inner_to, inner_from)
            for inner_to in 1:m, inner_from in 1:m
                mapvec = Vector{Int}(undef, N)
                @inbounds for (idx, param) in enumerate(parameters)
                    if !param.indexed_var
                        mapvec[idx] = idx
                        continue
                    end
                    # remap only if this param actually uses subspace s at inner_from
                    has_s = false
                    inners = Vector{Int}(undef, length(param.var_indexes))
                    for a in eachindex(param.var_indexes)
                        ix = param.var_indexes[a]
                        v = ix.inner
                        if ix.outer == s && ix.inner == inner_from
                            v = inner_to
                            has_s = true
                        end
                        inners[a] = v
                    end
                    if has_s
                        # dst = ensemble_index_maps[param.group_index][param.t_index+1][inners...]  ### unsorted variant 
                        inners_sorted = sort(inners)
                        dst = ensemble_index_maps[param.group_index][param.t_index+1][inners_sorted...]
                        mapvec[idx] = dst
                    else
                        mapvec[idx] = idx
                    end
                end
                
                pairmat[inner_to, inner_from] = mapvec
            end
            result[s] = pairmat
        else
            # non-ensemble subspace -> empty (0×0) matrix
            result[s] = Array{Vector{Int}}(undef, 0, 0)
        end
    end
    for s in 1:nsub
        for mat in result[s]
            if any(==(0), mat)
                error("build_subspace_index_maps: unmapped index (0) detected in subspace $s")
            end
        end
    end
    return result
end

function ParameterDefinitions2Parameters(vd::ParameterDefinitions, subspace_info::SubSpaceInfo, used_symbols::Set{Symbol}, max_t_ind::Int)::Tuple{Vector{Parameter}, ParameterInfo}
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
    ss_ensemble_present_by_group::Vector{Vector{Bool}} = Vector{Vector{Bool}}()
    expanded_is_indexed::Vector{Bool} = Bool[]

    # ---- build variables + index maps per group ----
    for (group_index, (var_name, of_t, index_strs)) in enumerate(var_param)
        var_name_sym::Symbol = Symbol(var_name)
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
            for outer_ind in 1:length(subspace_info.outer_labels_symbols)
                if index_sym in subspace_info.inner_labels_symbols[outer_ind]
                    subsystem_ind = outer_ind
                end
            end
            if subsystem_ind === nothing
                ensemble_indexes = subspace_info.inner_labels[subspace_info.where_ensembles]
                error("Index $index_sym is not affiliated with a subsystem. The defined ensemble subsystems have the indexes $ensemble_indexes.")
            end
            if !(subsystem_ind in subspace_info.where_ensembles)
                error("Index $index_sym is not an ensemble index.")
            end
            outer_subsystem_inds[i] = subsystem_ind
        end

        if !contiguous_blocks(outer_subsystem_inds)
            error("Ensemble indexes must be contiguous. Indexes belonging to the same ensemble must be grouped.")
        end

        blocks, block_lengths = find_blocks(outer_subsystem_inds)
        inner_label_symbols = subspace_info.inner_labels_symbols[blocks]
        ensemble_lengths = [length(inner_label_symbols[i]) for i in 1:length(inner_label_symbols) for _ in 1:block_lengths[i]]
        block_combinations = [collect(Combinatorics.with_replacement_combinations(1:ensemble_len, block_len))
                              for (ensemble_len, block_len) in zip(ensemble_lengths, block_lengths)]

        push!(ss_ensemble_indexes_by_group, unique(outer_subsystem_inds))
        ensemble_bool_vec::Vector{Bool} = in.(subspace_info.where_ensembles, Ref(unique(outer_subsystem_inds)))
        push!(ss_ensemble_present_by_group, ensemble_bool_vec)

        if !isempty(block_combinations)
            # -------- indexed group: real ensemble map, placeholder t-map --------
            index_map_vec::Vector{Array{Int}} = [zeros(Int, ensemble_lengths...) for _ in t_vals]

            for comb_comb in Iterators.product(block_combinations...)
                inner_subspace_inds = vcat(comb_comb...)
                symbol_comb = vcat([inner_labels[comb] for (inner_labels, comb) in zip(inner_label_symbols, comb_comb)]...)
                str_comb = string.(symbol_comb)
                var_name_str, var_name_latex = symbol2formatted(var_name, str_comb)
                curr_var_name = unformatted_var_name(var_name, str_comb)
                var_indexes::Vector{SubSpaceIndex} = [SubSpaceIndex(outer, inner, subspace_info)
                                                    for (outer, inner) in zip(outer_subsystem_inds, inner_subspace_inds)]
                for t_ind in t_vals
                    t_suff        = of_t ? "(" * t_suffix(t_ind) * ")" : ""
                    t_suff_latex  = of_t ? "(" * t_suffix(t_ind, do_latex=true) * ")" : ""
                    push!(parameters, Parameter(Symbol(var_name), curr_var_name*t_suff, var_name_str*t_suff,
                                                curr_var_name, var_name_str, var_name_latex*t_suff_latex, symbol_comb,
                                                nothing, of_t, false, t_ind, group_index, true, var_indexes))
                    index_map_vec[t_ind+1][inner_subspace_inds...] = length(parameters)
                    push!(expanded_is_indexed, true)
                end
            end
            push!(ensemble_index_maps, index_map_vec)         # real
            push!(t_index_maps, Vector{Int}())                # placeholder to keep g-alignment

        else
            # -------- non-indexed group: real t-map, placeholder ensemble map --------
            var_name_str, var_name_latex = symbol2formatted(var_name)
            t_index_map_vec::Vector{Int} = Vector{Int}(undef, length(t_vals))

            for t_ind in t_vals
                is_t = (var_name_sym==:t)
                t_suff        = (of_t && !is_t) ? "(" * t_suffix(t_ind) * ")" : ""
                t_suff_latex  = (of_t && !is_t) ? "(" * t_suffix(t_ind, do_latex=true) * ")" : ""
                
                t_suff        = (t_ind > 0 && is_t) ? str2sub(string(t_ind))  : ""
                t_suff_latex  = (t_ind > 0 && is_t) ? "_$t_ind" : ""
                push!(parameters, Parameter(Symbol(var_name), var_name*t_suff, var_name_str*t_suff,
                                            var_name, var_name_str, var_name_latex*t_suff_latex, Symbol[],
                                            nothing, of_t, is_t, t_ind, group_index, false, SubSpaceIndex[]))
                t_index_map_vec[t_ind+1] = length(parameters)
                push!(expanded_is_indexed, false)
            end

            push!(t_index_maps, t_index_map_vec)                              # real
            push!(ensemble_index_maps, [zeros(Int, 0) for _ in t_vals])       # placeholder to keep g-alignment
        end
        push!(used_symbols, Symbol(var_name))
    end

    # ---- build t_index_transform (t_to, t_from) -> mapping vector over full parameter space ----
    N = length(parameters)
    T = max_t_ind + 1
    t_index_transform = Array{Vector{Int}}(undef, T, T)
    for t_to in 0:max_t_ind
        for t_from in 0:max_t_ind
            mapvec = Vector{Int}(undef, N)
            @inbounds for (idx, param) in enumerate(parameters)
                if !param.var_of_t || param.t_index != t_from
                    mapvec[idx] = idx
                    continue
                end
                g = param.group_index
                if param.indexed_var

                    inners = [ix.inner for ix in param.var_indexes]
                    dst = ensemble_index_maps[g][t_to+1][inners...]
                    mapvec[idx] = dst
                else
                    mapvec[idx] = t_index_maps[g][t_to+1]
                end
            end
            t_index_transform[t_to+1, t_from+1] = mapvec
        end
    end

    # ---- build subspace_index_maps (always Matrix{Vector{Int}}, 0×0 for non-ensembles) ----
    subspace_index_maps = build_subspace_index_maps(parameters, ensemble_index_maps, subspace_info)
    #println(ensemble_index_maps)
    #println(t_index_maps)
    # ---- finalize ParameterInfo ----
    var_info = ParameterInfo(parameters, outer_labels_symbols, expanded_is_indexed,
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
        return M[i_to.inner, i_from.inner]
    end
end

# Return parameter mapping vector for switching from t_index2 to t_index1.
map_by_tindex(t_index1::Int, t_index2::Int, pinfo::ParameterInfo) = pinfo.t_index_transform[t_index1+1, t_index2+1]