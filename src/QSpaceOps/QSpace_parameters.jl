using Combinatorics

""" 
    Parameter(var_name::String, var_of_t::Bool, var_of_ensemble::Bool, var_ensemble_index::Int=0; var_val::Union{Nothing,Number,Vector{Number},Function}=nothing, var_suffix::String="")

Parameter is a struct that represents a variable in the state space, and information of how to access and print it.
"""
mutable struct Parameter
    var_symbol::Symbol
    var_name::String
    var_str::String
    var_latex::String
    var_val::Union{Nothing,Number,Vector{Number},Function}
    var_of_t::Bool
    t_index::Bool 
    group_index::Int 
    indexed_var::Bool
    var_indexes::Vector{SubSpaceIndex}
end

struct ParameterDefinitions
    var_param::Vector{Tuple{String, Bool, Vector{String}}}
    max_t_ind::Int
    function ParameterDefinitions(vars...; max_t_ind::Int=0)
        var_param::Vector{Tuple{String, Bool, Vector{String}}} = []
        for var in vars
            pre, indexes = underscore_separate(var)
            name, brace_elements = brace_separate(pre) 
            of_t = false
            if "t" in brace_elements 
                of_t = true
            end
            if length(brace_elements) > 1
                error("Currently only supports functions of t. Please send us your suggestions for what else you would want supported.")
            end
            push!(var_param, (name, of_t, indexes))
        end
        return new(var_param, max_t_ind)
    end
end

struct ParameterInfo 
    outer_labels_symbols::Vector{Symbol}
    inner_labels_symbols_flat::Vector{Symbol}
    
    outer_labels::Vector{String}
    inner_labels_flat::Vector{String}

    expanded_is_indexed::Vector{Bool}   
    outer_group_by_index::Vector{Int}    # continue here
    t_index_by_index::Vector{Int}
    ss_ensemble_indexes_by_group::Vector{Vector{Int}}    # which ss ensembles are used for indexing in each group. 
    ss_ensemble_present_by_group::Vector{Vector{Bool}}   # which ss ensembles are present in each group.

    ensemble_index_maps::Vector{Vector{Array{Int}}}    # maps different t_index and ensemble index combination to expanded variable indexes
    t_index_maps::Vector{Vector{Int}}                   # maps different t_index variations of a var group to expanded variable indexes
    function ParameterInfo(parameters::Vector{Parameter}, outer_labels_symbols::Vector{Symbol}, expanded_is_indexed::Vector{Bool}, ss_ensemble_indexes_by_group::Vector{Vector{Int}}, ensemble_index_maps::Vector{Vector{Array{Int}}}, t_index_maps::Vector{Vector{Int}})
        inner_labels_symbols_flat::Vector{Symbol} = Symbol[parameter.var_symbol]
        inner_labels_flat::Vector{String} = String[parameter.var_name]
        outer_labels::Vector{String} = String.(outer_labels_symbols)
        outer_group_by_index::Vector{Int} = zeros(Int, length(parameters) ) 
        t_index_by_index::Vector{Int} = zeros(Int, length(parameters) ) 
        for (i, param) in enumerate(parameters)
            outer_group_by_index[i] = param.group_index
            t_index_by_index[i] = param.t_index
        end
        return new(outer_labels_symbols, inner_labels_symbols_flat, outer_labels, inner_labels_flat, expanded_is_indexed, outer_group_by_index, t_index_by_index, ss_ensemble_indexes_by_group, ensemble_index_maps, t_index_maps)
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
function ParameterDefinitions2Parameters(vd::ParameterDefinitions, subspace_info::SubSpaceInfo, used_symbols::Vector{Symbol})::Tuple{Vector{Parameters}, ParameterInfo}
    var_param = vd.var_param    # elements defining the parameter group 
    max_t_ind = vd.max_t_ind
    
    outer_labels_symbols::Vector{Symbol}
    parameters::Vector{Parameter} = Parameter[]
    ensemble_index_maps::Vector{Vector{Array{Int}}} = Vector{Vector{Array{Int}}}()
    t_index_maps::Vector{Vector{Int}} = Vector{Vector{Int}}()
    ss_ensemble_indexes_by_group::Vector{Vector{Int}} = Vector{Vector{Int}}()
    ss_ensemble_present_by_group::Vector{Vector{Bool}} = Vector{Vector{Bool}}()
    expanded_is_indexed::Vector{Bool} = Vector{Bool}()
    for (group_index, (var_name, of_t, index_strs)) in enumerate(var_param) # <======= Create Variable Groups! 
        var_name_sym::Symbol = Symbol(var_name)
        push!(outer_labels_symbols, var_name_sym)   
        if var_name_sym in used_symbols            
            error("Variable name $var_name_sym is already used in the system! Please choose a distinct name.")
        end 
        if length(index_strs) > 0 
            var_of_ensemble = true
        end 
        t_vals = of_t ? collect(0:max_t_ind) : [0]
    
        # find ensemble indexes of the index_strs 
        index_str_syms = [Symbol(s) for s in index_strs]
        outer_subsystem_inds = -ones(Int, length(index_str_syms))
        for (i, index_sym) in enumerate(index_str_syms)
            subsystem_ind = nothing
            for outer_ind in 1:length(subspace_info.outer_labels_symbols)
                if index_sym in subspace_info.inner_label_symbols[outer_ind]
                    subsystem_ind = outer_ind    # inner_ind doesn't matter here! we iterate over them. 
                end
            end
            #subsystem_ind = findfirst(isequal(index_sym, subspace_info.inner_labels_symbols)) # <== Old way, using outer indexes only. but I think people prefer defining them via \omega_{i,j} for example. 
            if subsystem_ind === nothing
                ensemble_indexes = subspace_info.inner_labels[subspace_info.where_ensembles]
                error("Index $index_sym is not affiliated with a subsystem. 
                        The defined ensemble subsystems have the indexes $ensemble_indexes. 
                        Please choose one of those. ") 
            end

            if !(subsystem_ind in subspace_info.where_ensembles) 
                error("Index $index_sym is not an ensemble index.")
            end
            outer_subsystem_inds[i] = subsystem_ind
        end 
        if !contiguous_blocks(outer_subsystem_inds)
            error("Ensemble indexes must be contiguous. I.e. indexes belonging to the same ensemble must be grouped together. ")
        end
        blocks, block_lengths = find_blocks(outer_subsystem_inds)
        inner_label_symbols = subspace_info.inner_labels_symbols[blocks] 
        ensemble_lengths = length.(inner_label_symbols)
        block_combinations = [collect(Combinatorics.combinations(1:ensemble_len, block_len)) for (ensemble_len, block_len) in zip(ensemble_lengths, block_lengths)]
        push!(ss_ensemble_indexes_by_group, unique(outer_subsystem_inds))
        ensemble_bool_vec::Vector{Bool} = zeros(Bool, length(subspace_info.where_ensembles))
        ensemble_bool_vec[unique(outer_subsystem_inds)] .= true  # mark the ensembles that are
        push!(ss_ensemble_present_by_group, ensemble_bool_vec)

        # now construst all combinations of inner_label_symbols, by picking one from each bucket. 
        # post select only combinations in which each symbol appears at most once. 
        if length(block_combinations) > 0
            # create index maps 
            index_map_vec::Vector{Array{Int}} = [zeros(Int, ensemble_lengths...) for _ in t_vals]   # outer for t values, inner for ensemble indexes
            # create variables of different combinations:
            for comb_comb in Iterators.product(block_combinations...)
                inner_subspace_inds = vcat(comb_comb...)
                symbol_comb = vcat([inner_labels[comb] for (inner_labels, comb) in zip(inner_label_symbols, comb_comb)]...)
                str_comb = string.(symbol_comb)
                var_name_str, var_name_latex = symbol2formatted(var_name, str_comb)  
                curr_var_name = unformatted_var_name(var_name, str_comb) 
                var_indexes::Vector{SubSpaceIndex} = [SubSpaceIndex(outer, inner, subspace_info) for (outer, inner) in zip(outer_subsystem_inds, inner_subspace_inds)]
                for t_ind in t_vals
                    if of_t 
                        t_suff, t_suff_latex = "("* t_suffix(t_ind) *")", "("* t_suffix(t_ind, do_latex=true) *")" 
                    else 
                        t_suff, t_suff_latex = "", ""  
                    end
                    push!(parameters, Parameter(var_name_sym, curr_var_name*t_suff, var_name_str*t_suff, var_name_latex*t_suff_latex, 
                                                nothing, of_t, t_ind, true, var_indexes))
                    index_map_vec[t_ind+1][inner_subspace_inds...] = length(parameters)
                    push!(expanded_is_indexed, true)
                end
            end
            push!(ensemble_index_maps, index_map_vec)    
            
        else
            var_name_str, var_name_latex = symbol2formatted(var_name) 
            t_index_map_vec::Vector{Int} = Vector{Int}(undef, length(t_vals)
            for t_ind in t_vals
                if of_t 
                    t_suff, t_suff_latex = "("* t_suffix(t_ind) *")", "("* t_suffix(t_ind, do_latex=true) *")" 
                else 
                    t_suff, t_suff_latex = "", ""  
                end
                push!(parameters, Parameter(var_name_sym, var_name*t_suff, var_name_str*t_suff, var_name_latex*t_suff_latex, 
                                            nothing, of_t, t_ind, group_index, false, Vector{SubSpaceIndex}()))
                t_index_map_vec[t_ind+1] = length(parameters)
                push!(expanded_is_indexed, false)
            end
            push!(t_index_maps, t_index_map_vec)
        end
        push!(used_symbols, var_name_sym)
    end 
    var_info = ParameterInfo((parameters, outer_labels_symbols, ss_ensemble_indexes_by_group, ss_ensemble_present_by_group, ensemble_index_maps, t_index_maps)
    return parameters, var_info
end 
