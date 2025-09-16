""" 
    SubSpace(key::String, keys::Vector{String}, ss_outer_ind::Int, ss_inner_ind::Vector{Int}, op_set::OperatorSet, ensemble::Bool, fermion::Bool)

SubSpace defines a Subspace of a Hilbert space. It contains an operator set, aswell as additional information to reference and work with a subspace. 
Subspaces can be divided into of sub-subsystems (internally referred to as inner subsystems), multiple copies of the same subspace, so as to support ensemble descriptions.
"""
struct SubSpace
    key_symbol::Symbol
    keys_symbols::Vector{Symbol}
    key::String                     # Original input key
    keys::Vector{String}            # Allowed keys for this subspace 
    ss_outer_ind::Int            # Which Vector to use for ss_inner_ind  (this is for accessing the string elements)
    ss_inner_ind::Vector{Int}    # Indices to access operator values in the corresponding statespace main ind  (this is for accessing the string elements)
    is_ensemble_ss::Bool
    ensemble_size::Int
    num_operator_indexes::Int 
    num_sum_indexes::Int
    particle_type::String
    op_set::OperatorSet             # The operator set for this subspace.
end
# Define the custom show for SubSpace.
function Base.show(io::IO, statespace::SubSpace)
    # Print the subspace key and allowed keys.
    print(io, "SubSpace ", statespace.keys[1:statespace.num_operator_indexes], " , ∑", statespace.keys[statespace.num_operator_indexes+1:end], ": ")
    # Use the OperatorSet's show for the op_set field.
    show(io, statespace.op_set)
end

""" 
    SubSpaceDefinitions(; kwargs...)

SubSpaceDefinitions is a struct that processes the definition of the quantum subspaces, patth keyword arguments, the key being the Symbol used to identify the subspace, 
    and the argument either being an OperatorType or a Tuple specifying the number of ensemble indexes and summation indexes.
"""
struct SubSpaceDefinitions 
    subspaces::Vector{SubSpace}  # Vector of all sub
    used_symbols::Set{Symbol}
    I_op::Vector{Is}  # Vector of all neutral elements
    I_ensemble_op::Vector{Vector{Is}}  # Vector of all neutral elements within ensembe subsystems 

    function SubSpaceDefinitions(;kwargs...)
        subspaces = Vector{SubSpace}()
        used_symbols = Set{Symbol}()
        key_counter = 0
        core_keys::Vector{Symbol} = collect(keys(kwargs))
        for (outer_ind, (key_symbol, val)) in enumerate(kwargs) 
            is_ensemble_ss = false
            if isa(val, Tuple) 
                is_ensemble_ss = true
                num_operator_indexes, num_sum_indexes, op_set = val  # unpacking
                ensemble_size = num_operator_indexes + num_sum_indexes
            else
                ensemble_size, op_set = 1, val 
                num_operator_indexes = 1
                num_sum_indexes = 0
            end
            key = String(key_symbol) 
            if length(key) == 1
                key_char = key[1]
                keys = String[string(key_char+i) for i in 0:ensemble_size-1]
                keys_symbols = Symbol.(keys) 
                if any(x->x in used_symbols, keys_symbols) 
                    keys = String[key*string(i) for i in 1:ensemble_size]
                    keys_symbols = Symbol.(keys)
                end
            else 
                keys = String[key*string(i) for i in 1:ensemble_size]
                keys_symbols = Symbol.(keys)
            end
            if any(x->x in used_symbols, keys_symbols) 
                error("Symbol $key already used")
            end 
            curr_inds = key_counter .+ collect(1:ensemble_size)
            curr_subspace = SubSpace(key_symbol, keys_symbols, key, keys, outer_ind, curr_inds, is_ensemble_ss, 
                        ensemble_size, num_operator_indexes, num_sum_indexes, op_set.particle_type, op_set) 
            key_counter += ensemble_size
            push!(subspaces, curr_subspace)
            union!(used_symbols, keys_symbols)
        end
        I_op::Vector{Is} = [s.op_set.neutral_element for s in subspaces for _ in 1:s.ensemble_size]
        I_ensemble_op::Vector{Vector{Is}} = [I_op[s.ss_inner_ind] for s in subspaces]
        return new(subspaces, used_symbols, I_op, I_ensemble_op)
    end
end
function Base.show(io::IO, subspace_def::SubSpaceDefinitions)
    println(io, "SubSpaceDefinitions: ")
    # Then print each subspace on its own line.
    for ss in subspace_def.subspaces
        println(io, "   - ", string(ss))
    end
end

struct SubSpaceInfo
    outer_labels_symbols::Vector{Symbol}
    inner_labels_symbols::Vector{Vector{Symbol}}
    inner_labels_symbols_flat::Vector{Symbol}
    
    outer_labels::Vector{String}               # Letters identifying outer subspaces
    inner_labels::Vector{Vector{String}}       # Letters identifying the inder subspaces 
    param_names::Vector{String}         # Letters identifying the inder subspaces flattened

    subsystem_sizes::Vector{Int}                # How many inner subsystems in each other subsystem
    outer_ss_of_expanded::Vector{Int}          # Identify the outer subsystem for each expanded index
    inner_ss_of_expanded::Vector{Int}          # Identify the inner subsystem for each expanded index

    expanded_index_by_outer::Vector{Vector{Int}} # Which indexes belong to an outer index 
    where_ensembles::Vector{Int}               # which subspaces are ensembles? 
    ensemble_index_by_outer_index::Vector{Int}  # which index among ensemble indexes is an outer index (gives 0 if not an ensemble subspace)


    how_many_by_ensemble::Vector{Int}
    how_many_non_sum_by_ensemble::Vector{Int}
    how_many_sum_by_ensemble::Vector{Int}
    ensemble_indexes::Vector{Vector{Int}}      # The expanded indexes of the ensembles 
    summation_indexes::Vector{Vector{Int}}     # The expanded indexes of the ensembles, that are for summations  
    of_time::Bool
end
# Primary constructor from labels
function SubSpaceInfo(outer_labels_symbols::Vector{Symbol}, inner_labels_symbols::Vector{Vector{Symbol}}, are_ensemble_ss::Vector{Bool}, 
                      num_operator_indexes::Vector{Int}, num_summation_indexes::Vector{Int}, of_time::Bool=false)
    inner_labels_symbols_flat = vcat(inner_labels_symbols...)
    outer_labels = map(string, outer_labels_symbols)
    inner_labels = [string.(v) for v in inner_labels_symbols]
    param_names = map(string, inner_labels_symbols_flat)
    @assert length(outer_labels) == length(inner_labels)
    subsystem_sizes = Int[length(v) for v in inner_labels]
    @assert all(>=(1), subsystem_sizes) "Every subsystem must have at least one inner element."
    n_exp = sum(subsystem_sizes)

    outer_ss_of_expanded = Vector{Int}(undef, n_exp)
    inner_ss_of_expanded = Vector{Int}(undef, n_exp)
    e = 1
    for o in eachindex(subsystem_sizes)
        @inbounds for i in 1:subsystem_sizes[o]
            outer_ss_of_expanded[e] = o
            inner_ss_of_expanded[e] = i
            e += 1
        end
    end
    expanded_index_by_outer::Vector{Vector{Int}} = []
    counter = 0
    for (i, s) in enumerate(subsystem_sizes)
        push!(expanded_index_by_outer, [counter+j for j in 1:s])
        counter += s
    end

    where_ensembles::Vector{Int} = []
    ensemble_indexes::Vector{Vector{Int}} = []
    ensemble_index_by_outer_index::Vector{Int} = zeros(Int, length(subsystem_sizes))
    for o in eachindex(subsystem_sizes)
        if are_ensemble_ss[o] 
            push!(where_ensembles, o)
            push!(ensemble_indexes, expanded_index_by_outer[o])
            ensemble_index_by_outer_index[o] = length(where_ensembles)
        end
    end
    how_many_by_ensemble::Vector{Int} = [length(x) for x in ensemble_indexes]
    how_many_non_sum_by_ensemble::Vector{Int} = num_operator_indexes[where_ensembles]
    how_many_sum_by_ensemble::Vector{Int} =  num_summation_indexes[where_ensembles]
    summation_indexes::Vector{Vector{Int}} = [curr_ensemble_indexes[how_many_non+1:end] for (how_many_non, curr_ensemble_indexes) in zip(how_many_non_sum_by_ensemble, ensemble_indexes)]
    return SubSpaceInfo( outer_labels_symbols, inner_labels_symbols, inner_labels_symbols_flat, outer_labels, inner_labels, param_names, 
                         subsystem_sizes, outer_ss_of_expanded, inner_ss_of_expanded, expanded_index_by_outer,
                         where_ensembles, ensemble_index_by_outer_index, how_many_by_ensemble, how_many_non_sum_by_ensemble, how_many_sum_by_ensemble, 
                         ensemble_indexes, summation_indexes, of_time )
end

function Base.show(io::IO, info::SubSpaceInfo)
    print(io, "SubSpaceInfo:\n")
    print(io, "  outer_labels:            ", info.outer_labels, "\n")
    print(io, "  inner_labels:            ", info.inner_labels, "\n")
    print(io, "  subsystem_sizes:          ", info.subsystem_sizes, "\n")
    print(io, "  outer_ss_of_expanded:    ", info.outer_ss_of_expanded, "\n")
    print(io, "  inner_ss_of_expanded:    ", info.inner_ss_of_expanded, "\n")
    print(io, "  expanded_index_by_outer: ", info.expanded_index_by_outer, "\n")
end

# Convenience: build from your existing `SubSpace` vector
function SubSpaceInfo(subspaces::Vector{SubSpace})
    outers  = [s.key_symbol        for s in subspaces]
    inners  = [copy(s.keys_symbols) for s in subspaces]   # for non-ensemble subspaces this is length 1
    are_ensemble_ss = [s.is_ensemble_ss for s in subspaces] 
    num_operator_indexes = [s.num_operator_indexes for s in subspaces]
    num_summation_indexes = [s.num_sum_indexes for s in subspaces]
    return SubSpaceInfo(outers, inners, are_ensemble_ss, num_operator_indexes, num_summation_indexes)
end

@inline function outer_inner_2_expanded(info::SubSpaceInfo, outer::Int, inner::Int=1)
    return info.expanded_index_by_outer[outer][inner] 
end
@inline function expanded_2_outer_inner(info::SubSpaceInfo, expanded::Int)
    return (info.outer_ss_of_expanded[expanded], info.inner_ss_of_expanded[expanded])
end
@inline function label_2_outer_inner_expanded(info::SubSpaceInfo, label_symbol::Symbol)
    expanded_index = findfirst(x -> x == label_symbol, info.inner_labels_symbols_flat) 
    if expanded_index === nothing 
        error("Label $label_symbol not found in SubSpaceInfo.inner_labels_symbols")
        return nothing 
    else 
        return (info.outer_ss_of_expanded[expanded_index], info.inner_ss_of_expanded[expanded_index], expanded_index)
    end
    return nothing 
end
label_2_outer_inner_expanded(info::SubSpaceInfo, label::String) = label_2_outer_inner_expanded(info, Symbol(label))


# Specifies a subsystem location
struct SubSpaceIndex
    outer::Int      # subspace index
    inner::Int      # index in ensemble
    expanded::Int  
    function SubSpaceIndex(outer::Int, inner::Int, expanded::Int)
        return new(outer, inner, expanded)
    end
    function SubSpaceIndex(expanded::Int, info::SubSpaceInfo)
        outer,inner = expanded_2_outer_inner(info, expanded)
        return new(outer, inner, expanded)
    end
    function SubSpaceIndex(outer::Int, inner::Int, info::SubSpaceInfo)
        expanded = outer_inner_2_expanded(info, outer, inner)
        return new(outer, inner, expanded)
    end
    function SubSpaceIndex(label::Union{String,Symbol}, info::SubSpaceInfo)      
        outer, inner, expanded = label_2_outer_inner_expanded(info, label)
        return new(outer, inner, expanded) 
    end
end
@inline outer(i::SubSpaceIndex)    = i.outer
@inline inner(i::SubSpaceIndex)    = i.inner
@inline expanded(i::SubSpaceIndex) = i.expanded
@inline outer(is::Vector{SubSpaceIndex}) = [outer(i) for i in is]
@inline inner(is::Vector{SubSpaceIndex}) = [inner(i) for i in is]
@inline expanded(is::Vector{SubSpaceIndex}) = [expanded(i) for i in is]
@inline Index2Symbol(i::SubSpaceIndex, info::SubSpaceInfo) =  info.inner_labels_symbols_flat[i.expanded]
@inline Index2String(i::SubSpaceIndex, info::SubSpaceInfo) =  info.param_names[i.expanded]
@inline function Index2Ensemble(i::SubSpaceIndex, info::SubSpaceInfo) 
    ensemble = info.ensemble_index_by_outer_index[i.outer] # shouldn't be zero, otherwise not ensemble index
    @assert ensemble != 0 "index $i not an ensemble index" 
    return ensemble
end
# return index among current ensembles summation indexes 
@inline function Index2Summation(i::SubSpaceIndex, info::SubSpaceInfo) 
    ensemble = info.ensemble_index_by_outer_index[i.outer] 
    @assert ensemble != 0 "index $i not an ensemble index" 
    return i.inner - info.how_many_non_sum_by_ensemble[ensemble]
end
@inline function Index2Ensemble_and_Summation(i::SubSpaceIndex, info::SubSpaceInfo) 
    ensemble = info.ensemble_index_by_outer_index[i.outer] 
    @assert ensemble != 0 "index $i not an ensemble index" 
    return ensemble, i.inner - info.how_many_non_sum_by_ensemble[ensemble]
end

function Base.isless(a::SubSpaceIndex, b::SubSpaceIndex)::Bool
    return a.expanded < b.expanded
end
function Base.isequal(a::SubSpaceIndex, b::SubSpaceIndex)::Bool 
    return a.expanded == b.expanded  # is sufficient
end