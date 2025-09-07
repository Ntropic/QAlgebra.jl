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
    particle_type::String
    op_set::OperatorSet             # The operator set for this subspace.
end
function SubSpace(key_symbol::Symbol, ss_outer_ind::Int, index_counter::Int, is_ensemble_ss::Bool, ensemble_size::Int, op_set::OperatorSet)::Tuple{SubSpace, Vector{Symbol}}
    key::String = String(key_symbol)
    return SubSpace(key, ss_outer_ind, index_counter, is_ensemble_ss, ensemble_size, op_set) 
end
function SubSpace(key::String, ss_outer_ind::Int, index_counter::Int, is_ensemble_ss::Bool, ensemble_size::Int, op_set::OperatorSet, core_keys::Vector{Symbol})::Tuple{SubSpace, Vector{Symbol}}
    if length(key) > 1 && ensemble_size > 1
        error("Currently only supports single characters for spin baths!")
    end
    if key== "t"
        error("Letter t is reserved for time, cannot be used for baths. ")
    end
    if ensemble_size < 1 
        error("For subspace pattern key \"$key\", the ensemble_size must be a positive integer. ")
    end
    start_char = Char(key[1])
    # first attempt at keys -> subsequent letters 
    keys::String = [string(start_char + i) for i in 0:(ensemble_size-1)]
    keys_symbols::Vector{Symbol} = [Symbol(key) for key in keys]
    # check if any of the keys_symbols are in core_keys   =========> Alternative to i,j,k,... indexing i0, i1, i2 and so on.
    if any(keys_symbol -> keys_symbol in core_keys, keys_symbols) 
        keys = [string(start_char)*"_"*string(i)  for i in 0:(ensemble_size-1)]
        keys_symbols = [Symbol(key) for key in keys]
    end
    key_symbol::Symbol = keys_symbol[1]
    ss_inner_ind::Vector{Int} = [index_counter+i for i in 1:ensemble_size]
    particle_type::String = OperatorSet.particle_type 
    return SubSpace(key_symbol, keys_symbol, key, keys, ss_outer_ind, ss_inner_ind, is_ensemble_ss, ensemble_size, particle_type, op_set), keys_symbols 
end
# Define the custom show for SubSpace.
function Base.show(io::IO, statespace::SubSpace)
    # Print the subspace key and allowed keys.
    print(io, "SubSpace ", statespace.keys, ": ")
    # Use the OperatorSet's show for the op_set field.
    show(io, statespace.op_set)
end

""" 
    SubSpaceDefinitions(; kwargs...)

SubSpaceDefinitions is a struct that processes the definition of the quantum subspaces, 
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
                ensemble_size, op_set = val  # unpacking
            else
                ensemble_size, op_set = 1, val 
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
            curr_subspace = SubSpace(key_symbol, keys_symbols, key, keys, outer_ind, curr_inds, is_ensemble_ss, ensemble_size, op_set.particle_type, op_set) 
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
    ensemble_indexes::Vector{Vector{Int}}      # The expanded indexes of the ensembles 
    of_time::Bool
end
# Primary constructor from labels
function SubSpaceInfo(outer_labels_symbols::Vector{Symbol}, inner_labels_symbols::Vector{Vector{Symbol}}, are_ensemble_ss::Vector{Bool}, of_time::Bool=false)
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
    return SubSpaceInfo( outer_labels_symbols, inner_labels_symbols, inner_labels_symbols_flat, 
                         outer_labels, inner_labels, param_names, 
                         subsystem_sizes, outer_ss_of_expanded, inner_ss_of_expanded, expanded_index_by_outer,
                         where_ensembles, ensemble_index_by_outer_index, how_many_by_ensemble, ensemble_indexes, of_time )
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
    return SubSpaceInfo(outers, inners, are_ensemble_ss)
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

function Base.isless(a::SubSpaceIndex, b::SubSpaceIndex)::Bool
    return a.expanded < b.expanded
end
function Base.isequal(a::SubSpaceIndex, b::SubSpaceIndex)::Bool 
    return a.expanded == b.expanded  # is sufficient
end