import ..SubSpaceIndex
using ..EnsembleSamples: AbstractEnsembleSample, DiscreteSamples, ContinuousSamples

"""
    Ensemble(num_operator_indexes, num_sum_indexes, operator_set; kwargs...)
    Ensemble(num_operator_indexes, operator_set; kwargs...)
    Ensemble(; num_operator_indexes, num_sum_indexes=0, operator_set, kwargs...)

Create an `Ensemble`: a container for ensemble–subspace metadata.

# Positional arguments
- `num_operator_indexes::Int`: Number of operator indices reserved for the ensemble.
- `num_sum_indexes::Int`: Number of summation indices reserved for the ensemble.   (Omitted and defaulted to `0` in the 2-arg / keyword-only constructors.)
- `operator_set::OperatorSet`: The operator set associated with this ensemble.

# Keyword arguments
- `num_modes::Int = -1`: Physical number of instantiated modes. `-1` means the ensemble only reserves indices (no fixed system size).
- `as_continuum::Bool = false`: Whether a continuum approximation is planned.
- `max_operator_order::Int = -1`: Maximum operator order allowed (convention: `-1` = unbounded).
- `sample_method::Symbol = :default`: Preferred sampling method. `:default` resolves to `:random` for discrete ensembles and `:chebychev` for continuum ensembles. 
    Discrete sampling accepts `:random` or `:density`;
    Continuum sampling accepts `:chebychev`, `:uniform`, `:leja`, or `:fekete` (also `:chebyshev` alias).

- `sample_num_nodes::Int = 25`: Number of interpolation nodes used for CDF approximations.
- `sample_atol::Float64 = 1e-9` Absolute tolerance reused by sampling helpers.
- `sample_rtol::Float64 = 1e-7`: Relative tolerance reused by sampling helpers.
- `sample_max_iter::Int = 128`: Maximum refinement iterations for inverse-CDF halving-steps.
"""
mutable struct Ensemble
    num_operator_indexes::Int
    num_sum_indexes::Int
    operator_set::OperatorSet
    num_modes::Int
    max_operator_order::Int
    as_continuum::Bool
    parameter_groups::Vector{Symbol}
    parameter_group_indices::Vector{Int}
    parameter_function_group_indices::Vector{Int}
    sample_method::Symbol
    sample_num_nodes::Int
    sample_atol::Float64
    sample_rtol::Float64
    sample_max_iter::Int
    sampler::Union{Nothing,AbstractEnsembleSample}
    qspace_ref::Union{Nothing,WeakRef}
    function Ensemble(num_operator_indexes::Int, num_sum_indexes::Int, operator_set::OperatorSet;
                      num_modes::Int=-1,
                      max_operator_order::Int=-1,
                      as_continuum::Bool=false,
                      parameter_groups::Vector{Symbol}=Symbol[],
                      parameter_group_indices::Vector{Int}=Int[],
                      parameter_function_group_indices::Vector{Int}=Int[],
                      sample_method::Symbol=:default,
                      sample_num_nodes::Int=25,
                      sample_atol::Float64=1e-9,
                      sample_rtol::Float64=1e-7,
                      sample_max_iter::Int=128,
                      sampler::Union{Nothing,AbstractEnsembleSample}=nothing,
                      qspace_ref::Union{Nothing,WeakRef}=nothing)
        return new(num_operator_indexes, num_sum_indexes, operator_set, num_modes,
                   max_operator_order, as_continuum, copy(parameter_groups), copy(parameter_group_indices),
                   copy(parameter_function_group_indices), sample_method, sample_num_nodes,
                   sample_atol, sample_rtol, sample_max_iter, sampler, qspace_ref)
    end
    function Ensemble(num_operator_indexes::Int, operator_set::OperatorSet; kwargs...)
        return Ensemble(num_operator_indexes, 0, operator_set; kwargs...)
    end
    function Ensemble(; num_operator_indexes::Int, num_sum_indexes::Int=0, operator_set::OperatorSet, kwargs...)
        return Ensemble(num_operator_indexes, num_sum_indexes, operator_set; kwargs...)
    end
end

function Base.show(io::IO, ensemble::Ensemble)
    if get(io, :compact, false)
        # --- compact version ---
        kind = ensemble.sampler === nothing ?
            (ensemble.sample_method === :default ?
                (ensemble.as_continuum ? ":chebychev" : ":random") :
                string(ensemble.sample_method)) :
            (ensemble.sampler isa ContinuousSamples ? ":cont/" : ":disc/") * string(ensemble.sampler.method)
        print(io, "Ensemble(", ensemble.num_operator_indexes, "op,", ensemble.num_sum_indexes, "∑,",
              ensemble.num_modes, "modes, method=", kind, ")")
        return
    end
    print(io, "Ensemble: ")
    print(io, "operators=" , ensemble.num_operator_indexes)
    print(io, ", summations=" , ensemble.num_sum_indexes)
    print(io, ", num_modes=" , ensemble.num_modes)
    if ensemble.max_operator_order != -1
        print(io, ", max_order=" , ensemble.max_operator_order)
    end
    if !isempty(ensemble.parameter_groups)
        print(io, ", parameter_groups=" , ensemble.parameter_groups)
    end
    if !isempty(ensemble.parameter_group_indices)
        print(io, ", parameter_group_indices=" , ensemble.parameter_group_indices)
    end
    if !isempty(ensemble.parameter_function_group_indices)
        print(io, ", parameter_function_group_indices=" , ensemble.parameter_function_group_indices)
    end
    if ensemble.sampler !== nothing
        sample = ensemble.sampler
        kind = sample isa ContinuousSamples ? ":continuous" : ":discrete"
        print(io, ", sample_method=" , kind, "/", sample.method)
    else
        resolved = ensemble.sample_method === :default ?
            (ensemble.as_continuum ? :chebychev : :random) : ensemble.sample_method
        print(io, ", sample_method=" , resolved)
    end
end

"""
    SubSpace(...)

Internal representation of a named subspace inside a [`QSpace`](@ref). Each subspace
records the operator set it draws from, multiplicities for ensemble replication, and
links back to the parent `Ensemble` when applicable. The public API constructs these
objects through [`SubSpaceDefinitions`](@ref).
"""
struct SubSpace
    key_symbol::Symbol
    keys_symbols::Vector{Symbol}
    key::String                     # Original input key
    keys::Vector{String}            # Allowed keys for this subspace 
    keys_latex::Vector{String}
    ss_outer_ind::Int            # Which Vector to use for ss_inner_ind  (this is for accessing the string elements)
    ss_inner_ind::Vector{Int}    # Indices to access operator values in the corresponding qspace main ind  (this is for accessing the string elements)
    is_ensemble_ss::Bool
    ensemble_size::Int    # how many modes 
    as_continuum::Bool
    num_operator_indexes::Int 
    num_sum_indexes::Int
    particle_type::String
    op_set::OperatorSet             # The operator set for this subspace.
    ensemble::Union{Nothing,Ensemble}
    min_ints::Vector{Int}
    max_ints::Vector{Int}           # -1 entries signal unbounded axes
    max_operator_magnitude::Int
end
# Define the custom show for SubSpace.
function Base.show(io::IO, qspace::SubSpace)
    # Print the subspace key and allowed keys.
    print(io, "SubSpace ", qspace.keys[1:qspace.num_operator_indexes], " , ∑", qspace.keys[qspace.num_operator_indexes+1:end], ": ")
    # Use the OperatorSet's show for the op_set field.
    show(io, qspace.op_set)
end

const _ALPHABET = [c for c in 'a':'z']

function _numeric_labels(base::String, num_op::Int, num_sum::Int)
    total = num_op + num_sum
    labels_symbol = Symbol[]
    labels = String[]
    labels_latex = String[]
    for idx in 0:(total-1)
        label_str = base * string(idx)
        push!(labels_symbol, Symbol(label_str))
        push!(labels, label_str)
        push!(labels_latex, base * "_{" * string(idx) * "}") 
    end
    return labels_symbol, labels, labels_latex
end

function _alphabetic_labels(base_char::Char, total::Int)
    labels = String[]
    start_idx = findfirst(==(base_char), _ALPHABET)
    if start_idx === nothing
        error("Not a valid index: $base_char. ")
    end
    idx = start_idx
    while length(labels) < total
        push!(labels, string(_ALPHABET[idx]))
        idx += 1
        if idx > length(_ALPHABET)
            error("Exceeded index range by reaching $(_ALPHABET[idx-1]). ")
        end
    end
    return labels
end

"""
    SubSpaceDefinitions(; kwargs...)

Collect all subspace declarations for a [`QSpace`](@ref). Each keyword maps a symbolic
identifier to either an `OperatorSet` (single subspace) or an [`Ensemble`](@ref)
configuration. The constructor derives internal labels, ensemble metadata, and neutral
elements for each subspace.
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
        reserved_outer = Set{String}(String.(core_keys))
        for (outer_ind, (key_symbol, val)) in enumerate(kwargs) 
            is_ensemble_ss = false
            as_continuum = false
            ensemble_cfg::Union{Nothing,Ensemble} = nothing
            if isa(val, Ensemble)
                ensemble_cfg = val
                is_ensemble_ss = true
                as_continuum = val.as_continuum
            elseif isa(val, OperatorSet)
                # handled below
            else
                error("Invalid subspace definition for $key_symbol: expected an OperatorSet or Ensemble, got $(typeof(val)).")
            end

            if is_ensemble_ss
                @assert ensemble_cfg !== nothing
                num_operator_indexes = ensemble_cfg.num_operator_indexes
                num_sum_indexes = ensemble_cfg.num_sum_indexes
                op_set = ensemble_cfg.operator_set
                ensemble_size = num_operator_indexes + num_sum_indexes
                max_op_mag = ensemble_cfg.max_operator_order != -1 ? ensemble_cfg.max_operator_order : max_operator_magnitude(op_set)
            else
                ensemble_size = 1
                op_set = val
                num_operator_indexes = 1
                num_sum_indexes = 0
                max_op_mag = max_operator_magnitude(op_set)
            end

            isa(op_set, OperatorSet) || error("Invalid subspace definition for $key_symbol: expected an OperatorSet or Ensemble.")
            key = String(key_symbol) 
            reserved_current = Set{String}(reserved_outer)
            delete!(reserved_current, key)
            if length(key) == 1
                key_char = key[1]
                keys = _alphabetic_labels(key_char, ensemble_size)
                keys_symbols = Symbol.(keys) 
                keys_latex = keys
                if any(sym-> sym in used_symbols, keys_symbols) || any(lbl -> lbl in reserved_current, keys)
                    keys_symbols, keys, keys_latex = _numeric_labels(key, num_operator_indexes, num_sum_indexes)
                end
            else 
                keys_symbols, keys, keys_latex = _numeric_labels(key, num_operator_indexes, num_sum_indexes)
            end
            if any(x->x in used_symbols, keys_symbols) 
                error("Symbol $key already used")
            end 
            curr_inds = key_counter .+ collect(1:ensemble_size)
            curr_subspace = SubSpace(key_symbol, keys_symbols, key, keys, keys_latex, outer_ind, curr_inds, is_ensemble_ss, 
                        ensemble_size, as_continuum, num_operator_indexes, num_sum_indexes, op_set.particle_type, op_set, ensemble_cfg,
                        copy(op_set.min_ints), copy(op_set.max_ints), max_op_mag) 
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
function SubSpaceInfo(outer_labels_symbols::Vector{Symbol}, inner_labels_symbols::Vector{Vector{Symbol}}, are_ensemble_ss::BitVector, 
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
    are_ensemble_ss::BitVector = [s.is_ensemble_ss for s in subspaces] 
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
#SubSpaceIndex(outer::Int, inner::Int, expanded::Int) = SubSpaceIndex(outer, inner, expanded)
function SubSpaceIndex(expanded::Int, info::SubSpaceInfo)
    outer, inner = expanded_2_outer_inner(info, expanded)
    SubSpaceIndex(outer, inner, expanded)
end
function SubSpaceIndex(outer::Int, inner::Int, info::SubSpaceInfo)
    expanded = outer_inner_2_expanded(info, outer, inner)
    SubSpaceIndex(outer, inner, expanded)
end
function SubSpaceIndex(label::Union{String,Symbol}, info::SubSpaceInfo)
    outer, inner, expanded = label_2_outer_inner_expanded(info, label)
    SubSpaceIndex(outer, inner, expanded)
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
@inline function SummationIndex2SubSpaceIndex(outer_ind::Int, ensemble_ind::Int, summation_ind::Int, info::SubSpaceInfo)::SubSpaceIndex
    inner_ind = summation_ind + info.how_many_non_sum_by_ensemble[ensemble_ind]
    return SubSpaceIndex(outer_ind, inner_ind, outer_inner_2_expanded(info, outer_ind, inner_ind))
end

@inline function Base.isless(a::SubSpaceIndex, b::SubSpaceIndex)::Bool
    return a.expanded < b.expanded
end
@inline function Base.isequal(a::SubSpaceIndex, b::SubSpaceIndex)::Bool 
    return a.expanded == b.expanded  # is sufficient
end
