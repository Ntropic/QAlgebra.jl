import ..SubSpaceIndex
import ..EnsembleIndex
using ..ParameterGroups: AbstractSubSpace, AbstractEnsemble, AbstractSubSpaceInfo
using ..EnsembleSamples: AbstractEnsembleSample, DiscreteSamples, ContinuousSamples
using ..QAlgebra: SAMPLE_ATOL, SAMPLE_RTOL
using ..QIndexes: AbstractIndex
using ..StringUtils: str2sub, symbol2formatted

"""
    Ensemble(operator_set; sum_label=nothing, kwargs...)
    Ensemble(sum_label, operator_set; kwargs...)

Construct an ensemble subspace. The primary (non-summation) label is taken from the
owning subspace's key; providing `sum_label` supplies an explicit summation label.
When omitted, a sensible default is inferred.

# Arguments
- `operator_set::OperatorSet`: operator algebra for the ensemble.
- `sum_label::Union{Nothing,Symbol,String}`: optional summation label (positional or keyword).

# Keyword arguments
- `num_modes::Int = -1`: physical number of modes (`-1` reserves indices only).
- `as_continuum::Bool = false`: flag ensemble as continuum.
- `max_operator_order::Int = -1`: maximum operator order (`-1` for unbounded).
- `param_groups::Vector{Symbol} = Symbol[]`: associated parameter groups.
- `distribution_group_indices::Vector{Int} = Int[]`: distribution group indices.
- `ensemble_function_group_indices::Vector{Int} = Int[]`: ensemble function groups.

- `sample_method::Symbol = :default`: preferred sampling strategy.
- `sample_num_nodes::Int = 25`: interpolation nodes for sampling.
- `sample_atol::Float64 = SAMPLE_ATOL`: absolute sampling tolerance.
- `sample_rtol::Float64 = SAMPLE_RTOL`: relative sampling tolerance.
- `sample_max_iter::Int = 128`: maximum refinement iterations.
"""
mutable struct Ensemble <: AbstractEnsemble
    non_sum_symbol::Symbol
    non_sum_string::String
    sum_symbol::Symbol
    sum_string::String
    operator_set::OperatorSet
    num_modes::Int
    max_operator_order::Int
    as_continuum::Bool
    param_groups::Vector{Symbol}
    distribution_group_indices::Vector{Int}
    ensemble_function_group_indices::Vector{Int}
    sample_method::Symbol
    sample_num_nodes::Int
    sample_atol::Float64
    sample_rtol::Float64
    sample_max_iter::Int
    sampler::Union{Nothing,AbstractEnsembleSample}
    qspace_ref::Union{Nothing,WeakRef}
    function Ensemble(operator_set::OperatorSet; sum_label::Union{Nothing,Symbol,String}=nothing,
                      num_modes::Int=-1,
                      max_operator_order::Int=-1,
                      as_continuum::Bool=false,
                      param_groups::Vector{Symbol}=Symbol[],
                      distribution_group_indices::Vector{Int}=Int[],
                      ensemble_function_group_indices::Vector{Int}=Int[],
                      sample_method::Symbol=:default,
                      sample_num_nodes::Int=25,
                      sample_atol::Float64=SAMPLE_ATOL,
                      sample_rtol::Float64=SAMPLE_RTOL,
                      sample_max_iter::Int=128,
                      sampler::Union{Nothing,AbstractEnsembleSample}=nothing,
                      qspace_ref::Union{Nothing,WeakRef}=nothing)
        sum_symbol = isnothing(sum_label) ? Symbol("") : Symbol(sum_label)
        sum_string = isnothing(sum_label) ? "" : String(sum_label)
        return new(Symbol(""), "", sum_symbol, sum_string, operator_set, num_modes,
                   max_operator_order, as_continuum, copy(param_groups), copy(distribution_group_indices),
                   copy(ensemble_function_group_indices), sample_method, sample_num_nodes,
                   sample_atol, sample_rtol, sample_max_iter, sampler, qspace_ref)
    end
end


function Ensemble(sum_label::Union{Symbol,String}, operator_set::OperatorSet; kwargs...)
    return Ensemble(operator_set; sum_label=sum_label, kwargs...)
end


function Base.show(io::IO, ensemble::Ensemble)
    non_sum = isempty(ensemble.non_sum_string) ? "⟂" : ensemble.non_sum_string
    sum_lbl = isempty(ensemble.sum_string) ? "?" : ensemble.sum_string
    print(io, "Ensemble[", non_sum, ":", sum_lbl, "]")
    print(io, " ⟨", ensemble.operator_set.name, "⟩")
    print(io, ", modes=" , ensemble.num_modes)
    if ensemble.max_operator_order != -1
        print(io, ", max_order=" , ensemble.max_operator_order)
    end
    if !isempty(ensemble.param_groups)
        print(io, ", param_groups=" , ensemble.param_groups)
    end
    if !isempty(ensemble.distribution_group_indices)
        print(io, ", distribution_group_indices=" , ensemble.distribution_group_indices)
    end
    if !isempty(ensemble.ensemble_function_group_indices)
        print(io, ", ensemble_function_group_indices=" , ensemble.ensemble_function_group_indices)
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
struct SubSpace <: AbstractSubSpace
    key_symbol::Symbol
    key::String
    key_symbol_summation::Symbol
    key_summation::String
    is_ensemble_ss::Bool
    ensemble_size::Int
    as_continuum::Bool
    particle_type::String
    op_set::OperatorSet
    ensemble::Union{Nothing,Ensemble}
    min_ints::Vector{Int}
    max_ints::Vector{Int}
    max_operator_magnitude::Int
end

function Base.show(io::IO, sub::SubSpace)
    if sub.is_ensemble_ss && !isempty(sub.key_summation)
        print(io, "SubSpace: ", sub.key, " ∑ ", sub.key_summation, " (", sub.op_set.name, ")")
    else
        print(io, "SubSpace: ", sub.key, " (", sub.op_set.name, ")")
    end
end

function AbstractIndex2string(subspaces::Vector{SubSpace}, index::AbstractIndex; do_latex::Bool=false, as_index::Bool=true)::Bool 
    ensemble = subspaces[index.subspace].ensemble
    index_symbol = index.summation ? ensemble.sum_string : ensemble.non_sum_string
    subindex = index.slot == 0 ? "" : String(index.slot)
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



function _next_symbol(base::Symbol)
    str = String(base)
    if length(str) == 1
        ch = str[1]
        if 'a' <= ch < 'z'
            nxt = Char(ch + 1)
            return Symbol(nxt), string(nxt)
        end
    end
    error("Cannot infer summation label from $(base); please supply it explicitly.")
end

function _next_available_symbol(base::Symbol, used::Set{Symbol})
    sym, str = _next_symbol(base)
    while sym in used
        sym, str = _next_symbol(sym)
    end
    return sym, str
end

@inline has_summation(sub::SubSpace) = sub.is_ensemble_ss && !isempty(sub.key_summation)
@inline primary_symbol(sub::SubSpace) = sub.key_symbol
@inline secondary_symbol(sub::SubSpace) = sub.key_symbol_summation
@inline primary_label(sub::SubSpace) = sub.key
@inline secondary_label(sub::SubSpace) = sub.key_summation

function subspace_symbols(sub::SubSpace)
    if has_summation(sub)
        return Symbol[primary_symbol(sub), secondary_symbol(sub)]
    else
        return Symbol[primary_symbol(sub)]
    end
end

function subspace_labels(sub::SubSpace)
    if has_summation(sub)
        return String[primary_label(sub), secondary_label(sub)]
    else
        return String[primary_label(sub)]
    end
end

"""
    SubSpaceDefinitions(; kwargs...)

Collect all subspace declarations for a [`QSpace`](@ref). Each keyword maps a symbolic
identifier to either an `OperatorSet` (single subspace) or an [`Ensemble`](@ref)
configuration. The constructor derives internal labels, ensemble metadata, and neutral
elements for each subspace.
"""
struct SubSpaceDefinitions 
    subspaces::Vector{SubSpace}
    used_symbols::Set{Symbol}

    function SubSpaceDefinitions(; kwargs...)
        subspaces = SubSpace[]
        used_symbols = Set{Symbol}()

        for (key_symbol, val) in kwargs
            op_set, ensemble_cfg, is_ensemble_ss, as_continuum = begin
                if val isa Ensemble
                    (val.operator_set, val, true, val.as_continuum)
                elseif val isa OperatorSet
                    (val, nothing, false, false)
                else
                    error("Invalid subspace definition for $key_symbol: expected an OperatorSet or Ensemble, got $(typeof(val)).")
                end
            end

            isa(op_set, OperatorSet) || error("Invalid subspace definition for $key_symbol: expected an OperatorSet.")

            non_sym = key_symbol
            non_str = String(key_symbol)
            sum_sym = Symbol("")
            sum_str = ""
            ensemble_size = 1

            if is_ensemble_ss
                @assert ensemble_cfg !== nothing
                if isempty(ensemble_cfg.non_sum_string)
                    ensemble_cfg.non_sum_symbol = non_sym
                    ensemble_cfg.non_sum_string = non_str
                end
                non_sym = ensemble_cfg.non_sum_symbol
                non_str = ensemble_cfg.non_sum_string

                if isempty(ensemble_cfg.sum_string)
                    sum_sym, sum_str = _next_available_symbol(non_sym, used_symbols)
                    ensemble_cfg.sum_symbol = sum_sym
                    ensemble_cfg.sum_string = sum_str
                else
                    sum_sym = ensemble_cfg.sum_symbol
                    sum_str = ensemble_cfg.sum_string
                end
                ensemble_size = isempty(sum_str) ? 1 : 2
            end

            non_sym in used_symbols && error("Symbol $(non_sym) already used")
            push!(used_symbols, non_sym)
            if !isempty(sum_str)
                sum_sym in used_symbols && error("Symbol $(sum_sym) already used")
                push!(used_symbols, sum_sym)
            end

            max_op_mag = if is_ensemble_ss && ensemble_cfg.max_operator_order != -1
                ensemble_cfg.max_operator_order
            else
                max_operator_magnitude(op_set)
            end

            push!(subspaces, SubSpace(non_sym, non_str, sum_sym, sum_str, is_ensemble_ss, ensemble_size,
                                      as_continuum, op_set.particle_type, op_set, ensemble_cfg,
                                      copy(op_set.min_ints), copy(op_set.max_ints), max_op_mag))
        end

        return new(subspaces, used_symbols)
    end
end
function Base.show(io::IO, subspace_def::SubSpaceDefinitions)
    println(io, "SubSpaceDefinitions: ")
    # Then print each subspace on its own line.
    for ss in subspace_def.subspaces
        println(io, "   - ", string(ss))
    end
end

struct SubSpaceInfo <: AbstractSubSpaceInfo
    subspaces::Vector{SubSpace}
    where_ensembles::Vector{Int}
    ensemble_index_by_subspace_index::Vector{Int}
    ensemble_non_sum_labels::Vector{Symbol}
    ensemble_sum_labels::Vector{Symbol}
    ensemble_label_symbols::Vector{Vector{Symbol}}
    ensemble_labels::Vector{Vector{String}}
    of_time::Bool
end

function SubSpaceInfo(subspaces::Vector{SubSpace}; of_time::Bool=false)
    where_ensembles = Int[]
    ensemble_index_by_subspace_index = zeros(Int, length(subspaces))
    ensemble_non_sum_labels = Symbol[]
    ensemble_sum_labels = Symbol[]
    ensemble_label_symbols = Vector{Symbol}[]
    ensemble_labels = Vector{String}[]

    for (outer, subspace) in enumerate(subspaces)
        if subspace.is_ensemble_ss
            push!(where_ensembles, outer)
            ensemble_index_by_subspace_index[outer] = length(where_ensembles)
            ens = subspace.ensemble
            ens === nothing && error("Ensemble metadata missing for subspace $(outer).")

            non_sym = isempty(ens.non_sum_string) ? subspace.key_symbol : ens.non_sum_symbol
            sum_sym = isempty(ens.sum_string) ? non_sym : ens.sum_symbol
            non_str = isempty(ens.non_sum_string) ? subspace.key : ens.non_sum_string
            sum_str = isempty(ens.sum_string) ? non_str : ens.sum_string

            push!(ensemble_non_sum_labels, non_sym)
            push!(ensemble_sum_labels, sum_sym)
            push!(ensemble_label_symbols, Symbol[non_sym, sum_sym])
            push!(ensemble_labels, String[non_str, sum_str])
        end
    end

    return SubSpaceInfo(subspaces, where_ensembles, ensemble_index_by_subspace_index,
                        ensemble_non_sum_labels, ensemble_sum_labels, ensemble_label_symbols, ensemble_labels, of_time)
end

function Base.show(io::IO, info::SubSpaceInfo)
    println(io, "SubSpaceInfo:")
    println(io, "  subspaces:      ", length(info.subspaces))
    println(io, "  ensembles:      ", info.where_ensembles)
    println(io, "  non-sum labels: ", info.ensemble_non_sum_labels)
    println(io, "  sum labels:     ", info.ensemble_sum_labels)
end

# Legacy helpers are intentionally removed. Downstream code must transition to the
# new particle-based indexing before these utilities can be reintroduced.
function outer_inner_2_expanded(args...)
    error("outer_inner_2_expanded was removed during the subspace refactor. Update call sites to the new API.")
end

function expanded_2_outer_inner(args...)
    error("expanded_2_outer_inner was removed during the subspace refactor. Update call sites to the new API.")
end

function label_2_outer_inner_expanded(args...)
    error("label_2_outer_inner_expanded was removed during the subspace refactor. Update call sites to the new API.")
end

function SubSpaceIndex(args...)
    error("SubSpaceIndex helpers were removed during the subspace refactor. Update call sites to the new API.")
end
