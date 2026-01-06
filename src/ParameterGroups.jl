module ParameterGroups

export ParameterGroupKind, ParameterGroup, ParameterGroupLike, ParameterGroupScalar, ParameterGroupTimeScalar, ParameterGroupTimeFunction,
       ParameterGroupDistribution, ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction,
       ParameterGroupStorageUnion, parameter_group_input_type, parameter_group_storage_type, parameter_group_storage_target_type, parameter_group_value_type, parameter_group_kind_name, 
       WhereWhichParamGroup, AbstractEnsemble, AbstractSubSpaceInfo

using ..Sampler: QDistribution, QEnsembleFunction
using ..StringUtils: symbol2formatted, normalize_underscore_indices, format_normalized_indices
import Base: show

"""
    ParameterGroupKind

Enumerates the supported parameter-group forms used across the `QSpace`
pipeline.  The kind determines how a group's payload is validated and how its
values are materialised:

  * `ParameterGroupScalar` – literal value or array shared by all indices.
  * `ParameterGroupTimeScalar` – dedicated storage for explicit time parameters.
  * `ParameterGroupTimeFunction` – scalar function of time.
  * `ParameterGroupDistribution` – distribution backing ensemble sampling.
  * `ParameterGroupEnsembleFunction` – function that maps distribution samples to concrete values.
  * `ParameterGroupEnsembleTimeFunction` – ensemble function with explicit time
    dependence.
"""
@enum ParameterGroupKind::UInt8 begin
    ParameterGroupTimeScalar
    ParameterGroupScalar
    ParameterGroupTimeFunction
    ParameterGroupDistribution
    ParameterGroupEnsembleFunction
    ParameterGroupEnsembleTimeFunction
end

const PARAMETER_GROUP_PAYLOAD_INPUTS = (
    Union{Nothing, Number},                         # Time scalar (converted to vector)
    Union{Nothing, Number},                         # Scalar
    Union{Nothing, Function},                       # Time function
    Union{Nothing, QDistribution},                  # Distribution
    Union{Function, QEnsembleFunction},    # Ensemble function --> removed Nothing as an option, define at construction
    Union{Function, QEnsembleFunction})    # Ensemble time function --> removed Nothing as an option, define at construction


const PARAMETER_GROUP_PAYLOAD_TYPES = (
    Union{Nothing, Vector{Float64}},
    Union{Nothing, Float64},
    Union{Nothing, Function},
    Union{Nothing, QDistribution},
    Union{Nothing, QEnsembleFunction},
    Union{Nothing, QEnsembleFunction})

const PARAMETER_GROUP_PAYLOAD_TARGET_TYPES = (
    Vector{Float64},
    Float64,
    Function,
    QDistribution,
    QEnsembleFunction,
    QEnsembleFunction)

const PARAMETER_GROUP_VALUE_TYPES = (
    Vector{Float64},
    Float64,
    Vector{Float64},
    Vector{Float64},
    Array{Float64},
    Array{Float64})

const PARAMETER_GROUP_KIND_NAMES = (
    "TimeScalar",
    "Scalar",
    "TimeFunction",
    "Distribution",
    "EnsembleFunction",
    "EnsembleTimeFunction")

abstract type AbstractSubSpace end
abstract type AbstractEnsemble end
abstract type AbstractSubSpaceInfo end

const ParameterGroupStorageUnion = Union{PARAMETER_GROUP_PAYLOAD_TYPES...}

@inline parameter_group_input_type(kind::ParameterGroupKind) = PARAMETER_GROUP_PAYLOAD_INPUTS[Int(kind) + 1]
@inline parameter_group_storage_type(kind::ParameterGroupKind) = PARAMETER_GROUP_PAYLOAD_TYPES[Int(kind) + 1]
@inline parameter_group_storage_target_type(kind::ParameterGroupKind) = PARAMETER_GROUP_PAYLOAD_TARGET_TYPES[Int(kind) + 1]
@inline parameter_group_value_type(kind::ParameterGroupKind) = PARAMETER_GROUP_VALUE_TYPES[Int(kind) + 1]
@inline parameter_group_kind_name(kind::ParameterGroupKind) = PARAMETER_GROUP_KIND_NAMES[Int(kind) + 1]

"""
    ParameterGroup{T}

Shared metadata for a parameter group.  Skeleton instances are created during
calls so `ParameterDefinitions` and completed upon calls to QSpace.
"""
mutable struct ParameterGroup{T}
    param_symbol::Symbol
    param_raw::String
    param_str::String
    param_latex::String
    kind::ParameterGroupKind
    of_t::Bool
    indices::Vector{String}
    function_args::Vector{String}
    dependency_names::Vector{String}
    dependency_indices::Vector{Int}
    ensemble_indices::Vector{Int}
    unique_ensemble_indices::Vector{Int}
    subspace_indices::Vector{Int}
    index_symbol_pairs::Vector{Tuple{Symbol,Symbol}}
    index_string_pairs::Vector{Tuple{String,String}}
    index_sizes::Vector{Int}
    sample_sizes::Vector{Int}
    is_time_group::Bool
    payload::T
    function ParameterGroup(param_symbol::Symbol,
                            param_raw::String,
                            param_str::String,
                            param_latex::String,
                            kind::ParameterGroupKind,
                            of_t::Bool,
                            indices::Vector{String},
                            function_args::Vector{String},
                            dependency_names::Vector{String},
                            dependency_indices::Vector{Int},
                            ensemble_indices::Vector{Int},
                            unique_ensemble_indices::Vector{Int},
                            subspace_indices::Vector{Int},
                            index_symbol_pairs::Vector{Tuple{Symbol,Symbol}},
                            index_string_pairs::Vector{Tuple{String,String}},
                            index_sizes::Vector{Int},
                            sample_sizes::Vector{Int},
                            is_time_group::Bool,
                            payload)
        storage_type = parameter_group_storage_type(kind)
        payload isa storage_type || throw(ArgumentError("Payload for $(kind) must be of type $(storage_type), got $(typeof(payload))."))
        return new{storage_type}(param_symbol,
                                 param_raw,
                                 param_str,
                                 param_latex,
                                 kind,
                                 of_t,
                                 indices,
                                 function_args,
                                 dependency_names,
                                 dependency_indices,
                                 ensemble_indices,
                                 unique_ensemble_indices,
                                 subspace_indices,
                                 index_symbol_pairs,
                                 index_string_pairs,
                                 index_sizes,
                                 sample_sizes,
                                 is_time_group,
                                 payload)
    end
end

const _PARAMETER_GROUP_TYPES = ntuple(i -> ParameterGroup{PARAMETER_GROUP_PAYLOAD_TYPES[i]}, length(PARAMETER_GROUP_PAYLOAD_TYPES))
const ParameterGroupLike = Union{_PARAMETER_GROUP_TYPES...}

function _format_index_pairs(pairs::Vector{Tuple{String,String}})::String
    isempty(pairs) && return "[]"
    formatted = String[]
    for (non_str, sum_str) in pairs
        if isempty(sum_str)
            push!(formatted, "(" * non_str * ")")
        else
            push!(formatted, "(" * non_str * "," * sum_str * ")")
        end
    end
    return "[" * join(formatted, ", ") * "]"
end

@inline function _payload_type(group::ParameterGroup)
    payload = group.payload
    payload === nothing && return "unset"
    return String(nameof(typeof(payload)))
end

@inline function _format_group_argument(arg::String)::String
    arg == "t" && return "t"
    base, idxs = normalize_underscore_indices(arg)
    base_str = symbol2formatted(base)[1]
    return base_str * format_normalized_indices(idxs; do_latex=false)
end

function _format_group_signature(group::ParameterGroupLike)::String
    base_str = symbol2formatted(group.param_raw)[1]
    if !isempty(group.indices)
        base_str *= format_normalized_indices(group.indices; do_latex=false)
    end
    args = String[]
    if !isempty(group.function_args)
        for arg in group.function_args
            push!(args, _format_group_argument(arg))
        end
    elseif group.of_t && group.param_symbol != :t
        push!(args, "t")
    end
    if !isempty(args)
        base_str *= "(" * join(args, ",") * ")"
    end
    return base_str
end

function _single_line_summary(group::ParameterGroupLike)::String
    entries = String[]
    push!(entries, "indices=" * _format_index_pairs(group.index_string_pairs))
    push!(entries, "deps=" * "[" * join(group.dependency_names, ", ") * "]")
    push!(entries, "payload=" * _payload_type(group))
    if group.is_time_group && group.param_symbol != :t
        push!(entries, "time-group")
    end
    return string(parameter_group_kind_name(group.kind), ": ", _format_group_signature(group), " (", join(entries, ", "), ")")
end

function show(io::IO, group::ParameterGroup{T}) where {T}
    signature = _format_group_signature(group)
    if get(io, :compact, false)
        print(io, signature)
        return
    end
    print(io, "ParameterGroup(", parameter_group_kind_name(group.kind), "): ", signature)
    entries = Pair{String,String}[]
    push!(entries, "indices" => _format_index_pairs(group.index_string_pairs))
    push!(entries, "deps" => "[" * join(group.dependency_names, ", ") * "]")
    payload_str = _payload_type(group)
    push!(entries, "payload" => payload_str)
    if group.is_time_group && group.param_symbol != :t
        push!(entries, "time-group" => "true")
    end
    for (key, value) in entries
        print(io, "\n  ", key, "=", value)
    end
end

function show(io::IO, ::MIME"text/plain", groups::Vector{<:ParameterGroupLike})
    print(io, "ParameterGroups(", length(groups), ")")
    for group in groups
        print(io, "\n  ", _single_line_summary(group))
    end
end

struct WhereWhichParamGroup
    scalar_groups::Vector{Int}
    time_function_groups::Vector{Int}
    distribution_groups::Vector{Int}
    ensemble_function_groups::Vector{Int}
    time_group::Int
    first_qensemble_group::Int
end

function WhereWhichParamGroup(groups::Vector{ParameterGroupLike})
    scalar_groups = Int[]
    time_function_groups = Int[]
    distribution_groups = Int[]
    ensemble_function_groups = Int[]
    time_group = 0
    first_qensemble_group = 0
    for (idx, group) in enumerate(groups)
        group.is_time_group && (time_group = idx)
        kind = group.kind
        if kind == ParameterGroupScalar || kind == ParameterGroupTimeScalar
            push!(scalar_groups, idx)
        elseif kind == ParameterGroupTimeFunction
            push!(time_function_groups, idx)
        elseif kind == ParameterGroupDistribution
            push!(distribution_groups, idx)
        elseif kind == ParameterGroupEnsembleFunction || kind == ParameterGroupEnsembleTimeFunction
            push!(ensemble_function_groups, idx)
            first_qensemble_group == 0 && (first_qensemble_group = idx)
        end
    end
    first_qensemble_group == 0 && (first_qensemble_group = length(groups) + 1)
    return WhereWhichParamGroup(scalar_groups, time_function_groups, distribution_groups,
                                ensemble_function_groups, time_group, first_qensemble_group)
end

end # module ParameterGroups
