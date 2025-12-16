module ParameterGroups

export ParameterGroupKind, ParameterGroup, ParameterGroupLike, ParameterGroupScalar, ParameterGroupTimeScalar, ParameterGroupTimeFunction,
       ParameterGroupDistribution, ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction,
       ParameterGroupStorageUnion, parameter_group_input_type, parameter_group_storage_type, parameter_group_storage_target_type, parameter_group_value_type, parameter_group_kind_name, 
       WhereWhichParamGroup, AbstractEnsemble, AbstractSubSpaceInfo

using ..Sampler: QDistribution, QEnsembleFunction

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
`ParameterDefinitions` and finalised inside `ParameterDefinitions2Parameters`.
Once assembled they are referenced by `ParameterInfo`, `ParameterValues`, and
each `Ensemble`.  The mutable `payload` field stores the
user-provided definition (distribution, function, literal value, …) and may be
updated after `QSpace` construction via `resolve_param!`.  The type parameter
`T` records the storage type for `payload` (e.g. `Union{Nothing,Number}` for
scalar groups or `Union{Nothing,QEnsembleFunction}` for ensemble functions),
allowing the compiler to reason precisely about group contents.

Fields capture the group's declarative signature (`name`, `indices`,
`function_args`), dependency tracking (`dependency_names`/`dependency_indices`),
ensemble affiliation (`ensemble_indices`), and derived shape information
(`index_sizes`).
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
                                 index_sizes,
                                 sample_sizes,
                                 is_time_group,
                                 payload)
    end
end

const _PARAMETER_GROUP_TYPES = ntuple(i -> ParameterGroup{PARAMETER_GROUP_PAYLOAD_TYPES[i]}, length(PARAMETER_GROUP_PAYLOAD_TYPES))
const ParameterGroupLike = Union{_PARAMETER_GROUP_TYPES...}

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
