module ParameterGroups

export ParameterGroupKind, ParameterGroup, ParameterGroupScalar, ParameterGroupTimeScalar, ParameterGroupTimeFunction,
       ParameterGroupDistribution, ParameterGroupEnsembleFunction, ParameterGroupPayload, WhereWhichParamGroup

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
  * `ParameterGroupEnsembleFunction` – function that maps distribution samples
    (and possibly time) to concrete values.
"""
@enum ParameterGroupKind::UInt8 begin
    ParameterGroupScalar
    ParameterGroupTimeScalar
    ParameterGroupTimeFunction
    ParameterGroupDistribution
    ParameterGroupEnsembleFunction
end

const ParameterGroupPayload = Union{Nothing, Function, QDistribution, QEnsembleFunction, Number, AbstractArray}

"""
    ParameterGroup

Shared metadata for a parameter group.  Instances are created during
`ParameterDefinitions2Parameters` and then referenced by `ParameterInfo`,
`ParameterValues`, and each `Ensemble`.  The mutable `payload` field stores the
user-provided definition (distribution, function, literal value, …) and may be
updated after `QSpace` construction via `resolve_param!`.

Fields capture the group's declarative signature (`name`, `indexes`,
`function_args`), dependency tracking (`dependency_names`/`dependency_indices`),
ensemble affiliation (`ensemble_outer_indices`, `ensemble_presence`), the list
of concrete parameters created for the group (`parameter_indices`), and derived
shape information (`time_count`, `index_sizes`).
"""
mutable struct ParameterGroup
    name::Symbol
    display_signature::String
    kind::ParameterGroupKind
    of_t::Bool
    indexes::Vector{String}
    function_args::Vector{String}
    dependency_names::Vector{String}
    dependency_indices::Vector{Int}
    ensemble_outer_indices::Vector{Int}
    index_outer_subspaces::Vector{Int}
    ensemble_presence::BitVector
    parameter_indices::Vector{Int}
    time_count::Int
    index_sizes::Vector{Int}
    sample_sizes::Vector{Int}
    is_time_group::Bool
    payload::ParameterGroupPayload
end

struct WhereWhichParamGroup
    scalar_groups::Vector{Int}
    time_function_groups::Vector{Int}
    distribution_groups::Vector{Int}
    ensemble_function_groups::Vector{Int}
    time_group::Int
    first_qensemble_group::Int
end

function WhereWhichParamGroup(groups::Vector{ParameterGroup})
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
        elseif kind == ParameterGroupEnsembleFunction
            push!(ensemble_function_groups, idx)
            first_qensemble_group == 0 && (first_qensemble_group = idx)
        end
    end
    first_qensemble_group == 0 && (first_qensemble_group = length(groups) + 1)
    return WhereWhichParamGroup(scalar_groups, time_function_groups, distribution_groups,
                                ensemble_function_groups, time_group, first_qensemble_group)
end

end # module ParameterGroups
