module ParameterGroups

export ParameterGroupKind, ParameterGroup, ParameterGroupScalar, ParameterGroupTimeFunction,
       ParameterGroupDistribution, ParameterGroupEnsembleFunction, ParameterGroupPayload

using ..Sampler: QDistribution, QEnsembleFunction

"""
    ParameterGroupKind

Enumerates the supported parameter-group forms used across the `QSpace`
pipeline.  The kind determines how a group's payload is validated and how its
values are materialised:

  * `ParameterGroupScalar` – literal value or array shared by all indices.
  * `ParameterGroupTimeFunction` – scalar function of time.
  * `ParameterGroupDistribution` – distribution backing ensemble sampling.
  * `ParameterGroupEnsembleFunction` – function that maps distribution samples
    (and possibly time) to concrete values.
"""
@enum ParameterGroupKind::UInt8 begin
    ParameterGroupScalar
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

end # module ParameterGroups
