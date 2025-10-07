module EnsembleSamples

using ..QDistributions: QDistribution
using ..QInterpolations: Interpolator

abstract type AbstractEnsembleSample end

"""
    DiscreteSamples(method, group_indices, group_symbols, group_names, distributions, samples)

Container describing discrete ensemble samples. `samples` is an `n×m` matrix whose
rows correspond to joint sample tuples and columns align with `group_indices`.
`method` identifies the sampling strategy (e.g. `:random`, `:density`).
"""
struct DiscreteSamples <: AbstractEnsembleSample
    method::Symbol
    group_indices::Vector{Int}
    group_symbols::Vector{Symbol}
    group_names::Vector{String}
    distributions::Vector{QDistribution}
    samples::Matrix{Float64}
end

"""
    ContinuousSamples(method, group_indices, group_symbols, group_names, distributions, samples, interpolator)

Container describing continuum-oriented ensemble samples along an interpolation grid.
`samples` stores the node coordinates; `interpolator` retains the `Interpolator`
constructed from those nodes.
"""
struct ContinuousSamples <: AbstractEnsembleSample
    method::Symbol
    group_indices::Vector{Int}
    group_symbols::Vector{Symbol}
    group_names::Vector{String}
    distributions::Vector{QDistribution}
    samples::Matrix{Float64}
    interpolator::Interpolator
end

end # module EnsembleSamples
