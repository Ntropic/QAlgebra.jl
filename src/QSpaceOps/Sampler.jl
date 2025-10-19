module Sampler

using Random
using Base: CartesianIndices
using ..QAlgebra: SAMPLE_ATOL, SAMPLE_RTOL

include("../SampleHelpers/Interpolators.jl")
include("../SampleHelpers/Distributions.jl")
include("../SampleHelpers/Integrators.jl")
include("../SampleHelpers/SampleHelpers.jl")

module EnsembleSamples

export AbstractEnsembleSample, DiscreteSamples, ContinuousSamples

abstract type AbstractEnsembleSample end

"""
    DiscreteSamples(method, group_indices, group_symbols, group_names, distributions, samples)

Container describing discrete ensemble samples produced from inverse-CDF sampling.
`samples` is an `n×m` matrix whose rows enumerate joint samples and whose columns align
with `group_indices`/`group_symbols`.
"""
struct DiscreteSamples <: AbstractEnsembleSample
    method::Symbol
    group_indices::Vector{Int}
    group_symbols::Vector{Symbol}
    group_names::Vector{String}
    distributions::Vector{parentmodule(@__MODULE__).QDistribution}
    samples::Matrix{Float64}
end

"""
    ContinuousSamples(method, group_indices, group_symbols, group_names, distributions, samples, interpolator)

Container for grid-based samples used in continuous ensemble integrations.
`samples` stores tensor-product node coordinates; `interpolator` retains the
`QInterpolator` built from those nodes.
"""
struct ContinuousSamples <: AbstractEnsembleSample
    method::Symbol
    group_indices::Vector{Int}
    group_symbols::Vector{Symbol}
    group_names::Vector{String}
    distributions::Vector{parentmodule(@__MODULE__).QDistribution}
    samples::Matrix{Float64}
    interpolator::parentmodule(@__MODULE__).QInterpolator
end

end # module EnsembleSamples

using .EnsembleSamples: AbstractEnsembleSample, DiscreteSamples, ContinuousSamples

function nodes(sample::Union{EnsembleSamples.DiscreteSamples, EnsembleSamples.ContinuousSamples}, separate::Bool=false)
    if !separate
        return sample.samples
    else
        return [sample.samples[:, i] for i in 1:size(sample.samples, 2)]
    end
end

export build_discrete_samples, build_continuous_samples
export EnsembleSamples, AbstractEnsembleSample, DiscreteSamples, ContinuousSamples

# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

function _ensure_consistent_sample_count(dists::AbstractVector{<:QDistribution})
    isempty(dists) && return 0
    counts = unique(dist.num_samples for dist in dists)
    length(counts) == 1 || error("Distributions combined within an ensemble must share the same num_samples (got $(collect(counts))).")
    return counts[1]
end

@inline function _clamp_to_support!(vals::AbstractVector{<:Real}, dist::QDistribution)
    a = dist.minimum
    b = dist.maximum
    @inbounds for i in eachindex(vals)
        vi = vals[i]
        if vi < a
            vals[i] = a
        elseif vi > b
            vals[i] = b
        end
    end
    return vals
end

# ------------------------------------------------------------
# Sample matrix builders
# ------------------------------------------------------------

function _build_discrete_matrix(method::Symbol,
                                dists::AbstractVector{<:QDistribution};
                                rng::AbstractRNG=Random.default_rng(),
                                num_nodes::Int=25,
                                atol::Float64=SAMPLE_ATOL,
                                rtol::Float64=SAMPLE_RTOL,
                                max_iter::Int=128)
    isempty(dists) && return zeros(Float64, 0, 0)
    n = _ensure_consistent_sample_count(dists)
    cols = length(dists)
    samples = Matrix{Float64}(undef, n, cols)

    for (j, dist) in enumerate(dists)
        cdf_inter = pdf2cdf(dist;
                            num_nodes=num_nodes,
                            atol=atol,
                            rtol=rtol)
        inv_inter = cdf2inverse(cdf_inter;
                                num_nodes=num_nodes,
                                atol=atol,
                                rtol=rtol,
                                max_iter=max_iter)

        values = if method == :random
            us = sort!(rand(rng, n))
            [Float64(inv_inter(u)) for u in us]
        elseif method == :density
            us = ((1:n) .- 0.5) ./ n
            [Float64(inv_inter(u)) for u in us]
        else
            error("Unsupported discrete sampling method $method. Choose :random or :density.")
        end

        _clamp_to_support!(values, dist)
        samples[:, j] = values
    end
    return samples
end

function _build_continuous_matrix(method::Symbol,
                                  dists::AbstractVector{<:QDistribution};
                                  endpoints::Bool=true)
    isempty(dists) && return (zeros(Float64, 0, 0), nothing)
    params = [(max(dist.num_samples, 2), dist.minimum, dist.maximum) for dist in dists]
    pdfs = [x -> pdf(dist, x) for dist in dists]
    bounds, nodes = build_interpolation_nodes(method, params; pdfs=pdfs, endpoints=endpoints)
    dims = length(nodes)
    total = prod(length.(nodes))
    samples = Matrix{Float64}(undef, total, dims)
    idx = 1
    sz = ntuple(i -> length(nodes[i]), dims)
    for CI in CartesianIndices(sz)
        @inbounds for d in 1:dims
            samples[idx, d] = nodes[d][CI[d]]
        end
        idx += 1
    end
    inter = QInterpolator(nodes, bounds; method=method, endpoints=endpoints)
    return samples, inter
end

# ------------------------------------------------------------
# Public constructors
# ------------------------------------------------------------

"""
    build_discrete_samples(ensemble, group_indices, group_symbols, group_names, dists; kwargs...)

Assemble an `EnsembleSamples.DiscreteSamples` object for the provided ensemble
groups using inverse-CDF sampling of each `QDistribution`.

Arguments:
- `ensemble`: object describing the parent ensemble; its `sample_method` supplies the default sampling mode.
- `group_indices`, `group_symbols`, `group_names`: metadata aligning ensemble groups with matrix columns.
- `dists::Vector{<:QDistribution}`: per-column marginal distributions supplying `num_samples`, `minimum`, and `maximum`.
Keyword arguments:
- `method::Symbol = ensemble.sample_method`: `:random` sorts uniform draws; `:density` uses equispaced quantiles (resolves to `:random` when `:default`).
- `rng::AbstractRNG = Random.default_rng()`: source of randomness when `method == :random`.
- `num_nodes::Int = ensemble.sample_num_nodes`: number of Chebyshev nodes used when approximating the CDF and its inverse.
- `atol::Float64 = ensemble.sample_atol` / `rtol::Float64 = ensemble.sample_rtol`: tolerances for `quadgk` integrations and inverse refinement.
- `max_iter::Int = ensemble.sample_max_iter`: caps the inverse-CDF refinement iterations (forwarded to `cdf2inverse`).
"""
function build_discrete_samples(ensemble,
                                group_indices::Vector{Int},
                                group_symbols::Vector{Symbol},
                                group_names::Vector{String},
                                dists::AbstractVector{<:QDistribution};
                                method::Symbol=ensemble.sample_method,
                                rng::AbstractRNG=Random.default_rng(),
                                num_nodes::Int=ensemble.sample_num_nodes,
                                atol::Float64=ensemble.sample_atol,
                                rtol::Float64=ensemble.sample_rtol,
                                max_iter::Int=ensemble.sample_max_iter)
    requested_method = method === :default ? :random : method
    m = Symbol(lowercase(String(requested_method)))
    samples = _build_discrete_matrix(m, dists;
                                     rng=rng,
                                     num_nodes=num_nodes,
                                     atol=atol,
                                     rtol=rtol,
                                     max_iter=max_iter)
    return DiscreteSamples(m,
                           copy(group_indices),
                           copy(group_symbols),
                           copy(group_names),
                           copy(dists),
                           samples)
end

"""
    build_continuous_samples(ensemble, group_indices, group_symbols, group_names, dists; kwargs...)

Construct an `EnsembleSamples.ContinuousSamples` object by forming a tensor grid
over the per-group distributions and storing the associated `QInterpolator`.

Arguments:
- `ensemble`: parent ensemble definition; `sample_method` sets the default interpolation scheme.
- `group_indices`, `group_symbols`, `group_names`: metadata for each ensemble group.
- `dists::Vector{<:QDistribution}`: distributions providing bounds and `num_samples` used for node counts.

Keyword arguments:
- `method::Symbol = ensemble.sample_method`: interpolation node selector (`:chebyshev`, `:uniform`, `:leja`, ...) resolving to `:chebychev` when `:default`.
- `endpoints::Bool = true`: ensure interval endpoints participate in the grid where the method allows it.
"""
function build_continuous_samples(ensemble,
                                  group_indices::Vector{Int},
                                  group_symbols::Vector{Symbol},
                                  group_names::Vector{String},
                                  dists::AbstractVector{<:QDistribution};
                                  method::Symbol=ensemble.sample_method,
                                  endpoints::Bool=true)
    requested_method = method === :default ? :chebychev : method
    m = Symbol(lowercase(String(requested_method)))
    samples, interpolator = _build_continuous_matrix(m, dists; endpoints=endpoints)
    return ContinuousSamples(m,
                             copy(group_indices),
                             copy(group_symbols),
                             copy(group_names),
                             copy(dists),
                             samples,
                             interpolator)
end

end # module Sampler
