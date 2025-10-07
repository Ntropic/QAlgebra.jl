module Sampler

export build_discrete_samples, build_continuous_samples, attach_samples!

using Random
using QuadGK
using Base: CartesianIndices
using ...QDistributions: QDistribution, pdf
using ...CFunctions: ParameterValues
using ...EnsembleSamples
using ...EnsembleSamples: AbstractEnsembleSample, DiscreteSamples, ContinuousSamples
using ...QInterpolations: Interpolator, build_interpolation_nodes

# ------------------------------------------------------------
# Quantile utilities
# ------------------------------------------------------------

function _cdf(dist::QDistribution, x::Float64)
    if x <= dist.minimum
        return 0.0
    elseif x >= dist.maximum
        return 1.0
    end
    val, _ = QuadGK.quadgk(y -> pdf(dist, y), dist.minimum, x; atol=1e-10, rtol=1e-8)
    return min(max(val, 0.0), 1.0)
end

function _quantile(dist::QDistribution, u::Float64; tol::Float64=1e-9, max_iter::Int=60)
    lower = dist.minimum
    upper = dist.maximum
    while upper - lower > tol && max_iter > 0
        mid = 0.5 * (lower + upper)
        cdf_mid = _cdf(dist, mid)
        if cdf_mid > u
            upper = mid
        else
            lower = mid
        end
        max_iter -= 1
    end
    return 0.5 * (lower + upper)
end

function _quantile_vector(dist::QDistribution, us::AbstractVector{<:Real})
    [ _quantile(dist, Float64(u)) for u in us ]
end

function _random_quantiles(dist::QDistribution, n::Int, rng::AbstractRNG)
    us = rand(rng, n)
    sort!(us)
    return _quantile_vector(dist, us)
end

function _density_quantiles(dist::QDistribution, n::Int)
    if n == 1
        return [_quantile(dist, 0.5)]
    end
    us = ((1:n) .- 0.5) ./ n
    return _quantile_vector(dist, us)
end

function _ensure_consistent_sample_count(dists::AbstractVector{<:QDistribution})
    counts = unique(dist.num_samples for dist in dists)
    length(counts) == 1 ||
        error("Distributions combined within an ensemble must share the same num_samples (got $(collect(counts))).")
    return counts[1]
end

# ------------------------------------------------------------
# Sample matrix builders
# ------------------------------------------------------------

function _build_discrete_matrix(method::Symbol, dists::AbstractVector{<:QDistribution}; rng::AbstractRNG=Random.default_rng())
    isempty(dists) && return zeros(Float64, 0, 0)
    n = _ensure_consistent_sample_count(dists)
    cols = length(dists)
    samples = Matrix{Float64}(undef, n, cols)
    for (j, dist) in enumerate(dists)
        values = if method == :random
            _random_quantiles(dist, n, rng)
        elseif method == :density
            _density_quantiles(dist, n)
        else
            error("Unsupported discrete sampling method $method.")
        end
        samples[:, j] = values
    end
    return samples
end

function _build_continuous_matrix(method::Symbol, dists::AbstractVector{<:QDistribution}; endpoints::Bool=true, M::Int=0)
    isempty(dists) && return (zeros(Float64, 0, 0), nothing)
    params = [(max(dist.num_samples, 2), dist.minimum, dist.maximum) for dist in dists]
    pdfs = [x -> pdf(dist, x) for dist in dists]
    bounds, nodes = build_interpolation_nodes(method, params; pdfs=pdfs, endpoints=endpoints, M=M)
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
    inter = Interpolator(nodes, bounds; method=method, endpoints=endpoints)
    return samples, inter
end

# ------------------------------------------------------------
# Public constructors
# ------------------------------------------------------------

function build_discrete_samples(ensemble, group_indices::Vector{Int}, group_symbols::Vector{Symbol},
                                group_names::Vector{String}, dists::AbstractVector{<:QDistribution};
                                method::Symbol=ensemble.discrete_method,
                                rng::AbstractRNG=Random.default_rng())
    m = Symbol(lowercase(String(method)))
    samples = _build_discrete_matrix(m, dists; rng=rng)
    return DiscreteSamples(m, copy(group_indices), copy(group_symbols), copy(group_names),
                           copy(dists), samples)
end

function build_continuous_samples(ensemble, group_indices::Vector{Int}, group_symbols::Vector{Symbol},
                                  group_names::Vector{String}, dists::AbstractVector{<:QDistribution};
                                  method::Symbol=ensemble.continuous_method,
                                  endpoints::Bool=true, M::Int=0)
    m = Symbol(lowercase(String(method)))
    samples, interpolator = _build_continuous_matrix(m, dists; endpoints=endpoints, M=M)
    return ContinuousSamples(m, copy(group_indices), copy(group_symbols), copy(group_names),
                             copy(dists), samples, interpolator)
end

# ------------------------------------------------------------
# ParameterValues integration helpers
# ------------------------------------------------------------

function _assign_group_storage!(storage, values)
    if storage === nothing || storage isa Number
        return values
    elseif storage isa AbstractArray
        storage .= Ref(values)
        return storage
    else
        error("Unsupported storage type $(typeof(storage)) for ensemble sample assignment.")
    end
end

function attach_samples!(pv::ParameterValues, sample::AbstractEnsembleSample)
    for (col, group_idx) in enumerate(sample.group_indices)
        values = sample.samples[:, col]
        storage = pv.group_values[group_idx]
        pv.group_values[group_idx] = _assign_group_storage!(storage, values)
        pv.group_initialized[group_idx] = true
        pv.ensemble_group_samples[group_idx] = sample
    end
    return sample
end

end # module Sampler
