module QDistributions

export QDistribution, pdf, QNormal, QUniform, QEnsembleFunction

using ..QInterpolations: Interpolator, normalization_constant
import ..QInterpolations: integrate_node_funs, Integrator

const PDF_BOUND_ATOL = 1e-10
const PDF_BOUND_RTOL = 1e-8

"""
    QDistribution(minimum, maximum, pdf, num_samples)

Probability distribution describing an ensemble parameter.

# Arguments
- `minimum`, `maximum`: bounds of the support interval (must satisfy `minimum < maximum`).
- `pdf`: non-negative weight function defined over the interval.
- `num_samples`: number of samples along the associated ensemble direction.
- The provided density is always normalized via `normalization_constant` from the interpolation utilities.
"""
struct QDistribution{F<:Function}
    minimum::Float64
    maximum::Float64
    pdf::F
    num_samples::Int
    normalization_constant::Float64
end

"""
    (dist::QDistribution)(x)

Evaluate the normalized probability density associated with `dist` at `x`.
Outside the support interval the density is zero.
"""
function (dist::QDistribution)(x::Real)
    return pdf(dist, x)
end

"""
    pdf(dist::QDistribution, x)

Return the normalized density value of `dist` at `x`. Values outside the
support evaluate to `0.0`.
"""
function pdf(dist::QDistribution, x::Real)
    xv = Float64(x)
    if xv < dist.minimum || xv > dist.maximum
        return 0.0
    end
    return dist.normalization_constant * max(dist.pdf(xv), 0.0)
end

function QDistribution(minimum::Real, maximum::Real, pdf::Function, num_samples::Int)
    minv = Float64(minimum)
    maxv = Float64(maximum)
    minv < maxv || error("QDistribution bounds must satisfy minimum < maximum (got $minimum, $maximum).")
    num_samples > 0 || error("num_samples must be a positive integer (got $num_samples).")

    base_pdf = x -> Float64(pdf(Float64(x)))

    inter = Interpolator(:uniform, [(max(num_samples, 2), minv, maxv)])
    norm = normalization_constant([base_pdf], inter)

    return QDistribution{typeof(base_pdf)}(minv, maxv, base_pdf, num_samples, norm)
end

"""
    QUniform(minimum, maximum, num_samples)

Construct a uniform distribution over `[minimum, maximum]` with the given
number of samples.
"""
function QUniform(minimum::Real, maximum::Real, num_samples::Int)
    return QDistribution(minimum, maximum, x -> 1.0, num_samples)
end

"""
    QNormal(mean, std, how_many_stds, num_samples)

Construct a truncated Gaussian centred at `mean` with standard deviation `std`.
The support interval spans `mean ± how_many_stds * std`. The resulting
distribution is normalized over that interval.
"""
function QNormal(mean::Real, std::Real, how_many_stds::Real, num_samples::Int)
    std > 0 || error("QNormal requires std > 0 (got $std).")
    how_many_stds > 0 || error("QNormal requires how_many_stds > 0 (got $how_many_stds).")
    num_samples > 0 || error("num_samples must be a positive integer (got $num_samples).")

    μ = Float64(mean)
    σ = Float64(std)
    span = Float64(how_many_stds) * σ
    minv = μ - span
    maxv = μ + span
    pdf_fun = x -> exp(-0.5 * ((x - μ) / σ)^2)
    return QDistribution(minv, maxv, pdf_fun, num_samples)
end

"""
    QEnsembleFunction(name, argument_symbols, func)

Metadata wrapper for ensemble parameter functions that depend on other parameters
(or time) when sampling across ensembles.

- `name`: base name of the parameter group.
- `argument_symbols`: ordered symbols describing function arguments (e.g. `[:t, :alpha, :beta]`).
- `func`: callable evaluated with arguments matching `argument_symbols`.
"""
struct QEnsembleFunction
    name::String
    argument_symbols::Vector{Symbol}
    func::Function
    function QEnsembleFunction(name::String, argument_symbols::Vector{Symbol}, func::Function)
        isempty(argument_symbols) && error("QEnsembleFunction requires at least one argument symbol; include e.g. t for time or parameter names.")
        return new(name, argument_symbols, func)
    end
end

function _pdfs_and_scaling(inter::Interpolator,
                           distributions::Vector{QDistribution})
    d = inter.dims
    length(distributions) == d ||
        error("Expected one distribution per interpolation dimension (got $(length(distributions)), expected $d).")

    pdfs = Vector{Function}(undef, d)
    scale = 1.0

    @inbounds for i in 1:d
        dist = distributions[i]
        pdfs[i] = dist.pdf
        scale *= dist.normalization_constant

        ai, bi = inter.bounds[i]
        minv = dist.minimum
        maxv = dist.maximum

        if !(isapprox(ai, minv; atol=PDF_BOUND_ATOL, rtol=PDF_BOUND_RTOL) &&
             isapprox(bi, maxv; atol=PDF_BOUND_ATOL, rtol=PDF_BOUND_RTOL))
            error("Distribution bounds ($(minv), $(maxv)) do not match interpolation bounds ($(ai), $(bi)) along dimension $i.")
        end
    end

    return pdfs, scale
end

function integrate_node_funs(inter::Interpolator,
                             distributions::Vector{QDistribution};
                             f::Union{Nothing,Function}=nothing,
                             constant::Float64=1.0,
                             atol::Float64=1e-9,
                             rtol::Float64=1e-7)
    pdfs, scale = _pdfs_and_scaling(inter, distributions)
    total_constant = constant * scale
    return QInterpolations.integrate_node_funs(inter, pdfs;
                                               f=f,
                                               constant=total_constant,
                                               atol=atol,
                                               rtol=rtol)
end

function Integrator(inter::Interpolator,
                    distributions::Vector{QDistribution};
                    f::Union{Nothing,Function}=nothing,
                    atol::Float64=1e-9,
                    rtol::Float64=1e-7)
    pdfs, scale = _pdfs_and_scaling(inter, distributions)
    return QInterpolations.Integrator(inter, pdfs;
                                      f=f,
                                      constant=scale,
                                      atol=atol,
                                      rtol=rtol)
end

end # module QDistributions
