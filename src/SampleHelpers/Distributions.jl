export QDistribution, pdf, QNormal, QUniform, QEnsembleFunction

using ..StringUtils: underscore_separate

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

    inter = QInterpolator(:uniform, [(max(num_samples, 2), minv, maxv)])
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
    QEnsembleFunction(name, group_indexes, argument_signatures, func)

Metadata wrapper for ensemble parameter functions that depend on other parameters
(or time) when sampling across ensembles.  Besides the callable `func`, the
struct stores how each argument relates to the parent group's abstract indexes
and which parameter group supplies the samples.

# Arguments
- `name`: base name of the parameter group.
- `group_indexes`: abstract index labels of the parent group (e.g. `["i","j"]`).
- `argument_signatures`: ordered argument signatures (e.g. `["t", "delta_i"]`).
- `func`: callable evaluated with arguments matching `argument_signatures`.
"""
struct QEnsembleFunction
    name::String
    argument_symbols::Vector{Symbol}
    func::Function
    argument_group_names::Vector{String}
    argument_group_indices::Vector{Int}
    argument_self_index_positions::Vector{Vector{Int}}
    function QEnsembleFunction(name::String, argument_symbols::Vector{Symbol}, func::Function, argument_group_names::Vector{String}, argument_group_indices::Vector{Int}, argument_self_index_positions::Vector{Vector{Int}})
        isempty(argument_symbols) && error("QEnsembleFunction requires at least one argument symbol; include e.g. t for time or parameter names.")
        return new(name, argument_symbols, func, argument_group_names, argument_group_indices, argument_self_index_positions)
    end
end

function QEnsembleFunction(name::String, group_indexes::Vector{String}, argument_signatures::Vector{String}, func::Function)
    isempty(argument_signatures) && error("QEnsembleFunction requires at least one argument; include e.g. t for time or parameter names.")
    arg_symbols = Symbol.(argument_signatures)
    arg_group_names = Vector{String}(undef, length(argument_signatures))
    self_positions = Vector{Vector{Int}}(undef, length(argument_signatures))
    for (idx, spec) in enumerate(argument_signatures)
        if spec == "t"
            arg_group_names[idx] = "t"
            self_positions[idx] = Int[]
            continue
        end
        base, tokens = underscore_separate(spec)
        arg_group_names[idx] = base
        positions = Vector{Int}(undef, length(tokens))
        for (pos_idx, tok) in enumerate(tokens)
            pos = findfirst(==(tok), group_indexes)
            pos === nothing &&
                error("Argument \"$spec\" references index \"$tok\" which is not defined on ensemble function group \"$name\".")
            positions[pos_idx] = pos
        end
        self_positions[idx] = positions
    end
    return QEnsembleFunction(name, arg_symbols, func, arg_group_names, zeros(Int, length(argument_signatures)), self_positions)
end

function QEnsembleFunction(name::String, argument_symbols::Vector{Symbol}, func::Function)
    return QEnsembleFunction(name, String[], String.(argument_symbols), func)
end

function _pdfs_and_scaling(inter::QInterpolator,
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

function integrate_node_funs(inter::QInterpolator,
                             distributions::Vector{QDistribution};
                             f::Union{Nothing,Function}=nothing,
                             constant::Float64=1.0,
                             atol::Float64=1e-9,
                             rtol::Float64=1e-7)
    pdfs, scale = _pdfs_and_scaling(inter, distributions)
    total_constant = constant * scale
    return integrate_node_funs(inter, pdfs;
                               f=f,
                               constant=total_constant,
                               atol=atol,
                               rtol=rtol)
end

function Integrator(inter::QInterpolator, distributions::Vector{QDistribution}; f::Union{Nothing,Function}=nothing, atol::Float64=1e-9, rtol::Float64=1e-7)
    pdfs, scale = _pdfs_and_scaling(inter, distributions)
    return QIntegrator(inter, pdfs;  f=f, constant=scale, atol=atol, rtol=rtol)
end
