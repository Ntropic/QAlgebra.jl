using QuadGK
using ..QAlgebra: QUADGK_ATOL, QUADGK_RTOL

export QInterpolator, build_interpolation_nodes, eval_interpolation, nodes, basis_values, basis_values!
export QIntegrator, integrate_node_funs, normalization_constant, eval_integration
export QDistribution, pdf, QNormal, QUniform, QEnsembleFunction
export pdf2cdf, cdf2inverse

function _pdf2cdf_core(pdf_fn::Function, a::Float64, b::Float64;
                       num_nodes::Int=25,
                       atol::Float64=QUADGK_ATOL,
                       rtol::Float64=QUADGK_RTOL)
    n = max(num_nodes, 2)
    params = [(n, a, b)]
    bounds, nodes_per_dim = build_interpolation_nodes(:chebychev, params;
                                                      endpoints=true)
    nodes = copy(nodes_per_dim[1])
    sort!(nodes)

    cdf_vals = Vector{Float64}(undef, length(nodes))
    accum = 0.0
    prev = a
    for (i, xnode) in pairs(nodes)
        integral = xnode > prev ? QuadGK.quadgk(pdf_fn, prev, xnode;
                                                atol=atol,
                                                rtol=rtol)[1] : 0.0
        accum += integral
        cdf_vals[i] = clamp(accum, 0.0, 1.0)
        prev = xnode
    end

    if prev < b
        tail, _ = QuadGK.quadgk(pdf_fn, prev, b; atol=atol, rtol=rtol)
        accum += tail
    end

    accum > 0 || error("Computed CDF integral is non-positive.")
    cdf_vals ./= accum
    cdf_vals[1] = 0.0
    cdf_vals[end] = 1.0

    return QInterpolator([nodes], bounds; method=:chebychev, endpoints=true, values=cdf_vals)
end

"""
    pdf2cdf(dist::QDistribution; kwargs...)

Approximate the cumulative distribution function for `dist` by sampling its pdf on a
one-dimensional interpolation grid. Returns a `QInterpolator` whose stored values are
the normalised CDF.

Keyword arguments:
- `num_nodes::Int = 25`: number of Chebyshev interpolation nodes spanning the support.
- `atol::Float64 = QUADGK_ATOL` / `rtol::Float64 = QUADGK_RTOL`: absolute/relative tolerances for `QuadGK`.
"""
function pdf2cdf(dist::QDistribution;
                 num_nodes::Int=25,
                 atol::Float64=QUADGK_ATOL,
                 rtol::Float64=QUADGK_RTOL)
    return _pdf2cdf_core(x -> pdf(dist, x), dist.minimum, dist.maximum;
                         num_nodes=num_nodes,
                         atol=atol,
                         rtol=rtol)
end

"""
    pdf2cdf(pdf_inter::QInterpolator; kwargs...)

Build a CDF interpolator directly from a stored pdf interpolator (`dims == 1`). The
input interpolator must carry stored values; they are integrated and normalised before
returning the new `QInterpolator`.
"""
function pdf2cdf(pdf_inter::QInterpolator;
                 num_nodes::Int=25,
                 atol::Float64=QUADGK_ATOL,
                 rtol::Float64=QUADGK_RTOL)
    pdf_inter.dims == 1 || error("pdf2cdf currently supports only 1D pdf interpolators.")
    values = pdf_inter.default_values
    values === nothing && error("pdf2cdf(pdf_inter) requires the interpolator to store pdf samples.")
    a, b = pdf_inter.bounds[1]
    return _pdf2cdf_core(x -> pdf_inter(x), a, b;
                         num_nodes=num_nodes,
                         atol=atol,
                         rtol=rtol)
end

"""
    cdf2inverse(cdf_inter::QInterpolator; kwargs...)

Construct an interpolator for the inverse CDF using a monotone one-dimensional CDF.
Sampling against the returned object enables inverse-transform draws and deterministic
density nodes.

Keyword arguments:
- `num_nodes::Int = 25`: number of Chebyshev probability nodes spanning `[0, 1]`.
- `atol::Float64 = QUADGK_ATOL` / `rtol::Float64 = QUADGK_RTOL`: tolerances for refinement of the inverse search.
- `max_iter::Int = 128`: maximum bisection iterations used per probability value.
"""
function cdf2inverse(cdf_inter::QInterpolator;
                     num_nodes::Int=25,
                     atol::Float64=QUADGK_ATOL,
                     rtol::Float64=QUADGK_RTOL,
                     max_iter::Int=128)
    cdf_inter.dims == 1 || error("cdf2inverse currently supports only 1D CDFs.")
    values = cdf_inter.default_values
    values === nothing && error("cdf2inverse requires the CDF interpolator to carry stored values.")
    cdf_vals = vec(copy(values))
    nodes = copy(cdf_inter.nodes[1])
    ord = sortperm(nodes)
    nodes = nodes[ord]
    cdf_vals = cdf_vals[ord]

    n = max(num_nodes, 2)
    bounds_prob, prob_nodes_dim = build_interpolation_nodes(:chebychev, [(n, 0.0, 1.0)];
                                                            endpoints=true)
    prob_nodes = copy(prob_nodes_dim[1])
    sort!(prob_nodes)

    inv_values = Vector{Float64}(undef, length(prob_nodes))
    xmin, xmax = nodes[1], nodes[end]
    tol_from_p(p) = max(atol, rtol * max(abs(p), 1.0))

    for (i, p) in pairs(prob_nodes)
        if p <= cdf_vals[1] + tol_from_p(p)
            inv_values[i] = xmin
            continue
        elseif p >= cdf_vals[end] - tol_from_p(p)
            inv_values[i] = xmax
            continue
        end
        idx = searchsortedlast(cdf_vals, p)
        idx = clamp(idx, 1, length(cdf_vals) - 1)
        lo = nodes[idx]
        hi = nodes[idx + 1]
        f_lo = cdf_vals[idx]
        f_hi = cdf_vals[idx + 1]

        for _ in 1:max_iter
            mid = 0.5 * (lo + hi)
            f_mid = cdf_inter(mid)
            if abs(f_mid - p) <= tol_from_p(p) || abs(hi - lo) <= max(atol, rtol * max(abs(hi), abs(lo), 1.0))
                lo = hi = mid
                f_lo = f_hi = f_mid
                break
            elseif f_mid < p
                lo = mid
                f_lo = f_mid
            else
                hi = mid
                f_hi = f_mid
            end
        end
        inv_values[i] = 0.5 * (lo + hi)
    end

    return QInterpolator([prob_nodes], bounds_prob; method=:chebychev, endpoints=true, values=inv_values)
end

"""
    cdf2inverse(nodes, values; kwargs...)

Convenience wrapper accepting monotone CDF samples directly. Wraps the input into a
`QInterpolator` before calling [`cdf2inverse(::QInterpolator)`](@ref).
"""
function cdf2inverse(nodes::AbstractVector{<:Real},
                     values::AbstractVector{<:Real}; kwargs...)
    length(nodes) == length(values) ||
        error("nodes and values must have the same length.")
    n = Float64.(nodes)
    v = Float64.(values)
    ord = sortperm(n)
    n = n[ord]
    v = v[ord]
    bounds = [(n[1], n[end])]
    cdf_inter = QInterpolator([n], bounds; method=:custom, endpoints=true, values=v)
    return cdf2inverse(cdf_inter; kwargs...)
end

cdf2inverse(cdf_data::NamedTuple; kwargs...) =
    cdf2inverse(cdf_data.nodes, cdf_data.values; kwargs...)
