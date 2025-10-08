module QIntegrators

using ..QInterpolators: QInterpolator, basis_values!

import QuadGK

export normalization_constant, integrate_node_funs, QIntegrator, eval_integration

const _PDF_BOUND_ATOL = 1e-10
const _PDF_BOUND_RTOL = 1e-8

@inline _unity_function(::Any) = 1.0

function normalization_constant(pdfs::AbstractVector{<:Function}, inter::QInterpolator;
                                atol::Float64=1e-9, rtol::Float64=1e-7)::Float64
    d = inter.dims
    length(pdfs) == d || error("normalization_constant requires one pdf per interpolation dimension (got $(length(pdfs)), expected $d).")

    Z = 1.0
    for i in 1:d
        ai, bi = inter.bounds[i]
        ρi = pdfs[i]
        Zi, _ = QuadGK.quadgk(x -> max(Float64(ρi(x)), 0.0), ai, bi; atol=atol, rtol=rtol, norm=abs)
        Zi <= 0 && error("Normalization failed: integral of pdf[$i] over [$ai, $bi] is non-positive.")
        Z *= Zi
    end
    return 1.0 / Z
end

function _integrate_node_funs(inter::QInterpolator,
                              pdfs::Vector{<:Function},
                              f::Function;
                              constant::Float64=1.0,
                              atol::Float64=1e-9,
                              rtol::Float64=1e-7)
    d = inter.dims
    length(pdfs) == d || error("Expected one pdf per interpolation dimension (got $(length(pdfs)), expected $d).")

    sz  = ntuple(i -> length(inter.nodes[i]), d)
    Lx  = Array{Float64}(undef, sz)
    xpt = zeros(Float64, d)
    arr_norm = x -> sqrt(sum(abs2, x))

    function _integrate_dim(i::Int, ρaccum::Float64)
        ai, bi = inter.bounds[i]
        ρi = pdfs[i]

        if i < d
            g = function (xi::Float64)
                xpt[i] = xi
                ρnew = ρaccum * max(Float64(ρi(xi)), 0.0)
                return _integrate_dim(i+1, ρnew)
            end
            val, _ = QuadGK.quadgk(g, ai, bi; atol=atol, rtol=rtol, norm=arr_norm)
            return val
        else
            glast = function (xd::Float64)
                xpt[d] = xd
                basis_values!(Lx, inter, xpt)
                ρx = ρaccum * max(Float64(ρi(xd)), 0.0)
                fx = f(xpt)
                return Lx .* (ρx * fx)
            end
            val, _ = QuadGK.quadgk(glast, ai, bi; atol=atol, rtol=rtol, norm=arr_norm)
            return val
        end
    end

    W = _integrate_dim(1, 1.0)
    W .*= constant
    return W
end

"""
    integrate_node_funs(inter::QInterpolator,
                        pdfs;
                        f::Union{Nothing,Function}=nothing,
                        constant::Float64=1.0,
                        atol::Float64=1e-9,
                        rtol::Float64=1e-7)

Compute `W[j...] = ∫ ρ(x) f(x) L_j(x) dx` over the interpolation domain.
`pdfs` must be a `Vector{<:Function}` containing one density per dimension.
When `f` is omitted it defaults to the constant-one function, yielding
integration weights only.
"""
function integrate_node_funs(inter::QInterpolator,
                             pdfs::Vector{<:Function};
                             f::Union{Nothing,Function}=nothing,
                             constant::Float64=1.0,
                             atol::Float64=1e-9,
                             rtol::Float64=1e-7)
    integrand = f === nothing ? _unity_function : f
    return _integrate_node_funs(inter, pdfs, integrand;
                                constant=constant, atol=atol, rtol=rtol)
end

function _evaluate_on_nodes(inter::QInterpolator, f::Function)
    d = inter.dims
    sz = ntuple(i -> length(inter.nodes[i]), d)
    values = Array{Float64}(undef, sz)
    xpt = zeros(Float64, d)
    @inbounds for I in CartesianIndices(values)
        @inbounds for j in 1:d
            xpt[j] = inter.nodes[j][I[j]]
        end
        values[I] = Float64(f(xpt))
    end
    return values
end

struct QIntegrator{I<:QInterpolator, W<:AbstractArray{Float64}, V}
    interpolator::I
    node_weights::W
    node_values::V
    function QIntegrator(inter::QInterpolator,
                        node_weights::AbstractArray{<:Real},
                        node_values)
        weights = node_weights isa AbstractArray{Float64} ? node_weights : Array{Float64}(node_weights)

        d = inter.dims
        ndims(weights) == d ||
            error("node_weights must be $d-D (got $(ndims(weights))-D).")
        @inbounds for i in 1:d
            size(weights, i) == length(inter.nodes[i]) ||
                error("node_weights size mismatch along dim $i: got $(size(weights,i)), expected $(length(inter.nodes[i])).")
        end

        values = node_values
        if node_values === nothing
            # no values stored
        elseif node_values isa AbstractArray{<:Real}
            ndims(node_values) >= d ||
                error("node_values must have at least $d dimensions (got $(ndims(node_values))).")
            @inbounds for i in 1:d
                size(node_values, i) == size(weights, i) ||
                    error("node_values size mismatch along dim $i: got $(size(node_values,i)), expected $(size(weights,i)).")
            end
            values = node_values isa AbstractArray{Float64} ? node_values : Array{Float64}(node_values)
        else
            error("node_values must be either nothing or an array of real numbers.")
        end

        return new{typeof(inter), typeof(weights), typeof(values)}(inter, weights, values)
    end
end

function QIntegrator(inter::QInterpolator,
                    pdfs::Vector{<:Function};
                    f::Union{Nothing,Function}=nothing,
                    constant::Float64=1.0,
                    atol::Float64=1e-9,
                    rtol::Float64=1e-7)
    weights = integrate_node_funs(inter, pdfs;
                                  constant=constant,
                                  atol=atol,
                                  rtol=rtol)
    stored = f === nothing ? nothing : _evaluate_on_nodes(inter, f)
    return QIntegrator(inter, weights, stored)
end

"""
    QIntegrator(inter::QInterpolator,
                distributions::Vector;
                f::Union{Nothing,Function}=nothing,
                atol::Float64=1e-9,
                rtol::Float64=1e-7)

Accepts a vector of distribution-like objects (e.g. `Vector{QDistribution}`).
Each entry must supply `:pdf` and `:normalization_constant` properties. Their
normalization constants are multiplied into the resulting weights.
"""
function QIntegrator(inter::QInterpolator,
                    distributions::Vector;
                    f::Union{Nothing,Function}=nothing,
                    atol::Float64=1e-9,
                    rtol::Float64=1e-7)
    isempty(distributions) && error("distributions vector must contain at least one entry.")

    pdfs = Vector{Function}(undef, length(distributions))
    scale = 1.0
    @inbounds for i in 1:length(distributions)
        dist = distributions[i]
        hasproperty(dist, :pdf) && hasproperty(dist, :normalization_constant) ||
            error("Each distribution must provide :pdf and :normalization_constant properties.")
        pdfs[i] = getproperty(dist, :pdf)
        scale *= Float64(getproperty(dist, :normalization_constant))
    end

    return QIntegrator(inter, pdfs;
                      f=f,
                      constant=scale,
                      atol=atol,
                      rtol=rtol)
end

"""
    eval_integration(inter::QInterpolator,
                     values::Array{<:Number},
                     weights::Array{<:Number})

Return the weighted sum `∑ weights[j] * values[j]`. The arrays must be dense
and aligned with the interpolation grid.
"""
function eval_integration(::QInterpolator,
                          values::Array{<:Number},
                          weights::Array{<:Number})
    T = promote_type(eltype(values), eltype(weights))
    total = zero(T)
    @inbounds @simd for idx in eachindex(values, weights)
        total += T(values[idx]) * T(weights[idx])
    end
    return total
end

function (int::QIntegrator)(values::AbstractArray{<:Number})
    values isa Array ||
        error("Integrator expects a dense Array of node values matching the interpolation grid.")
    return eval_integration(int.interpolator, values, int.node_weights)
end

function (int::QIntegrator)()
    values = int.node_values
    values === nothing &&
        error("Integrator constructed without stored node values; supply a values array or provide f at construction.")
    return eval_integration(int.interpolator, values, int.node_weights)
end

function (int::QIntegrator)(f::Function)
    values = _evaluate_on_nodes(int.interpolator, f)
    return eval_integration(int.interpolator, values, int.node_weights)
end

end # module QIntegrators
