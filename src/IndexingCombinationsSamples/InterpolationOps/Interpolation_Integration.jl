import QuadGK
export normalization_constant, node_functionals_quadgk_array

function normalization_constant(pdfs::AbstractVector{<:Function}, inter::Interpolator;
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

"""
    node_functionals_quadgk_array(inter::Interpolator,
                                  pdfs=Vector{Function};
                                  f::Function,
                                  constant::Float64=1.0,
                                  atol::Float64=1e-9,
                                  rtol::Float64=1e-7)

Compute `W[j...] = ∫ ρ(x) f(x) L_j(x) dx` over the box, **all j at once**,
using nested QuadGK with an **array-valued** integrand built by `basis_values!`.

- `f(xvec)`: nonseparable function, `xvec::Vector{Float64}` of length `inter.dims`
- `pdfs`: 1D pdf or vector of pdfs per dimension (defaults to 1)
- `normalize=true`: divides by `∏_i ∫ ρ_i` so `ρ` acts as a normalized pdf

Returns `W::Array{Float64,N}` where `N == inter.dims` and `size(W) == (length(nodes[i])...)`.
"""
function node_functionals_quadgk_array(inter::Interpolator,
                                       pdfs::AbstractVector{<:Function};
                                       f::Function,
                                       constant::Float64,
                                       atol::Float64=1e-9,
                                       rtol::Float64=1e-7)

    d = inter.dims
    @assert d == length(pdfs) "Require one pdf per integration dimension."

    # buffers: basis tensor at current x, per-dim point vector
    sz  = ntuple(i -> length(inter.nodes[i]), d)
    Lx  = Array{Float64}(undef, sz)   # will hold L(x)
    xpt = zeros(Float64, d)

    # norm for array-valued integrands
    arr_norm = x -> sqrt(sum(abs2, x))

    # recursive nested integration; returns Arrays so QuadGK integrates arrays
    function _integrate_dim(i::Int, ρaccum::Float64)
        ai, bi = inter.bounds[i]
        ρi = pdfs[i]

        if i < d
            g = function (xi::Float64)
                xpt[i] = xi
                ρnew = ρaccum * max(ρi(xi), 0.0)
                # IMPORTANT: don't integrate the next dimension here;
                # just return the inner integral value.
                return _integrate_dim(i+1, ρnew)
            end
            val, _ = QuadGK.quadgk(g, ai, bi; atol=atol, rtol=rtol, norm=arr_norm)
            return val
        else
            glast = function (xd::Float64)
                xpt[d] = xd
                basis_values!(Lx, inter, xpt)  # fill Lx in-place
                ρx = ρaccum * max(ρi(xd), 0.0)
                fx = f(xpt)
                return Lx .* (ρx * fx)         # array-valued integrand
            end
            val, _ = QuadGK.quadgk(glast, ai, bi; atol=atol, rtol=rtol, norm=arr_norm)
            return val
        end
    end

    W = _integrate_dim(1, 1.0)
    W .*= constant 
    return W
end
