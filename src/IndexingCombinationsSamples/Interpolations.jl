module QInterpolations

export Interpolator, eval_interpolation, nodes, basis_values, basis_values!

using LoopVectorization
using LinearAlgebra

# --- helpers ---
uniform_nodes(n::Int, a::Float64, b::Float64, endpoints::Bool) =
    n == 1 ? [0.5*(a+b)] :
    endpoints ? collect(range(a, stop=b, length=n)) :
                [a + (k - 0.5) * ((b - a) / n) for k in 1:n]

chebyshev_nodes(n::Int, a::Float64, b::Float64, endpoints::Bool) = begin
    if n == 1
        return [0.5*(a+b)]
    end
    xs = endpoints ?
        [cos(pi * (k - 1) / (n - 1)) for k in 1:n] :        # 2nd-kind (includes endpoints)
        [cos(pi * (2k - 1) / (2n)) for k in 1:n]            # 1st-kind (no endpoints)
    sort!(xs)
    c, r = 0.5*(a+b), 0.5*(b-a)
    [c + r*x for x in xs]
end

# Denominators: ∏_{k≠j} (x[j] - x[k])
function basis_denoms(x::Vector{Float64})
    n = length(x)
    n == 1 && return [1.0]
    d = Vector{Float64}(undef, n)
    @tturbo for j in 1:n
        xj = x[j]
        p  = 1.0
        for k in 1:n
            p *= ifelse(k == j, 1.0, xj - x[k])
        end
        d[j] = p
    end
    d
end

# Per-dimension scratch buffers to avoid allocations
mutable struct Scratch
    diff::Vector{Float64}
    left::Vector{Float64}
    right::Vector{Float64}
end

# In-place: fills scratch.left with ∏_{k≠j}(xval - x[k])
@inline function basis_products!(sc::Scratch, x::Vector{Float64}, xval::Float64)
    n = length(x)
    if n == 1
        sc.left[1] = 1.0
        return sc.left
    end
    d, left, right = sc.diff, sc.left, sc.right
    @turbo for i in 1:n
        d[i] = xval - x[i]
    end
    left[1] = 1.0
    @inbounds for i in 2:n
        left[i] = left[i-1] * d[i-1]
    end
    right[n] = 1.0
    @inbounds for i in n-1:-1:1
        right[i] = right[i+1] * d[i+1]
    end
    @turbo for i in 1:n
        left[i] *= right[i]      # left now holds the products
    end
    left
end

# ---------- helpers for PDF-driven nodes ----------
@inline _sqrtpdf(pdf::Function, x::Float64) = sqrt(max(float(pdf(x)), eps()))

# candidate grid for searching nodes (dense, uniform or midpoint)
@inline function _candidate_grid(n::Int, a::Float64, b::Float64; M::Int=0, endpoints::Bool=true)
    M == 0 && (M = max(32*n, 1024))              # dense by default
    if endpoints
        collect(range(a, b; length=M))           # includes endpoints
    else
        h = (b - a) / M
        [a + (k - 0.5)*h for k in 1:M]           # midpoints (no endpoints)
    end
end

"""
    leja_nodes(pdf, n, a, b; M=0, endpoints=true, tol=1e-12)

Weighted **Leja** nodes on `[a,b]` for weight `pdf(x)`.

Greedy rule over a dense candidate grid:
maximize `sqrt(pdf(x)) * ∏|x - x_j|` (in the log-domain).
Returns `n` sorted nodes. If `endpoints=true` and `n≥2`, the interval endpoints are forced in.
"""
function leja_nodes(pdf::Function, n::Int, a::Float64, b::Float64;
                    M::Int=0, endpoints::Bool=true, tol::Float64=1e-12)
    xs = _candidate_grid(n, a, b; M=M, endpoints=endpoints)
    m  = length(xs)
    w  = [ _sqrtpdf(pdf, xi) for xi in xs ]
    lw = log.(w .+ eps())

    chosen   = Float64[]
    selected = falses(m)
    logs     = zeros(Float64, m)  # log of product terms

    if endpoints && n >= 2
        ia = argmin(abs.(xs .- a)); ib = argmin(abs.(xs .- b))
        for idx in (ia, ib)
            if !selected[idx]
                xv = xs[idx]
                push!(chosen, xv); selected[idx] = true
                @. logs += log(abs(xs - xv) + eps())
            end
        end
    else
        i0 = argmax(w)
        push!(chosen, xs[i0]); selected[i0] = true
        @. logs += log(abs(xs - xs[i0]) + eps())
    end

    while length(chosen) < n
        scores = lw .+ logs
        @inbounds for i in 1:m
            selected[i] && (scores[i] = -Inf)
        end
        inext = argmax(scores)
        xnew  = xs[inext]
        push!(chosen, xnew); selected[inext] = true
        @. logs += log(abs(xs - xnew) + eps())
    end

    sort!(chosen)
    return chosen
end

"""
    fekete_nodes(pdf, n, a, b; M=0, endpoints=true, tol=1e-12)

Approximate **weighted Fekete** nodes on `[a,b]` for weight `pdf(x)` using
pivoted-QR on a **weighted Vandermonde** over a dense candidate grid.
If `endpoints=true` and `n≥2`, the interval endpoints are forced in.
"""
function fekete_nodes(pdf::Function, n::Int, a::Float64, b::Float64;
                      M::Int=0, endpoints::Bool=true, tol::Float64=1e-12)
    xs_all = _candidate_grid(n, a, b; M=M, endpoints=endpoints)

    if endpoints && n >= 2
        interior = [x for x in xs_all if !(isapprox(x, a; atol=tol) || isapprox(x, b; atol=tol))]
        nin = n - 2
        chosen = Vector{Float64}(undef, n)
        chosen[1] = a; chosen[end] = b

        Mx = length(interior)
        c, r = (a + b)/2, (b - a)/2
        V = Array{Float64}(undef, Mx, nin)
        for j in 1:Mx
            t = (interior[j] - c) / r
            V[j,1] = 1.0
            @inbounds for k in 2:nin
                V[j,k] = V[j,k-1] * t
            end
        end
        w = [ _sqrtpdf(pdf, x) for x in interior ]
        @inbounds for j in 1:Mx
            @views V[j, :] .*= w[j]
        end
        F = transpose(V)                 # nin × Mx
        Ffac = qr(F, Val(true))          # column-pivoted QR
        pidx = Vector(Ffac.p)
        chosen[2:end-1] = sort(interior[pidx[1:nin]])
        sort!(chosen)
        return chosen
    else
        xs = xs_all
        Mx = length(xs)
        c, r = (a + b)/2, (b - a)/2
        V = Array{Float64}(undef, Mx, n)
        for j in 1:Mx
            t = (xs[j] - c) / r
            V[j,1] = 1.0
            @inbounds for k in 2:n
                V[j,k] = V[j,k-1] * t
            end
        end
        w = [ _sqrtpdf(pdf, x) for x in xs ]
        @inbounds for j in 1:Mx
            @views V[j, :] .*= w[j]
        end
        F = transpose(V)                 # n × Mx
        Ffac = qr(F, Val(true))
        pidx = Vector(Ffac.p)
        chosen = xs[pidx[1:n]]
        sort!(chosen)
        return chosen
    end
end

# default flat pdf
_default_pdf(x::Float64) = 1.0
_default_pdf(x) = 1.0

# normalize a pdf or vector of pdfs to a vector length == dims
function _pdfvec(pdfs, dims::Int)
    if pdfs === nothing
        return [ _default_pdf for _ in 1:dims ]
    elseif pdfs isa Function
        return [ pdfs for _ in 1:dims ]
    else
        v = collect(pdfs)
        length(v) == dims || error("length(pdfs)=$(length(v)) must equal number of dimensions=$dims")
        return v
    end
end

# ---------- Unified constructor ----------
"""
    Interpolator(method::Symbol,
                 params::AbstractVector{<:Tuple{<:Integer,<:Real,<:Real}};
                 pdfs=nothing,
                 endpoints::Bool=true,
                 M::Int=0)

Build a tensor grid across modes:

- `method` ∈ `:uniform`, `:chebychev`, `:Leja`, `:Fekete` (case-insensitive).
- `params` is per-dimension `(n, a, b)`.
- `pdfs`: optional 1D pdf or vector of pdfs (used only for `:Leja`/`:Fekete`; defaults to flat ρ≡1).
- `endpoints`: include interval endpoints (`true`) or not (`false`) consistently in all modes.
- `M`: candidate grid size (pdf modes; `0` ⇒ auto `max(32n, 1024)`).

Precomputes Lagrange denominators and allocates scratch buffers.
"""
struct Interpolator
    dims::Int
    bounds::Vector{Tuple{Float64,Float64}}
    nodes::Vector{Vector{Float64}}
    denom::Vector{Vector{Float64}}
    work::Vector{Scratch}
    method::Symbol
    endpoints::Bool
end

function Interpolator(method::Symbol,
                      params::AbstractVector{<:Tuple{<:Integer,<:Real,<:Real}};
                      pdfs=nothing,
                      endpoints::Bool=true,
                      M::Int=0)
    m = Symbol(lowercase(String(method)))
    m ∈ (:uniform, :chebychev, :leja, :fekete) || error("method must be one of :uniform, :chebychev, :leja, :fekete")

    pnorm = [(Int(n), Float64(a), Float64(b)) for (n,a,b) in params]
    dims  = length(pnorm)

    pdfv = pdfs === nothing ? [nothing for _ in 1:dims] :
           (pdfs isa Function ? [pdfs for _ in 1:dims] : collect(pdfs))
    length(pdfv) == dims || error("length(pdfs) must equal number of dimensions")

    bounds = Vector{Tuple{Float64,Float64}}(undef, dims)
    nodesv = Vector{Vector{Float64}}(undef, dims)
    denom  = Vector{Vector{Float64}}(undef, dims)
    work   = Vector{Scratch}(undef, dims)

    for i in 1:dims
        n, a, b = pnorm[i]
        bounds[i] = (a, b)
        ρ = pdfv[i] === nothing ? (x->1.0) : pdfv[i]

        xi = if m === :uniform
            uniform_nodes(n, a, b, endpoints)
        elseif m === :chebychev
            chebyshev_nodes(n, a, b, endpoints)
        elseif m === :leja
            leja_nodes(ρ, n, a, b; M=M, endpoints=endpoints)
        else # :fekete
            fekete_nodes(ρ, n, a, b; M=M, endpoints=endpoints)
        end

        nodesv[i] = xi
        denom[i]  = basis_denoms(xi)
        work[i]   = Scratch(zeros(n), ones(n), ones(n))
    end

    return Interpolator(dims, bounds, nodesv, denom, work, method, endpoints)
end

# ---------------------------
# nodes (with separate grids option)
# ---------------------------
"""
    nodes(inter::Interpolator; separate::Bool=false)

Return the grid nodes.

- If `separate=false` (default): returns a `Vector{Vector{Float64}}` with each entry a
  coordinate vector `[x₁, x₂, …, x_d]` over the full Cartesian grid.

- If `separate=true`: returns a `NTuple{d,Array{Float64,d}}` giving per-dimension coordinate
  grids suitable for plotting (`X, Y, ...`). Each array has the same size as the grid.
"""
function nodes(inter::Interpolator; separate::Bool=false)
    d   = inter.dims
    axs = inter.nodes
    sz  = ntuple(i -> length(axs[i]), d)

    if !separate
        total = prod(sz)
        out = Vector{Vector{Float64}}(undef, total)
        k = 1
        @inbounds for I in CartesianIndices(sz)
            p = Vector{Float64}(undef, d)
            @inbounds for i in 1:d
                p[i] = axs[i][I[i]]
            end
            out[k] = p
            k += 1
        end
        return out
    else
        out = ntuple(i -> Array{Float64}(undef, sz), d)
        @inbounds for I in CartesianIndices(sz)
            for j in 1:d
                out[j][I] = axs[j][I[j]]
            end
        end
        return out
    end
end

# ---------------------------
# evaluator
# ---------------------------
"""
    eval_interpolation(inter::Interpolator,
                       pos::AbstractVector{<:Real},
                       values::AbstractArray{<:Real,N}) where {N}

Evaluate the interpolant at `pos` using node values `values`.

- `values` must be an `N`-D array with `N == inter.dims`, and `size(values,i) == length(inter.nodes[i])`.
- Allocation-free w.r.t. problem size (uses internal scratch buffers).
"""
function eval_interpolation(inter::Interpolator,
                            pos::AbstractVector{<:Real},
                            values::AbstractArray{<:Real,N}) where {N}
    d = inter.dims
    N == d || error("values must be $d-D (got $N-D)")
    length(pos) == d || error("pos must have length $d")
    @inbounds for i in 1:d
        size(values, i) == length(inter.nodes[i]) ||
            error("values size($(size(values))) does not match grid sizes " *
                  "along dim $i: got $(size(values,i)), expected $(length(inter.nodes[i]))")
    end

    # per-dimension Lagrange weights
    Wref = Vector{Vector{Float64}}(undef, d)
    @inbounds for i in 1:d
        xi  = inter.nodes[i]
        num = basis_products!(inter.work[i], xi, Float64(pos[i]))  # writes into work[i].left
        di  = inter.denom[i]
        @turbo for j in 1:length(xi)
            num[j] = num[j] / di[j]
        end
        Wref[i] = num
    end

    s = 0.0
    @inbounds for I in CartesianIndices(values)
        w = 1.0
        @inbounds for i in 1:d
            w *= Wref[i][I[i]]
        end
        s += Float64(values[I]) * w
    end
    s
end

# ---------------------------
# Per-node basis evaluation
# ---------------------------
"""
    basis_values!(out::AbstractArray{Float64,N},
                  inter::Interpolator,
                  pos::AbstractVector{<:Real}) where {N}

Compute **all tensor-product Lagrange basis functions** at location `pos` and
store them into `out` (same shape as the node grid).

- `out` must have `N == inter.dims` and `size(out,i) == length(inter.nodes[i])`.

Allocation-free aside from small temporaries.
"""
function basis_values!(out::AbstractArray{Float64,N},
                       inter::Interpolator,
                       pos::AbstractVector{<:Real}) where {N}
    d = inter.dims
    N == d || error("out must be $d-D (got $N-D)")
    length(pos) == d || error("pos must have length $d")
    @inbounds for i in 1:d
        size(out, i) == length(inter.nodes[i]) ||
            error("out shape mismatch along dim $i")
    end

    # per-dimension weights
    Wref = Vector{Vector{Float64}}(undef, d)
    @inbounds for i in 1:d
        xi  = inter.nodes[i]
        num = basis_products!(inter.work[i], xi, Float64(pos[i]))
        di  = inter.denom[i]
        @turbo for j in 1:length(xi)
            num[j] = num[j] / di[j]
        end
        Wref[i] = num
    end

    @inbounds for I in CartesianIndices(out)
        w = 1.0
        @inbounds for i in 1:d
            w *= Wref[i][I[i]]
        end
        out[I] = w
    end
    return out
end

"""
    basis_values(inter::Interpolator, pos::AbstractVector{<:Real})

Allocate and return an array containing **all basis function values** at `pos`.
This is a convenience wrapper around [`basis_values!`](@ref).
"""
function basis_values(inter::Interpolator, pos::AbstractVector{<:Real})
    sz = ntuple(i -> length(inter.nodes[i]), inter.dims)
    out = Array{Float64}(undef, sz)
    return basis_values!(out, inter, pos)
end

include("InterpolationOps/Interpolation_Integration.jl")

end # module
