export expand

# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────
"""
    expand(f::CFunction, mode::Symbol, target::Type{T}, args...) where T<:CFunction

Walk the symbolic AST `f` and apply a mode-specific expansion whenever a node of type `target` is encountered.

# Modes
- `:Taylor`   — Taylor-series expansion on `CExp` or `CLog` nodes (needs `order::Int`).
- `:Rational` — Distribute `CRational` over sums in the numerator.
- `:Log`      — Distribute `log` over products/quotients (you can omit `target` with `expand(f, Val(:Log))`).
"""
function expand(f::CFunction, mode::Symbol, target::Type{T}, args...) where T<:CFunction
    newf, _ = _expand(f, Val(mode), target, args...)
    return newf
end

"Shorthand for `expand(f, :Log, CLog)`."
function expand(f::CFunction, ::Val{:Log})
    newf, _ = _expand(f, Val(:Log), CLog)
    return newf
end

# ──────────────────────────────────────────────────────────────────────────────
# Generic tree walk (default: recurse & rebuild if children changed)
# ──────────────────────────────────────────────────────────────────────────────

_expand(s::CAtom,     ::Val{M}, ::Type{T}, args...) where {M,T} = (s, false)
_expand(f::CFunction, ::Val{M}, ::Type{T}, args...) where {M,T} = (f, false)  # fallback

# Sums
function _expand(s::CSum, ::Val{M}, ::Type{T}, args...) where {M,T}
    terms2 = CFunction[]
    any_exp = false
    @inbounds for term in s.terms
        t2, e = _expand(term, Val(M), T, args...)
        if t2 isa CSum
            append!(terms2, t2.terms)
        else
            push!(terms2, t2)
        end
        any_exp |= e
    end
    return any_exp ? (_CSum(terms2), true) : (s, false)
end

# Products
function _expand(p::CProd, ::Val{M}, ::Type{T}, args...) where {M,T}
    terms2 = Vector{CFunction}(undef, length(p.terms))
    any_exp = false
    @inbounds for (i, term) in enumerate(p.terms)
        t2, e = _expand(term, Val(M), T, args...)
        terms2[i] = t2
        any_exp |= e
    end
    return any_exp ? (CProd(p.coeff, terms2, Val(:nosimp)), true) : (p, false)
end

# Rationals
function _expand(r::CRational, ::Val{M}, ::Type{T}, args...) where {M,T}
    n2, e1 = _expand(r.numer, Val(M), T, args...)
    d2, e2 = _expand(r.denom, Val(M), T, args...)
    any_exp = e1 || e2
    return any_exp ? (CRational(n2, d2, Val(:nosimp)), true) : (r, false)
end

# Exp / Log
function _expand(e::CExp, ::Val{M}, ::Type{T}, args...) where {M,T}
    x2, ch = _expand(e.x, Val(M), T, args...)
    return ch ? (CExp(e.coeff, x2, Val(:nosimp)), true) : (e, false)
end
function _expand(l::CLog, ::Val{M}, ::Type{T}, args...) where {M,T}
    x2, ch = _expand(l.x, Val(M), T, args...)
    return ch ? (CLog(l.coeff, x2, Val(:nosimp)), true) : (l, false)
end

# Power
function _expand(p::CPower, ::Val{M}, ::Type{T}, args...) where {M,T}
    x2, ch = _expand(p.x, Val(M), T, args...)
    return ch ? (CPower(p.coeff, x2, p.exponent, Val(:nosimp)), true) : (p, false)
end

# Vector / Matrix (descend elementwise; keep coeff/orientation/shape)
function _expand(v::CVector, ::Val{M}, ::Type{T}, args...) where {M,T}
    ent2 = Vector{CFunction}(undef, length(v.expr))
    any_exp = false
    @inbounds for i in eachindex(v.expr)
        t2, e = _expand(v.expr[i], Val(M), T, args...)
        ent2[i] = t2
        any_exp |= e
    end
    return any_exp ? (CVector(v.coeff, ent2; row=v.row), true) : (v, false)
end
function _expand(A::CMatrix, ::Val{M}, ::Type{T}, args...) where {M,T}
    m, n = size(A.expr)
    flat = Vector{CFunction}(undef, length(A.expr))
    any_exp = false
    @inbounds for (k, e) in enumerate(A.expr)
        t2, ch = _expand(e, Val(M), T, args...)
        flat[k] = t2
        any_exp |= ch
    end
    return any_exp ? (CMatrix(A.coeff, reshape(flat, m, n)), true) : (A, false)
end

# ──────────────────────────────────────────────────────────────────────────────
# Specialized modes
# ──────────────────────────────────────────────────────────────────────────────

# ---- :Taylor on CExp ----
# coeff * exp(x) ≈ sum_{n=0}^order coeff * x^n / n!
function _expand(e::CExp, ::Val{:Taylor}, ::Type{CExp}, order::Int)
    order > 0 || error("Order of Taylor expansion must be positive!")
    x2, _ = _expand(e.x, Val(:Taylor), CExp, order)

    terms = CFunction[]
    # n = 0 term: coeff * 1
    push!(terms, CAtom(ComplexRational(e.coeff.a, e.coeff.b, e.coeff.c), zeros(Int, dims(e))))

    # higher-order terms
    fac = 1
    for n in 1:order
        fac *= n
        cn = ComplexRational(e.coeff.a, e.coeff.b, e.coeff.c * fac)
        xn = (x2 ^ n)
        termn = xn * cn
        if termn isa CSum
            append!(terms, termn.terms)
        else
            push!(terms, termn)
        end
    end
    return (_CSum(terms), true)
end

# ---- :Taylor on CLog ----
# log(1+Δ) = sum_{n=1}^order (-1)^(n+1) Δ^n / n, scaled by l.coeff
function _expand(l::CLog, ::Val{:Taylor}, ::Type{CLog}, order::Int)
    order > 0 || error("Order of Taylor expansion must be positive!")
    x2, _ = _expand(l.x, Val(:Taylor), CLog, order)

    one = CAtom(ComplexRational(1,0,1), zeros(Int, dims(x2)))
    Δ = x2 - one
    curr = Δ

    terms = CFunction[]
    # n = 1
    push!(terms, curr * l.coeff)
    # n >= 2
    for n in 2:order
        curr = curr * Δ
        sgn = isodd(n) ? 1 : -1
        coef = l.coeff * ComplexRational(sgn, 0, n)
        termn = curr * coef
        if termn isa CSum
            append!(terms, termn.terms)
        else
            push!(terms, termn)
        end
    end
    return (_CSum(terms), true)
end

# ---- :Rational on CRational ----
# (A+B+…)/D  ⇒  A/D + B/D + …
function _expand(r::CRational, ::Val{:Rational}, ::Type{CRational}, args...)
    n2, _ = _expand(r.numer, Val(:Rational), CRational, args...)
    d2, _ = _expand(r.denom, Val(:Rational), CRational, args...)
    if n2 isa CSum
        out = CFunction[]
        @inbounds for t in n2.terms
            push!(out, t / d2)
        end
        return (_CSum(out), true)
    else
        return (r === CRational(n2, d2) ? r : CRational(n2, d2), (n2 !== r.numer) || (d2 !== r.denom))
    end
end

# ---- :Log on CLog ----
# log(prod) + … - log(den) + q*log(x) for powers; also includes product coeff
@inline _qcr(q::Rational{Int}) = ComplexRational(numerator(q), 0, denominator(q))

function _log_collect_terms(x::CFunction, acc::Vector{Tuple{ComplexRational,CFunction}})
    if x isa CProd
        # coefficient factor
        if !isone(x.coeff)
            push!(acc, (ComplexRational(1,0,1), CAtom(x.coeff, zeros(Int, dims(x)))))
        end
        @inbounds for t in x.terms
            _log_collect_terms(t, acc)
        end
    elseif x isa CRational
        _log_collect_terms(x.numer, acc)
        # subtract logs of denominator
        tmp = Tuple{ComplexRational,CFunction}[]
        _log_collect_terms(x.denom, tmp)
        for (c,f) in tmp
            push!(acc, (c * ComplexRational(-1,0,1), f))
        end
    elseif x isa CPower
        # q * log(base) + log(coeff) if coeff != 1
        if !isone(x.coeff)
            push!(acc, (ComplexRational(1,0,1), CAtom(x.coeff, zeros(Int, dims(x)))))
        end
        push!(acc, (_qcr(x.exponent), x.x))
    else
        push!(acc, (ComplexRational(1,0,1), x))
    end
    return acc
end

function _expand(l::CLog, ::Val{:Log}, ::Type{CLog}, args...)
    x2, _ = _expand(l.x, Val(:Log), CLog, args...)
    pairs = _log_collect_terms(x2, Tuple{ComplexRational,CFunction}[])
    logs = CFunction[]
    for (mul, f) in pairs
        push!(logs, CLog(l.coeff * mul, f))
    end
    return (_CSum(logs), true)
end
