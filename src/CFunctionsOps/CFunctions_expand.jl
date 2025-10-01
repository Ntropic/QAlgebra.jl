export expand

# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────
"""
    expand(f::CFunction, mode::Symbol, target::Type{T}; order=0) where T<:CFunction

Walk the symbolic AST `f` and apply a mode-specific expansion whenever a node of type `target` is encountered.

# Modes
- `:Taylor`   — Taylor-series expansion on `CExp` or `CLog` nodes (needs `order`).
- `:Rational` — Distribute `CRational` over sums in the numerator.
- `:Log`      — Distribute `log` over products/quotients (you can omit `target` with `expand(f, Val(:Log))`).
"""
function expand(f::CFunction, mode::Symbol, target::Type{T}; order::Int=0) where T<:CFunction
    newf, _ = _expand(f, Val(mode), Val(Symbol(target)), order)
    return newf
end
function expand(f::CFunction, mode::Symbol, target::Symbol; order::Int=0)
    newf, _ = _expand(f, Val(mode), Val(target), order)
    return newf
end

function expand(f::CFunction, mode::Symbol; order::Int=0)
    if mode in (:CExp, :Exp)
        return expand(f, Val(:Exp); order=order)
    elseif mode in (:CLog, :Log)
        return expand(f, Val(:Log); order=order)
    elseif mode in (:CRational, :Rational)
        return expand(f, Val(:Rational); order=order)
    end
    newf, _ = _expand(f, Val(mode), Val(:Any), order)
    return newf
end

function expand(f::CFunction, ::Val{:Taylor}, ::Val{:Log}; order::Int=2)
    ord = order > 0 ? order : 2
    newf, _ = _expand(f, Val(:Taylor), Val(:Log), ord)
    return newf
end
function expand(f::CFunction, ::Val{:Log}; order::Int=2)
    ord = order > 0 ? order : 2
    newf, _ = _expand(f, Val(:Log), Val(:Log), ord)
    return newf
end
function expand(f::CFunction, ::Val{:Exp}; order::Int=2)
    ord = order > 0 ? order : 2
    newf, _ = _expand(f, Val(:Taylor), Val(:CExp), ord)
    return newf
end
function expand(f::CFunction, ::Val{:Taylor}, ::Val{:Exp}; order::Int=2)
    ord = order > 0 ? order : 2
    newf, _ = _expand(f, Val(:Taylor), Val(:CExp), ord)
    return newf
end
function expand(f::CFunction, ::Val{:Rational}; order::Int=0)
    newf, _ = _expand(f, Val(:Rational), Val(:CRational), order)
    return newf
end

# ──────────────────────────────────────────────────────────────────────────────
# Generic tree walk (default: recurse & rebuild if children changed)
# ──────────────────────────────────────────────────────────────────────────────

_expand(f::CFunction, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T} = (f, false)  # fallback includes CAtom and CAbstracht

# Sums
function _expand(s::CSum, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T}
    terms2 = CFunction[]
    any_exp = false
    @inbounds for term in s.terms
        t2, e = _expand(term, Val(M), Val(T), order, args...)
        if t2 isa CSum
            append!(terms2, t2.terms)
        else
            push!(terms2, t2)
        end
        any_exp |= e
    end
    return any_exp ? (_CSum(s.param_info, terms2), true) : (s, false)
end

# Products
function _expand(p::CProd, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T}
    terms2 = Vector{CFunction}(undef, length(p.terms))
    any_exp = false
    @inbounds for (i, term) in enumerate(p.terms)
        t2, e = _expand(term, Val(M), Val(T), order, args...)
        terms2[i] = t2
        any_exp |= e
    end
    return any_exp ? (CProd(p.param_info, p.coeff, terms2, Val(:nosimp)), true) : (p, false)
end

# Rationals
function _expand(r::CRational, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T}
    n2, e1 = _expand(r.numer, Val(M), Val(T), order, args...)
    d2, e2 = _expand(r.denom, Val(M), Val(T), order, args...)
    any_exp = e1 || e2
    return any_exp ? (CRational(r.param_info, n2, d2, Val(:nosimp)), true) : (r, false)
end

# Exp / Log
function _expand(e::CExp, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T}
    x2, ch = _expand(e.expr, Val(M), Val(T), order, args...)
    return ch ? (CExp(e.param_info, e.coeff, x2, Val(:nosimp)), true) : (e, false)
end
function _expand(l::CLog, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T}
    x2, ch = _expand(l.expr, Val(M), Val(T), order, args...)
    return ch ? (CLog(l.param_info, l.coeff, x2, Val(:nosimp)), true) : (l, false)
end

# Power
function _expand(p::CPower, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T}
    x2, ch = _expand(p.expr, Val(M), Val(T), order, args...)
    return ch ? (CPower(p.param_info, p.coeff, x2, p.exponent, Val(:nosimp)), true) : (p, false)
end

# Vector / Matrix (descend elementwise; keep coeff/orientation/shape)
function _expand(v::CVector, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T}
    ent2 = Vector{CFunction}(undef, length(v.expr))
    any_exp = false
    @inbounds for i in eachindex(v.expr)
        t2, e = _expand(v.expr[i], Val(M), Val(T), order, args...)
        ent2[i] = t2
        any_exp |= e
    end
    return any_exp ? (CVector(v.param_info, v.coeff, ent2; row=v.row), true) : (v, false)
end
function _expand(A::CMatrix, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T}
    m, n = size(A.expr)
    flat = Vector{CFunction}(undef, length(A.expr))
    any_exp = false
    @inbounds for (k, e) in enumerate(A.expr)
        t2, ch = _expand(e, Val(M), Val(T), order, args...)
        flat[k] = t2
        any_exp |= ch
    end
    return any_exp ? (CMatrix(A.param_info, A.coeff, reshape(flat, m, n)), true) : (A, false)
end

########################################################################################################################################################################
# just expand the arguments 
const DefinitionAliases = (:Def, :def, :Self, :self)
function insert(a::CCustomType)
    return insert(a, a.expr)
end
function insert(A::CCustomType, expr::Vector{CFunction}) 
    fun = A.ctype_def.fun
    index_map = A.ctype_def.index_map
    return insert(fun, index_map, expr)
end

insert(fun::CAtom, index_map::Vector{Int}, expr::Vector{CFunction}) = fun, false 
function insert(fun::CAbstract, index_map::Vector{Int}, expr::Vector{CFunction})
    if fun.index > length(index_map) || index_map[fun.index] == 0
        return fun, false 
    else
        return expr[index_map[fun.index]], true
    end
end
function insert(fun::CComposite, index_map::Vector{Int}, expr::Vector{CFunction})
    n, c = insert(fun.expr, index_map, expr) 
    if c 
        return modify_expr(fun, n), c 
    else
        return fun, false 
    end
end
function insert(fun::CMultiComposite, index_map::Vector{Int}, expr::Vector{CFunction})
    new_fun::Vector{CFunction} = []
    any_changed = false 
    for expr in fun.expr
        n, c = insert(fun.expr, index_map, expr) 
        push!(new_fun, n)
        any_changed |= c 
    end
    if any_changed 
        return modify_expr(fun, new_fun), c 
    else
        return fun, false 
    end
end
function insert(fun::CCustomType, index_map::Vector{Int}, expr::Vector{CFunction})
    new_fun::Vector{CFunction} = []
    any_changed = false 
    for expr in fun.expr
        n, c = insert(fun.expr, index_map, expr) 
        push!(new_fun, n)
        any_changed |= c 
    end
    if any_changed 
        return modify_expr(fun, new_fun), c 
    else
        return fun, false 
    end
end
########################################################################################################################################################################

function _expand(A::CCustomType, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T}
    new_expr::Vector{CFunction} = []>
    if M ∈ DefinitionAliases 
        if T ∈ A.ctype_def.type_symbols
            return insert(A, A.expr)
        else
            if A.ctype_def.has_abstract  # expand the substituted expression 
                any_changed = false
                for expr in new_expr 
                    n, c = _expand(expr, Val(M), Val(T), order, args...)
                    push!(new_expr, n)
                    any_changed |= c
                end
                return A.coeff*new_expr, any_changed
            else
                # has only its expression no substitutions involved. 
                n, c = _expand(expr, Val(M), Val(T), order, args...)
                if c 
                    return A.coeff*n, c 
                else
                    return A, false
                end
            end
        end
    elseif A.ctype_def.has_abstract  # expand the substituted expression 
        any_changed = false
        for expr in new_expr 
            n, c = _expand(expr, Val(M), Val(T), order, args...)
            push!(new_expr, n)
            any_changed |= c
        end
        return new_expr, any_changed
    else
        # has only its expression no substitutions involved. 
        n, c = _expand(expr, Val(M), Val(T), order, args...)
        if c 
            return n, c 
        else
            return A, false
        end
    end
end
        

# ──────────────────────────────────────────────────────────────────────────────
# Specialized modes
# ──────────────────────────────────────────────────────────────────────────────
const ExpAliases = (:exp, :Exp, :CExp, :Any, :any)
const LogAliases = (:log, :Log, :CLog, :Any, :any)
const RationalAliases = (:rational, :Rational, :CRational, :Any, :any)


# ---- :Taylor on CExp ----
# coeff * exp(x) ≈ sum_{n=0}^order coeff * x^n / n!
function _expand(e::CExp, ::Val{:Taylor}, ::Val{S}, order::Int) where {S}
    S ∈ ExpAliases || return (e, false)
    ord = order > 0 ? order : 2

    # recurse with the SAME alias S
    x2, _ = _expand(e.expr, Val(:Taylor), Val(S), ord)

    terms = CFunction[]
    # n = 0 term: coeff * 1
    curr_coeff = e.coeff
    push!(terms, CAtom(e.param_info, curr_coeff, spzeros(Int, dims(e))))

    # higher-order terms
    #fac = 1
    for n in 1:ord
        #fac *= n
        curr_coeff /= n
        xn = (x2 ^ n)
        termn = xn * curr_coeff
        if termn isa CSum
            append!(terms, termn.expr)
        else
            push!(terms, termn)
        end
    end
    return (_CSum(e.param_info, terms, Val(:nosimp)), true)
end

# ---- :Taylor on CLog ----
# log(1+Δ) = sum_{n=1}^order (-1)^(n+1) Δ^n / n, scaled by l.coeff
function _expand(l::CLog, ::Val{:Taylor}, ::Val{S}, order::Int) where {S}
    S ∈ LogAliases || return (l, false)
    ord = order > 0 ? order : 2

    # recurse with SAME alias
    x2, _ = _expand(l.expr, Val(:Taylor), Val(S), ord)

    one = CAtom(l.param_info, ComplexRational(1,0,1), spzeros(Int, dims(x2)))
    Δ = x2 - one
    curr = Δ

    terms = CFunction[]
    # n = 1
    push!(terms, curr * l.coeff)
    # n >= 2
    for n in 2:ord
        curr = curr * Δ
        sgn = isodd(n) ? 1 : -1
        coef = l.coeff * ComplexRational(sgn, 0, n)
        termn = curr * coef
        if termn isa CSum
            append!(terms, termn.expr)
        else
            push!(terms, termn)
        end
    end
    return (_CSum(l.param_info, terms, Val(:nosimp)), true)
end

# ---- :Rational on CRational ----
# (A+B+…)/D  ⇒  A/D + B/D + …
function _expand(r::CRational, ::Val{:Rational}, ::Val{S}, order::Int, args...) where {S}
    S ∈ RationalAliases || return (r, false)

    n2, _ = _expand(r.numer, Val(:Rational), Val(S), order, args...)
    d2, _ = _expand(r.denom, Val(:Rational), Val(S), order, args...)
    if n2 isa CSum
        out = CFunction[]
        @inbounds for t in n2.expr
            push!(out, t / d2)
        end
        return (_CSum(r.param_info, out, Val(:nosimp)), true)
    else
        return (CRational(r.param_info, n2, d2, Val(:nosimp)), 
                (n2 !== r.numer) || (d2 !== r.denom))
    end
end

# ---- :Log on CLog ----
# log(prod) + … - log(den) + q*log(x) for powers; also includes product coeff
@inline _qcr(q::Rational{Int}) = ComplexRational(numerator(q), 0, denominator(q))

function _log_collect_terms(x::CFunction, acc::Vector{Tuple{ComplexRational,CFunction}})
    if x isa CProd
        # coefficient factor
        if !isone(x.coeff)
            push!(acc, (ComplexRational(1,0,1), CAtom(x.param_info, x.coeff, spzeros(Int, dims(x)))))
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
            push!(acc, (ComplexRational(1,0,1), CAtom(x.param_info, x.coeff, spzeros(Int, dims(x)))))
        end
        push!(acc, (_qcr(x.exponent), x.expr))
    else
        push!(acc, (ComplexRational(1,0,1), x))
    end
    return acc
end

function _expand(l::CLog, ::Val{:Log}, ::Val{S}, args...) where {S}
    S ∈ LogModeAliases || return (l, false)

    # recurse with SAME alias
    x2, _ = _expand(l.expr, Val(:Log), Val(S), args...)
    pairs = _log_collect_terms(x2, Tuple{ComplexRational,CFunction}[])
    logs = CFunction[]
    for (mul, f) in pairs
        push!(logs, CLog(l.param_info, l.coeff * mul, f))
    end
    return (_CSum(l.param_info, logs, Val(:nosimp)), true)
end
