import ..CFunctions: expand
using ..CFunctions: CExp, CLog, CRational, CFunction

const Q_EXP_ALIASES = (:QExp, :Exp, :exp)
const Q_LOG_ALIASES = (:QLog, :Log, :log)
const Q_POWER_ALIASES = (:QPower, :Power, :power)
const Q_ROOT_ALIASES = (:QRoot, :Root, :root, :sqrt)
const Q_COMMUTATOR_ALIASES = (:QCommutator, :Commutator, :commutator)
const C_EXP_ALIASES = (:CExp,)
const C_LOG_ALIASES = (:CLog,)
const C_RATIONAL_ALIASES = (:CRational, :Rational, :rational)

@inline _is_alias(target::Symbol, aliases)::Bool = target in aliases

@inline function _expand_coeff_fun(coeff::CFunction, mode::Symbol, target::Symbol, order::Int)
    target == :Any && return coeff
    if _is_alias(target, Q_EXP_ALIASES) || target in C_EXP_ALIASES
        ord = order > 0 ? order : 2
        return expand(coeff, :Taylor, CExp; order=ord)
    elseif _is_alias(target, Q_LOG_ALIASES) || target in C_LOG_ALIASES
        ord = order > 0 ? order : 2
        return expand(coeff, :Taylor, CLog; order=ord)
    elseif target in C_RATIONAL_ALIASES || mode == :Rational
        return expand(coeff, :Rational, CRational; order=order)
    else
        return coeff
    end
end

"""
    expand(q::QObj, mode::Symbol, target::Type{T}; order=2) where {T<:QObj}

Walk the quantum-expression AST `q` and apply a mode/target-specific expansion
whenever a node of type `target` is encountered. The interface mirrors the
`expand` API provided for coefficient functions, while dispatching on quantum
objects (`QObj`).

# Common modes
- `:Taylor`    — series expansion for quantum exponentials/logarithms/powers (uses `order`).
- `:Rational`  — distributive expansion inside coefficient functions containing `CRational` nodes.
- `:Expand`    — commutator expansion without explicitly passing the mode (see `Q_COMMUTATOR_ALIASES`).

If no rule matches the requested combination the original expression is
returned unchanged. Unless stated otherwise, expansions default to using order
two.
"""
function expand(q::QObj, mode::Symbol, target::Type{T}, args...; order::Int=2) where {T<:QObj}
    return expand(q, mode, Symbol(target), args...; order=order)
end

function expand(q::QObj, mode::Symbol, target::Symbol, args...; order::Int=2)
    if target in Q_EXP_ALIASES
        return expand(q, Val(mode), Val(:Exp), args...; order=order)
    elseif target in Q_LOG_ALIASES
        return expand(q, Val(mode), Val(:Log), args...; order=order)
    elseif target in Q_POWER_ALIASES
        return expand(q, Val(mode), Val(:Power), args...; order=order)
    elseif target in Q_ROOT_ALIASES
        return expand(q, Val(mode), Val(:Root), args...; order=order)
    elseif target in Q_COMMUTATOR_ALIASES
        return expand(q, Val(mode), Val(:QCommutator), args...; order=order)
    elseif target in C_RATIONAL_ALIASES
        return expand(q, Val(mode), Val(:CRational), args...; order=order)
    elseif target in C_EXP_ALIASES
        return expand(q, Val(mode), Val(:CExp), args...; order=order)
    elseif target in C_LOG_ALIASES
        return expand(q, Val(mode), Val(:CLog), args...; order=order)
    end
    new_q, _ = _expand(q, Val(mode), Val(target), order, args...)
    return new_q
end

function expand(q::QObj, mode::Symbol, args...; order::Int=2)
    if mode in Q_EXP_ALIASES
        return expand(q, Val(:Exp), args...; order=order)
    elseif mode in Q_LOG_ALIASES
        return expand(q, Val(:Log), args...; order=order)
    elseif mode in Q_POWER_ALIASES
        return expand(q, Val(:Power), args...; order=order)
    elseif mode in Q_ROOT_ALIASES
        return expand(q, Val(:Root), args...; order=order)
    elseif mode in Q_COMMUTATOR_ALIASES
        return expand(q, Val(:QCommutator), args...; order=order)
    elseif mode in C_RATIONAL_ALIASES
        return expand(q, Val(:CRational), args...; order=order)
    elseif mode in C_EXP_ALIASES
        return expand(q, Val(:CExp), args...; order=order)
    elseif mode in C_LOG_ALIASES
        return expand(q, Val(:CLog), args...; order=order)
    end
    new_q, _ = _expand(q, Val(mode), Val(:Any), order, args...)
    return new_q
end

function expand(q::QObj, ::Val{mode}, ::Val{target}, args...; order::Int=2) where {mode,target}
    new_q, _ = _expand(q, Val(mode), Val(target), order, args...)
    return new_q
end

function expand(q::QObj, ::Val{mode}, target::Symbol, args...; order::Int=2) where {mode}
    return expand(q, mode, target, args...; order=order)
end

function expand(q::QObj, ::Val{mode}, target::Type{T}, args...; order::Int=2) where {mode,T}
    return expand(q, mode, Symbol(target), args...; order=order)
end

function expand(q::QObj, ::Val{mode}, args...; order::Int=2) where {mode}
    new_q, _ = _expand(q, Val(mode), Val(:Any), order, args...)
    return new_q
end

function expand(q::QObj, ::Val{:Exp}, args...; order::Int=2)
    ord = order > 0 ? order : 2
    new_q, _ = _expand(q, Val(:Taylor), Val(:QExp), ord, args...)
    return new_q
end

function expand(q::QObj, ::Val{:Taylor}, ::Val{:Exp}, args...; order::Int=2)
    ord = order > 0 ? order : 2
    new_q, _ = _expand(q, Val(:Taylor), Val(:QExp), ord, args...)
    return new_q
end

function expand(q::QObj, ::Val{:Log}, args...; order::Int=2)
    ord = order > 0 ? order : 2
    new_q, _ = _expand(q, Val(:Taylor), Val(:QLog), ord, args...)
    return new_q
end

function expand(q::QObj, ::Val{:Taylor}, ::Val{:Log}, args...; order::Int=2)
    ord = order > 0 ? order : 2
    new_q, _ = _expand(q, Val(:Taylor), Val(:QLog), ord, args...)
    return new_q
end

function expand(q::QObj, ::Val{:Power}, args...; order::Int=2)
    ord = order > 0 ? order : 2
    new_q, _ = _expand(q, Val(:Taylor), Val(:QPower), ord, args...)
    return new_q
end

function expand(q::QObj, ::Val{:Taylor}, ::Val{:Power}, args...; order::Int=2)
    ord = order > 0 ? order : 2
    new_q, _ = _expand(q, Val(:Taylor), Val(:QPower), ord, args...)
    return new_q
end

function expand(q::QObj, ::Val{:Root}, args...; order::Int=2)
    ord = order > 0 ? order : 2
    new_q, _ = _expand(q, Val(:Taylor), Val(:QRoot), ord, args...)
    return new_q
end

function expand(q::QObj, ::Val{:Taylor}, ::Val{:Root}, args...; order::Int=2)
    ord = order > 0 ? order : 2
    new_q, _ = _expand(q, Val(:Taylor), Val(:QRoot), ord, args...)
    return new_q
end

function expand(q::QObj, ::Val{:QCommutator}, args...; order::Int=2)
    new_q, _ = _expand(q, Val(:Expand), Val(:QCommutator), order, args...)
    return new_q
end

function expand(q::QObj, ::Val{:CRational}, args...; order::Int=2)
    new_q, _ = _expand(q, Val(:Rational), Val(:CRational), order, args...)
    return new_q
end

function expand(q::QObj, ::Val{:CExp}, args...; order::Int=2)
    ord = order > 0 ? order : 2
    new_q, _ = _expand(q, Val(:Taylor), Val(:CExp), ord, args...)
    return new_q
end

function expand(q::QObj, ::Val{:CLog}, args...; order::Int=2)
    ord = order > 0 ? order : 2
    new_q, _ = _expand(q, Val(:Taylor), Val(:CLog), ord, args...)
    return new_q
end

function expand(q::QObj, ::Val{:Taylor}, ::Val{:CExp}, args...; order::Int=2)
    ord = order > 0 ? order : 2
    new_q, _ = _expand(q, Val(:Taylor), Val(:CExp), ord, args...)
    return new_q
end

function expand(q::QObj, ::Val{:Taylor}, ::Val{:CLog}, args...; order::Int=2)
    ord = order > 0 ? order : 2
    new_q, _ = _expand(q, Val(:Taylor), Val(:CLog), ord, args...)
    return new_q
end

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

@inline _zero_expr(qspace::QSpace)::QExpr = QExpr(qspace, QComposite[])

@inline function _one_expr(qspace::QSpace)::QExpr
    one_term = QAtomProduct(qspace, qspace.c_one, QAtom[QTerm(qspace.I_op)])
    return QExpr(qspace, QComposite[one_term])
end

# -----------------------------------------------------------------------------
# Default recursive traversal
# -----------------------------------------------------------------------------

_expand(q::QObj, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T} = (q, false)

function _expand(vec::Vector{T}, ::Val{M}, ::Val{Tgt}, order::Int, args...) where {M,Tgt,T<:QObj}
    any_changed = false
    new_vec = Vector{T}(undef, length(vec))
    @inbounds for i in eachindex(vec)
        new_entry, changed = _expand(vec[i], Val(M), Val(Tgt), order, args...)
        new_vec[i] = new_entry
        any_changed |= changed
    end
    return any_changed ? (new_vec, true) : (vec, false)
end

function _expand(q::QExpr, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T}
    terms2 = QComposite[]
    any_changed = false
    for term in q.terms
        term2, changed = _expand(term, Val(M), Val(T), order, args...)
        if term2 isa QExpr
            append!(terms2, term2.terms)
            any_changed = true
        else
            push!(terms2, term2)
            any_changed |= changed
        end
    end
    if any_changed
        return QExpr(q.qspace, terms2), true
    else
        return q, false
    end
end

function _expand(q::QSum, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T}
    inner, changed = _expand(q.expr, Val(M), Val(T), order, args...)
    if changed
        return modify_expr(q, inner, Val(:nodecollision)), true
    else
        return q, false
    end
end

function _expand(q::QComposite, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T}
    q_work = q
    coeff_changed = false
    if hasproperty(q, :coeff_fun) && hasmethod(modify_coeff, Tuple{typeof(q), CFunction})
        coeff = getproperty(q, :coeff_fun)
        new_coeff = _expand_coeff_fun(coeff, M, T, order)
        if new_coeff !== coeff
            q_work = modify_coeff(q_work, new_coeff)
            coeff_changed = true
        end
    end
    inner, changed = _expand(getproperty(q_work, :expr), Val(M), Val(T), order, args...)
    if changed
        return modify_expr(q_work, inner), true
    elseif coeff_changed
        return q_work, true
    else
        return q, false
    end
end

function _expand(comm::QCommutator, ::Val{:Expand}, ::Val{T}, order::Int, args...) where {T}
    @assert length(comm.expr) == 2 "QCommutator must contain exactly two expressions."
    left_expanded, _ = _expand(comm.expr[1], Val(:Expand), Val(T), order, args...)
    right_expanded, _ = _expand(comm.expr[2], Val(:Expand), Val(T), order, args...)
    coeff = _expand_coeff_fun(comm.coeff_fun, :Expand, T, order)
    expanded_comm = _comm(left_expanded, right_expanded, _NOCHK)
    expanded_comm = multiply_coeff(expanded_comm, coeff)
    return expanded_comm, true
end

function _expand(q::QMultiComposite, ::Val{M}, ::Val{T}, order::Int, args...) where {M,T}
    q_work = q
    coeff_changed = false
    if hasmethod(modify_coeff, Tuple{typeof(q), CFunction})
        coeff = q.coeff_fun
        new_coeff = _expand_coeff_fun(coeff, M, T, order)
        if new_coeff !== coeff
            q_work = modify_coeff(q_work, new_coeff)
            coeff_changed = true
        end
    end
    inner, changed = _expand(q_work.expr, Val(M), Val(T), order, args...)
    if changed
        return modify_expr(q_work, inner), true
    elseif coeff_changed
        return q_work, true
    else
        return q, false
    end
end

# -----------------------------------------------------------------------------
# Taylor expansions for specific composites
# -----------------------------------------------------------------------------

function _expand(q::QExp, ::Val{:Taylor}, ::Val{T}, order::Int=2, use_power::Bool=false) where {T}
    order < 0 && error("Taylor expansion order must be non-negative.")
    (T == :QExp || T == :Any) || return (q, false)
    ord = order > 0 ? order : 2
    inner, _ = _expand(q.expr, Val(:Taylor), Val(T), ord, use_power)
    qspace = q.qspace
    one = _one_expr(qspace)
    terms = QExpr[one]
    coeff = 1//1
    current = nothing
    for k in 1:ord
        coeff /= k
        term_expr = if use_power
            power(inner, k; force_symbolic=true)
        elseif k == 1
            inner
        else
            current * inner
        end
        current = use_power ? current : term_expr
        push!(terms, term_expr * coeff)
    end
    series = reduce(+, terms)
    coeff_target = T == :Any ? :QExp : T
    coeff_fun = _expand_coeff_fun(q.coeff_fun, :Taylor, coeff_target, ord)
    series = multiply_coeff(series, coeff_fun)
    return series, true
end

function _expand(q::QLog, ::Val{:Taylor}, ::Val{T}, order::Int=2, use_power::Bool=false) where {T}
    order < 0 && error("Taylor expansion order must be non-negative.")
    (T == :QLog || T == :Any) || return (q, false)
    ord = order > 0 ? order : 2
    ord == 0 && return (_zero_expr(q.qspace), true)
    inner, _ = _expand(q.expr, Val(:Taylor), Val(T), ord, use_power)
    qspace = q.qspace
    one = _one_expr(qspace)
    delta = inner - one
    terms = QExpr[]
    coeff = 1//1
    current = nothing
    for k in 1:ord
        if k == 1
            coeff = 1//1
        else
            coeff = -coeff * (k-1) / k
        end
        term_expr = if use_power
            power(delta, k; force_symbolic=true)
        elseif k == 1
            delta
        else
            current * delta
        end
        current = use_power ? current : term_expr
        push!(terms, term_expr * coeff)
    end
    series = isempty(terms) ? _zero_expr(qspace) : reduce(+, terms)
    coeff_target = T == :Any ? :QLog : T
    coeff_fun = _expand_coeff_fun(q.coeff_fun, :Taylor, coeff_target, ord)
    series = multiply_coeff(series, coeff_fun)
    return series, true
end

function _expand(q::QPower, ::Val{:Taylor}, ::Val{T}, order::Int=2, use_power::Bool=false) where {T}
    order > 0 || error("Taylor expansion order must be positive.")
    (T == :Any || _is_alias(T, Q_POWER_ALIASES)) || return (q, false)
    inner, _ = _expand(q.expr, Val(:Taylor), Val(T), order, use_power)
    qspace = q.qspace
    one = _one_expr(qspace)
    delta = inner - one
    terms = QExpr[one]
    coeff = 1//1
    current = nothing
    for k in 1:order
        term_expr = if use_power
            power(delta, k; force_symbolic=true)
        elseif k == 1
            delta
        else
            current * delta
        end
        current = use_power ? current : term_expr
        coeff = coeff * (q.n - (k-1)) / k
        iszero(coeff) && continue
        push!(terms, term_expr * coeff)
    end
    series = isempty(terms) ? _zero_expr(qspace) : reduce(+, terms)
    coeff_fun = _expand_coeff_fun(q.coeff_fun, :Taylor, T, order)
    series = multiply_coeff(series, coeff_fun)
    return series, true
end

function _expand(q::QRoot, ::Val{:Taylor}, ::Val{T}, order::Int=2, use_power::Bool=false) where {T}
    order > 0 || error("Taylor expansion order must be positive.")
    (T == :Any || _is_alias(T, Q_ROOT_ALIASES)) || return (q, false)
    inner, _ = _expand(q.expr, Val(:Taylor), Val(T), order, use_power)
    qspace = q.qspace
    one = _one_expr(qspace)
    delta = inner - one
    exponent = 1//q.n
    terms = QExpr[one]
    coeff = 1//1
    current = nothing
    for k in 1:order
        term_expr = if use_power
            power(delta, k; force_symbolic=true)
        elseif k == 1
            delta
        else
            current * delta
        end
        current = use_power ? current : term_expr
        coeff = coeff * (exponent - (k-1)) / k
        iszero(coeff) && continue
        push!(terms, term_expr * coeff)
    end
    series = isempty(terms) ? _zero_expr(qspace) : reduce(+, terms)
    coeff_fun = _expand_coeff_fun(q.coeff_fun, :Taylor, T, order)
    series = multiply_coeff(series, coeff_fun)
    return series, true
end
