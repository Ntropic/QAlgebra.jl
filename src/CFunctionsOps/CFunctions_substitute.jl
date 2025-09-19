"""
    substitute(node, abstract, value) -> (newnode, changed::Bool)

Replace every `AbstractCAbstract` in `node` that matches (by `.index`) one
of `abstract` with the corresponding `value[j]`. On any change,
rebuild using outer constructors (which call simplify) before returning.
- `abstract` and `value` must have the same length and order.
"""
function substitute(f::CFunction, abstract::AbstractCAbstract, value::CFunction)::CFunction
    return _sub(f, abstract, value)[1]
end
function substitute(f::CFunction, abstracts::Vector{AbstractCAbstract}, values::AbstractVector{<:CFunction})::CFunction
    length(abstracts) <= length(values) || throw(ArgumentError("need at least as many values as abstracts"))
    isempty(abstracts) && return f
    node = f
    for (abstract, value) in zip(abstracts, values)
        node, _ = _sub(node, abstract, value)
    end
    return node
end

# CAtom: nothing to do
_sub(a::CAtom,  ::AbstractCAbstract, ::CFunction) = (a, false)
function _sub(a::CAbstract, abstract::AbstractCAbstract, value::CFunction)
    if a.index == abstract.index
        new_value = a.dag ? conj(value) : value 
        if a.exponent != 1
            new_value = new_value^a.exponent
        end
        new_value *= a.coeff
        return (new_value, true)
    end
    return (a, false)
end

# Each unary composite: recurse on .expr; if changed, rebuild via outer ctor (simplifies)
function _sub(e::CExp, abstract::AbstractCAbstract, value::CFunction)
    newchild, ch = _sub(e.expr, abstract, value)
    ch || return (e, false)
    # outer constructor calls simplify_CExp
    return (CExp(e.param_info, e.coeff, newchild), true)
end

function _sub(l::CLog, abstract::AbstractCAbstract, value::CFunction)
    newchild, ch = _sub(l.expr, abstract, value)
    ch || return (l, false)
    return (CLog(l.param_info, l.coeff, newchild), true)
end

function _sub(p::CPower, abstract::AbstractCAbstract, value::CFunction)
    newchild, ch = _sub(p.expr, abstract, value)
    ch || return (p, false)
    return (CPower(p.param_info, p.coeff, newchild, p.exponent), true)
end

# ---------- multi composites ---------------------------------------------------

# CSum: map children; on change, rebuild via _CSum (simplifies) 
function _sub(s::CSum, abstract::AbstractCAbstract, value::CFunction)
    any_changed = false
    newkids = similar(s.expr)
    @inbounds for i in eachindex(s.expr)
        newkids[i], ch = _sub(s.expr[i], abstract, value)
        any_changed |= ch
    end
    any_changed || return (s, false)
    return (_CSum(s.param_info, newkids), true)  # simplifies/flat/possibly returns single child
end

# CProd: map children; on change, rebuild via outer CProd (simplifies)
function _sub(p::CProd, abstract::AbstractCAbstract, value::CFunction)
    any_changed = false
    newkids = similar(p.expr)
    @inbounds for i in eachindex(p.expr)
        newkids[i], ch = _sub(p.expr[i], abstract, value)
        any_changed |= ch
    end
    any_changed || return (p, false)
    return (CProd(p.param_info, p.coeff, newkids), true)
end

# CCustomType: map its arguments; if any changed, rebuild node (no auto-simplify here)
function _sub(c::CCustomType, abstract::AbstractCAbstract, value::CFunction)
    any_changed = false
    newkids = similar(c.expr)
    @inbounds for i in eachindex(c.expr)
        newkids[i], ch = _sub(c.expr[i], abstract, value)
        any_changed |= ch
    end
    any_changed || return (c, false)
    # Rebuild with same def; if you want extra rules, hook them here
    return (CCustomType(c.param_info, c.coeff, newkids, c.ctype_def), true)
end

# ---------- special binary -----------------------------------------------------

function _sub(r::CRational, abstract::AbstractCAbstract, value::CFunction)
    newN, chN = _sub(r.numer, abstract, value)
    newD, chD = _sub(r.denom, abstract, value)
    (chN || chD) || return (r, false)
    return (CRational(r.param_info, newN, newD), true)  # outer ctor simplifies
end

# ---------- containers ---------------------------------------------------------

function _sub(v::CVector, abstract::AbstractCAbstract, value::CFunction)
    any_changed = false
    newkids = similar(v.entries)
    @inbounds for i in eachindex(v.entries)
        newkids[i], ch = _sub(v.entries[i], abstract, value)
        any_changed |= ch
    end
    any_changed || return (v, false)
    return (CVector(v.param_info, v.coeff, newkids; row=v.row), true)
end

function _sub(M::CMatrix, abstract::AbstractCAbstract, value::CFunction)
    any_changed = false
    entries = M.entries
    flat = Vector{CFunction}(undef, length(entries))
    # linear index over entries
    k = 1
    for j in axes(entries, 2), i in axes(entries, 1)
        flat[k], ch = _sub(entries[i,j], abstract, value)
        any_changed |= ch
        k += 1
    end
    any_changed || return (M, false)
    return (CMatrix(M.param_info, M.coeff, reshape(flat, size(entries))), true)
end
