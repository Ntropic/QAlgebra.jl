# ---------- helpers ------------------------------------------------------------

@inline function _find_abs_position(abs_params::Vector{AbstractCAbstract}, idx::Int)::Int
    @inbounds for j in eachindex(abs_params)
        # stored type is AbstractCAbstract, but it's concretely CAbstract, to not have to check we specify
        if (abs_params[j]::CAbstract).index == idx
            return j
        end
    end
    return 0
end

"""
    substitute_and_simplify(node, abs_params, values) -> (newnode, changed::Bool)

Replace every `CAbstract` in `node` that matches (by `.index`) one
of `abs_params` with the corresponding `values[j]`. On any change,
rebuild using outer constructors (which call simplify) before returning.
- `abs_params` and `values` must have the same length and order.
"""
function substitute_and_simplify(f::CFunction,
                                 abs_params::Vector{AbstractCAbstract},
                                 values::Vector{CFunction})::Tuple{CFunction,Bool}
    @assert length(abs_params) == length(values)
    _sub(f, abs_params, values)
end


# CAtom: nothing to do
_sub(a::CAtom,    ::Vector{AbstractCAbstract}, ::Vector{CFunction}) = (a, false)
function _sub(a::CAbstract, abs_params::Vector{AbstractCAbstract}, values::Vector{CFunction})
    j = _find_abs_position(abs_params, a.index)
    return j == 0 ? (a, false) : (values[j], true)
end

# Each unary composite: recurse on .expr; if changed, rebuild via outer ctor (simplifies)
function _sub(e::CExp, abs_params, values)
    newchild, ch = _sub(e.expr, abs_params, values)
    ch || return (e, false)
    # outer constructor calls simplify_CExp
    return (CExp(e.param_info, e.coeff, newchild), true)
end

function _sub(l::CLog, abs_params, values)
    newchild, ch = _sub(l.expr, abs_params, values)
    ch || return (l, false)
    return (CLog(l.param_info, l.coeff, newchild), true)
end

function _sub(p::CPower, abs_params, values)
    newchild, ch = _sub(p.expr, abs_params, values)
    ch || return (p, false)
    return (CPower(p.param_info, p.coeff, newchild, p.exponent), true)
end

# ---------- multi composites ---------------------------------------------------

# CSum: map children; on change, rebuild via _CSum (simplifies) 
function _sub(s::CSum, abs_params, values)
    any_changed = false
    newkids = similar(s.expr)
    @inbounds for i in eachindex(s.expr)
        newkids[i], ch = _sub(s.expr[i], abs_params, values)
        any_changed |= ch
    end
    any_changed || return (s, false)
    return (_CSum(s.param_info, newkids), true)  # simplifies/flat/possibly returns single child
end

# CProd: map children; on change, rebuild via outer CProd (simplifies)
function _sub(p::CProd, abs_params, values)
    any_changed = false
    newkids = similar(p.expr)
    @inbounds for i in eachindex(p.expr)
        newkids[i], ch = _sub(p.expr[i], abs_params, values)
        any_changed |= ch
    end
    any_changed || return (p, false)
    return (CProd(p.param_info, p.coeff, newkids), true)
end

# CCustomType: map its arguments; if any changed, rebuild node (no auto-simplify here)
function _sub(c::CCustomType, abs_params, values)
    any_changed = false
    newkids = similar(c.expr)
    @inbounds for i in eachindex(c.expr)
        newkids[i], ch = _sub(c.expr[i], abs_params, values)
        any_changed |= ch
    end
    any_changed || return (c, false)
    # Rebuild with same def; if you want extra rules, hook them here
    return (CCustomType(c.param_info, c.coeff, newkids, c.ctype_def), true)
end

# ---------- special binary -----------------------------------------------------

function _sub(r::CRational, abs_params, values)
    newN, chN = _sub(r.numer, abs_params, values)
    newD, chD = _sub(r.denom, abs_params, values)
    (chN || chD) || return (r, false)
    return (CRational(r.param_info, newN, newD), true)  # outer ctor simplifies
end

# ---------- containers ---------------------------------------------------------

function _sub(v::CVector, abs_params, values)
    any_changed = false
    newkids = similar(v.entries)
    @inbounds for i in eachindex(v.entries)
        newkids[i], ch = _sub(v.entries[i], abs_params, values)
        any_changed |= ch
    end
    any_changed || return (v, false)
    return (CVector(v.param_info, v.coeff, newkids; row=v.row), true)
end

function _sub(M::CMatrix, abs_params, values)
    any_changed = false
    entries = M.entries
    flat = Vector{CFunction}(undef, length(entries))
    # linear index over entries
    k = 1
    for j in axes(entries, 2), i in axes(entries, 1)
        flat[k], ch = _sub(entries[i,j], abs_params, values)
        any_changed |= ch
        k += 1
    end
    any_changed || return (M, false)
    return (CMatrix(M.param_info, M.coeff, reshape(flat, size(entries))), true)
end
