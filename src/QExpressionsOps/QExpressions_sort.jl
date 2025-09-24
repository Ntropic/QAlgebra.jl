import Base: isless, sort, sort!

# -------------------------------
# Small helpers (boolean style)
# -------------------------------

@inline function less_vec_int(a::Vector{Vector{Int}}, b::Vector{Vector{Int}})
    na = length(a)
    nb = length(b)
    n = min(na, nb)
    @inbounds for offset in 0:n-1
        ai = a[na - offset]
        bi = b[nb - offset]
        if ai != bi
            return ai < bi
        end
    end
    return na < nb
end

@inline function compare_isless(x, y)
    isless(x, y) && return true
    isless(y, x) && return false
    return false  # equal
end

# -------------------------------
# Tags for cross-type ordering
# -------------------------------
qobj_tag(::QAtomProduct)      = 0
qobj_tag(::QSum)              = 1
qobj_tag(::QCompositeProduct) = 2
qobj_tag(::QExp)              = 3
qobj_tag(::QLog)              = 4
qobj_tag(::QCommutator)       = 5
qobj_tag(::QPower)            = 6
qobj_tag(::QRoot)             = 7

qatom_tag(::QTerm)     = 0
qatom_tag(::QAbstract) = 1

# -------------------------------
# Atom-level isless
# -------------------------------
# QTerm: compare op_indices lexicographically
function isless(a::QTerm, b::QTerm)::Bool
    a.time_index == b.time_index || return a.time_index < b.time_index
    return less_vec_int(a.op_indices, b.op_indices)
end

# QAbstract: compare (key_index, sub_index, exponent, dag)
function isless(a::QAbstract, b::QAbstract)
    return (a.time_index, a.key_index, a.sub_index, a.dag) < (b.time_index, b.key_index, b.sub_index, b.dag)
end

# Cross-type atoms
function isless(a::Union{QTerm,QAbstract}, b::Union{QTerm,QAbstract})
    return qatom_tag(a) < qatom_tag(b)
end

# -------------------------------
# QExpr isless: by length, then pairwise terms
# -------------------------------
function isless(a::QExpr, b::QExpr)
    la = length(a); lb = length(b)
    la != lb && return la < lb
    n = la  # == lb
    @inbounds for i in 1:n
        ai = a.terms[i]; bi = b.terms[i]
        if !(ai == bi)
            return isless(ai, bi)  # QComposite
        end
    end
    return false
end

# -------------------------------
# QComposite isless (cross-type tag, then same-kind)
# -------------------------------
function isless(a::QComposite, b::QComposite)
    ta = qobj_tag(a); tb = qobj_tag(b)
    ta != tb && return ta < tb
    return isless_same(a, b)
end

# Fallback for simple composites: compare inner exprs
isless_same(a::QComposite, b::QComposite) = isless(a.expr, b.expr)

# QSum: subsystem_index, |element_indexes|, element_indexes, expr, neq
function isless_same(a::QSum, b::QSum)
    # total number of indexes
    na_ind, nb_ind = length(a.eq_indexes), length(b.eq_indexes)
    if na_ind != nb_ind
        return na_ind < nb_ind
    end
    na_ind, nb_ind = length(a.neq_blocks), length(b.neq_blocks)
    if na_ind != nb_ind
        return na_ind < nb_ind
    end

    for (blk_a, blk_b) in zip(a.neq_blocks, b.neq_blocks)
        la, lb = length(blk_a), length(blk_b)
        if la != lb
            return la < lb
        end
        if blk_a != blk_b
            return blk_a < blk_b
        end
    end
    return a.expr < b.expr
end

# QAtomProduct: coeff first, then atoms (length + pairwise)
function isless_same(a::QAtomProduct, b::QAtomProduct)
    if !(a.coeff_fun == b.coeff_fun)
        lt = compare_isless(a.coeff_fun, b.coeff_fun)
        gt = compare_isless(b.coeff_fun, a.coeff_fun)
        lt != gt && return lt
    end
    la = length(a.expr); lb = length(b.expr)
    la != lb && return la < lb
    n = la
    @inbounds for i in 1:n
        ai = a.expr[i]; bi = b.expr[i]  # QTerm/QAbstract
        if !(ai == bi)
            return isless(ai, bi)
        end
    end
    return false
end

# QCompositeProduct: length + pairwise child QExprs
function isless_same(a::QCompositeProduct, b::QCompositeProduct)
    la = length(a.expr); lb = length(b.expr)
    la != lb && return la < lb
    n = la
    @inbounds for i in 1:n
        ai = a.expr[i]; bi = b.expr[i]  # QExpr
        if !(ai == bi)
            return isless(ai, bi)
        end
    end
    return false
end

# QPower / QRoot
function isless_same(a::QPower, b::QPower)
    a.n != b.n && return a.n < b.n
    return isless(a.expr, b.expr)
end
function isless_same(a::QRoot, b::QRoot)
    a.n != b.n && return a.n < b.n
    return isless(a.expr, b.expr)
end

# -------------------------------
# Sort APIs
# -------------------------------

# QExpr — new sorted value
function sort(qeq::QExpr; kwargs...)
    terms = copy(qeq.terms)
    sort!(terms; kwargs...)             # uses isless(::QComposite)
    return QExpr(qeq.qspace, terms)
end

# QExpr — in place
function sort!(qeq::QExpr; kwargs...)
    sort!(qeq.terms; kwargs...)
    return qeq
end

# QAtomProduct — leave outer as-is
sort(qprod::QAtomProduct; kwargs...)  = copy(qprod)
sort!(qprod::QAtomProduct; kwargs...) = qprod

# QMultiComposite — recursively sort children, not the outer
function sort(q::T; kwargs...) where {T<:QMultiComposite}
    sorted_children = map(only, sort.(q.expr; kwargs...))
    return only(modify_expr(q, sorted_children))
end
function sort!(q::T; kwargs...) where {T<:QMultiComposite}
    modify_expr!(q, sort!.(q.expr; kwargs...))
    return q
end

# QComposite — sort inner expr only
function sort(q::T; kwargs...) where {T<:QComposite}
    return only(modify_expr(q, sort(q.expr; kwargs...)))
end
function sort!(q::T; kwargs...) where {T<:QComposite}
    modify_expr!(q, sort!(q.expr; kwargs...))
    return q
end
