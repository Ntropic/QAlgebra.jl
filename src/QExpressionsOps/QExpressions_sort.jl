import Base: isless, sort, sort!

# -------------------------------
# Small helpers (boolean style)
# -------------------------------

@inline function less_vec_int(a::AbstractVector{<:Integer}, b::AbstractVector{<:Integer})
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        ai = a[i]; bi = b[i]
        if ai != bi
            return ai < bi
        end
    end
    return length(a) < length(b)
end

@inline function less_vec(a::AbstractVector, b::AbstractVector)
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        ai = a[i]; bi = b[i]
        if !(ai == bi)   # cheap eq; assumes well-behaved == for your types
            return isless(ai, bi)
        end
    end
    return length(a) < length(b)
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
isless(a::QTerm, b::QTerm) = less_vec_int(a.op_indices, b.op_indices)

# QAbstract: compare (key_index, sub_index, exponent, dag)
function isless(a::QAbstract, b::QAbstract)
    a.key_index != b.key_index && return a.key_index < b.key_index
    a.sub_index != b.sub_index && return a.sub_index < b.sub_index
    a.exponent  != b.exponent  && return a.exponent  < b.exponent
    dag_a = Int(a.dag); dag_b = Int(b.dag)
    dag_a != dag_b && return dag_a < dag_b
    return false
end

# Cross-type atoms
function isless(a::Union{QTerm,QAbstract}, b::Union{QTerm,QAbstract})
    ta = qatom_tag(a); tb = qatom_tag(b)
    ta != tb && return ta < tb
    return a isa QTerm ? isless(a::QTerm, b::QTerm) :
                         isless(a::QAbstract, b::QAbstract)
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
    a.subsystem_index != b.subsystem_index && return a.subsystem_index < b.subsystem_index
    la = length(a.element_indexes); lb = length(b.element_indexes)
    la != lb && return la < lb
    if !(a.element_indexes == b.element_indexes)
        return less_vec_int(a.element_indexes, b.element_indexes)
    end
    if !(a.expr == b.expr)
        return isless(a.expr, b.expr)
    end
    if !(a.neq == b.neq)
        return isless(a.neq, b.neq)   # assumes isless for its type
    end
    return false
end

# QAtomProduct: coeff first, then atoms (length + pairwise)
function isless_same(a::QAtomProduct, b::QAtomProduct)
    if !(a.coeff_fun == b.coeff_fun)
        return compare_isless(a.coeff_fun, b.coeff_fun)  # relies on CFunction.isless
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
    return QExpr(qeq.statespace, terms)
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
    return modify_expr(q, sort.(q.expr; kwargs...))
end
function sort!(q::T; kwargs...) where {T<:QMultiComposite}
    modify_expr!(q, sort!.(q.expr; kwargs...))
    return q
end

# QComposite — sort inner expr only
function sort(q::T; kwargs...) where {T<:QComposite}
    return modify_expr(q, sort(q.expr; kwargs...))
end
function sort!(q::T; kwargs...) where {T<:QComposite}
    modify_expr!(q, sort!(q.expr; kwargs...))
    return q
end
