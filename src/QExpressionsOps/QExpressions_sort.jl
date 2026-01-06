import Base: isless, sort, sort!

# -------------------------------
# Small helpers (boolean style)
# -------------------------------

@inline function less_vec_int(a::AbstractVector{<:Int}, b::AbstractVector{<:Int})
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        ai = a[i]
        bi = b[i]
        if ai != bi
            return ai < bi
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
# Particle comparisons
# -------------------------------

@inline function Base.isequal(a::QParticle, b::QParticle)
    return a.operator == b.operator && isequal(a.index, b.index)
end

@inline function Base.:(==)(a::QParticle, b::QParticle)
    return a.operator == b.operator && a.index == b.index
end

@inline function Base.isless(a::QParticle, b::QParticle)
    if a.operator != b.operator
        return less_vec_int(a.operator, b.operator)
    end
    return isless(a.index, b.index)
end

# -------------------------------
# Tags for cross-type ordering
# -------------------------------
qobj_tag(::QAtomProduct)      = 0
qobj_tag(::QSum)              = 1
qobj_tag(::QInt)              = 2
qobj_tag(::QCompositeProduct) = 3
qobj_tag(::QExp)              = 4
qobj_tag(::QLog)              = 5
qobj_tag(::QCommutator)       = 6
qobj_tag(::QPower)            = 7
qobj_tag(::QRoot)             = 8

qatom_tag(::QTerm)     = 0
qatom_tag(::QAbstract) = 1

@inline aggregator_order_key(::Type{SumAggregator}) = 0
@inline aggregator_order_key(::Type{IntegralAggregator}) = 1
@inline aggregator_order_key(::Type{T}) where {T<:AbstractQAggregator} = 10

# -------------------------------
# Atom-level isless
# -------------------------------
# QTerm: time index first, then length, then particles lexicographically
function isless(a::QTerm, b::QTerm)::Bool
    a.time_index == b.time_index || return a.time_index < b.time_index
    la = length(a.op_indices)
    lb = length(b.op_indices)
    la != lb && return la < lb
    @inbounds for i in 1:la
        ai = a.op_indices[i]
        bi = b.op_indices[i]
        if !isequal(ai, bi)
            return isless(ai, bi)
        end
    end
    return false
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

# AbstractQSum: compare block metadata and inner expression
@inline _total_indices(q::AbstractQSum) = sum(length(block.indices) for block in q.blocks)

function isless_same(a::AbstractQSum, b::AbstractQSum)
    ta = aggregator_type(a)
    tb = aggregator_type(b)
    key_a = aggregator_order_key(ta)
    key_b = aggregator_order_key(tb)
    if key_a != key_b
        return key_a < key_b
    elseif ta != tb
        return String(nameof(ta)) < String(nameof(tb))
    end
    na_idx = _total_indices(a)
    nb_idx = _total_indices(b)
    _total_indices(a) != _total_indices(b) && return na_idx < nb_idx

    @inbounds for (blk_a, blk_b) in zip(a.blocks, b.blocks)
        len_a = length(blk_a.indices); len_b = length(blk_b.indices)
        if len_a != len_b
            return len_a < len_b
        end

        @inbounds for (idx_a, idx_b) in zip(blk_a.indices, blk_b.indices)
            if idx_a != idx_b 
                return idx_a < idx_b 
            end
        end

        @inbounds for (row_a, row_b) in zip(blk_a.constraints, blk_b.constraints)
            row_a == row_b || return row_a < row_b 
        end
    end
    return isless(a.expr, b.expr)
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
