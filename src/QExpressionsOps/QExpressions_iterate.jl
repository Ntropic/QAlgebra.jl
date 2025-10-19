export tree_iter_composite, modify_coeff_funs_tree_composite
const OrderedTypes = Union{QCumulantOrdered, QAtomOrdered, QAtomIndexed, QNeutral} 
"""
    tree_iter_composite(obj::QObj)

Depth-first, pre-order iterator over every `QComposite` nested inside `obj`.
The iterator yields each composite before its children.
"""
tree_iter_composite(q::diffQEq) = Iterators.flatten((tree_iter_composite(q.left_hand_side), tree_iter_composite(q.expr)))
tree_iter_composite(::diffQEqOrdered) = error("tree_iter_composite is not defined for diffQEqOrdered; convert to diffQEq first.")

function tree_iter_composite(q::QExpr)
    return Iterators.flatten(tree_iter_composite(term) for term in q.terms)
end
tree_iter_composite(q::QAtomProduct) = (q,)
tree_iter_composite(::T) where T<:OrderedTypes = error("tree_iter_composite is not defined for $T or other ordered types; convert to a supported type first.")

function tree_iter_composite(q::T) where T<:QMultiComposite
    return Iterators.flatten(((q,), (tree_iter_composite(child) for child in q.expr)))
end

function tree_iter_composite(q::T) where T<:QComposite
    hasproperty(q, :expr) || error("tree_iter_composite requires composites with an `expr` field; got $(T).")
    return Iterators.flatten(((q,), tree_iter_composite(q.expr)))
end


# ==============================================> Reconstructing Terms <==============================================
@inline function _next_coeff!(coeffs::Vector{CFunction}, pos::Base.RefValue{Int})::CFunction
    idx = pos[]
    if idx > length(coeffs)
        throw(ArgumentError("Not enough coefficient functions provided; expected at least $(idx) but only got $(length(coeffs))."))
    end
    pos[] = idx + 1
    return coeffs[idx]
end

function _modify_coeff_funs_tree_composite(q::QExpr, coeffs::Vector{CFunction}, pos::Base.RefValue{Int})::QExpr
    new_terms = Vector{QComposite}(undef, length(q.terms))
    @inbounds for i in eachindex(q.terms)
        new_terms[i] = _modify_coeff_funs_tree_composite(q.terms[i], coeffs, pos)
    end
    return QExpr(q.qspace, new_terms, Val(:nosimp))
end

function _modify_coeff_funs_tree_composite(eq::diffQEq, coeffs::Vector{CFunction}, pos::Base.RefValue{Int})::diffQEq
    lhs = _modify_coeff_funs_tree_composite(eq.left_hand_side, coeffs, pos)
    rhs = _modify_coeff_funs_tree_composite(eq.expr, coeffs, pos)
    return diffQEq(eq.qspace, lhs, rhs, Val(:raw))
end

_modify_coeff_funs_tree_composite(::T, ::Vector{CFunction}, ::Base.RefValue{Int}) where T<:OrderedTypes =  error("modify_coeff_funs_tree_composite is not implemented for $T or other Ordered types.")

function _modify_coeff_funs_tree_composite(q::T, coeffs::Vector{CFunction}, pos::Base.RefValue{Int}) where T<:QMultiComposite
    new_children = similar(q.expr)
    @inbounds for i in eachindex(q.expr)
        new_children[i] = _modify_coeff_funs_tree_composite(q.expr[i], coeffs, pos)
    end
    coeff = _next_coeff!(coeffs, pos)
    return modify_coeff_expr(q, coeff, new_children)
end

function _modify_coeff_funs_tree_composite(q::T, coeffs::Vector{CFunction}, pos::Base.RefValue{Int}) where T<:QComposite
    hasproperty(q, :expr) || error("modify_coeff_funs_tree_composite requires composites with an `expr` field; got $(T).")
    expr_field = q.expr

    new_expr = _modify_coeff_funs_tree_composite(expr_field, coeffs, pos)

    coeff = hasproperty(q, :coeff_fun) ? _next_coeff!(coeffs, pos) : q.coeff_fun
    return modify_coeff_expr(q, coeff, new_expr)
end

function _modify_coeff_funs_tree_composite(::QCumulantOrdered, ::Vector{CFunction}, ::Base.RefValue{Int})
    error("modify_coeff_funs_tree_composite is not implemented for QCumulantOrdered; convert to a supported type first.")
end

"""
    modify_coeff_funs_tree_composite(obj::QObj, coeffs::Vector{CFunction})

Return a copy of `obj` where each `QComposite` that stores a `coeff_fun`
is replaced by one using coefficient functions from `coeffs`.
Coefficients are assigned in the same order that `tree_iter_composite(obj)`
would visit the composites that expose a `coeff_fun` field.

Throws an `ArgumentError` if the number of provided coefficients does
not match the number of such composites.
"""
function modify_coeff_funs_tree_composite(obj::QObj, coeffs::Vector{CFunction})
    pos = Base.RefValue{Int}(1)
    new_obj = _modify_coeff_funs_tree_composite(obj, coeffs, pos)
    used = pos[] - 1
    if used != length(coeffs)
        throw(ArgumentError("Too many coefficient functions provided; expected $(used) but got $(length(coeffs))."))
    end
    return new_obj
end
