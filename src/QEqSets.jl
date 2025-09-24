module QEqSets

using ..QExpressions: diffQEq, Expectation, QObj, QComposite, QAtomProduct, QExpr

export diffQEqSet, ExpectQAtomProducts

"""
    ExpectQAtomProducts(obj)

Return a copy of `obj` where every `QAtomProduct` is wrapped in an expectation value.
Works on individual products, composite expressions, differential equations, and
collections thereof. Non-composite quantum objects are returned unchanged.
"""
ExpectQAtomProducts(expr::QExpr) = Expectation(expr)
ExpectQAtomProducts(prod::QAtomProduct) = Expectation(prod)
ExpectQAtomProducts(comp::QComposite) = Expectation(comp)
ExpectQAtomProducts(obj::QObj) = obj
ExpectQAtomProducts(arr::AbstractVector{<:QObj}) = ExpectQAtomProducts.(arr)
function ExpectQAtomProducts(eq::diffQEq)
    lhs = ExpectQAtomProducts(eq.left_hand_side)
    rhs = ExpectQAtomProducts(eq.expr)
    return diffQEq(eq.qspace, lhs, rhs, Val(:raw))
end

struct diffQEqSet
    equations::Vector{diffQEq}
    loss::Function
    function diffQEqSet(equations::Vector{diffQEq}, loss::Function, ::Val{:raw})
        new(equations, loss)
    end
end

function diffQEqSet(; kwargs...)
    @warn "diffQEqSet construction is not implemented yet; returning an empty set."
    return diffQEqSet(diffQEq[], () -> nothing, Val(:raw))
end

function diffQEqSet(equations::Vector{diffQEq}; loss::Function = () -> nothing)
    return diffQEqSet(equations, loss, Val(:raw))
end

function diffQEqSet(equations::Vector{diffQEq}, loss::Function)
    return diffQEqSet(equations, loss, Val(:raw))
end

function ExpectQAtomProducts(eqset::diffQEqSet)
    equations = ExpectQAtomProducts(eqset.equations)
    return diffQEqSet(equations, eqset.loss, Val(:raw))
end

ExpectedValues(obj) = ExpectQAtomProducts(obj)

end
