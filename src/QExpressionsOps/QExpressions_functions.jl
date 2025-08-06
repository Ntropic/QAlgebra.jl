import Base: exp, log, sqrt

export QCommutator, QExp, QLog, QPower, power, QRoot, root

mutable struct QCommutator <: QMultiComposite
    statespace::StateSpace
    coeff_fun::CFunction
    expr::Vector{QExpr}
end
function QCommutator(q1::QExpr, q2::QExpr)::QExpr
    return QExpr(QCommutator(q1.statespace, q1.statespace.fone, QExpr[q1, q2]))
end
copy(q::QCommutator) = QCommutator(q.statespace, copy(q.coeff_fun), [copy(x) for x in q.expr])

""" 
    QCompositeProduct

Represents a product of QComposites. 
"""
mutable struct QCompositeProduct <: QMultiComposite
    statespace::StateSpace         # State space of the product.
    coeff_fun::CFunction
    expr::Vector{QExpr} 
    function QCompositeProduct(statespace::StateSpace, expr::Vector{QExpr} )
        new(statespace, statespace.fone, expr)
    end
    function QCompositeProduct(expr::Vector{QExpr} )
        new(expr[1].statespace, expr[1].statespace.fone, expr)
    end
    function QCompositeProduct(expr::Vector{QComposite} )
        new(expr[1].statespace, expr[1].statespace.fone, [QExpr([x]) for x in expr])
    end
end
copy(q::QCompositeProduct) = QCompositeProduct(q.statespace, copy(q.coeff_fun), [copy(x) for x in q.expr])

### Non-simple qFunctions 
mutable struct QExp <: QComposite
    statespace::StateSpace
    coeff_fun::CFunction
    expr::QExpr
end
function exp(q::QExpr)::QExpr
    return QExpr([QExp(q.statespace, q.statespace.fone, q)])s
end
copy(q::QExp) = QExp(q.statespace, copy(q.coeff_fun), copy(q.expr))
iszero(q::QExp) = iszero(q.coeff_fun) 
"""function simplify_rules(q::QExp)::Vector{QComposite}
    if iszero(q)
        return [QAtomProduct(q.statespace, q.statespace.fone, QTerm(q.statespace.neutral_op))]
    end
    if length(q.expr) == 1
        if isa(q.expr[1], QLog)
            return # Continue here
    return [q] 
end"""

mutable struct QLog <: QComposite
    statespace::StateSpace
    coeff_fun::CFunction
    expr::QExpr
end
function log(q::QExpr)::QExpr
    return QExpr([QLog(q.statespace, q.statespace.fone, q)])
end
copy(q::QLog) = QLog(q.statespace, copy(q.coeff_fun), copy(q.expr))
iszero(q::QLog) = iszero(q.coeff_fun) || isone(q.expr)

mutable struct QPower <: QComposite
    statespace::StateSpace
    coeff_fun::CFunction
    n::Int
    expr::QExpr
end
copy(q::QPower) = QPower(q.statespace, copy(q.coeff_fun), q.n, copy(q.expr))

""" 
    power(q::QExpr, n::Int)::QExpr

Returns the q^n of a quantum expression, without expanding the expression.
"""
function power(q::QExpr, n::Int)::QExpr
    if n == 1 
         return q 
    end
    return QExpr([QPower(q.statespace, q.statespace.fone, n, q)])
end


mutable struct QRoot <: QComposite
    statespace::StateSpace
    coeff_fun::CFunction
    n::Int
    expr::QExpr
end

""" 
    root(q::QExpr, n::Int)::QExpr

Returns the n'th root of a quantum expression q^{1/n}, without expanding the expression.
"""
function root(q::QExpr, n::Int=2)::QExpr
    if n == 1 
        return q 
    end
    return QExpr([QRoot(q.statespace, q.statespace.fone, n, q)])
end
function sqrt(q::QExpr)::QExpr
    return QExpr([QRoot(q.statespace, q.statespace.fone, 2, q)])
end
copy(q::QRoot) = QRoot(q.statespace, copy(q.coeff_fun), q.n, copy(q.expr))
