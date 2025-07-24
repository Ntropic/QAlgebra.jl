import Base: exp, log

export QCommutator, QExp, QLog, QPower, power, QRoot, root

mutable struct QCommutator <: QMultiComposite
    statespace::StateSpace
    expr::Vector{qExpr}
end
function QCommutator(q1::qExpr, q2::qExpr)::qExpr
    return qExpr(QCommutator(q1.statespace, qExpr[q1, q2]))
end
copy(q::QCommutator) = QCommutator(q.statespace, [copy(x) for x in q.expr])

""" 
    QCompositeProduct

Represents a product of qComposites. 
"""
mutable struct QCompositeProduct <: QMultiComposite
    statespace::StateSpace         # State space of the product.
    expr::Vector{qExpr} 
end
function QCompositeProduct(expr::Vector{qExpr} )
    new(expr[1].statespace, expr)
end
copy(q::QCompositeProduct) = QCompositeProduct(q.statespace, [copy(x) for x in q.expr])

### Non-simple qFunctions 
mutable struct QExp <: QComposite
    statespace::StateSpace
    expr::qExpr
end
function exp(q::qExpr)::qExpr
    return qExpr(QExp(q.statespace, q))
end
copy(q::QExp) = QExp(q.statespace, copy(q.expr))


mutable struct QLog <: QComposite
    statespace::StateSpace
    expr::qExpr
end
function log(q::qExpr)::qExpr
    return qExpr(QLog(q.statespace, q))
end
copy(q::QLog) = QLog(q.statespace, copy(q.expr))


mutable struct QPower <: QComposite
    statespace::StateSpace
    n::Int
    expr::qExpr
end
copy(q::QPower) = QPower(q.statespace, q.n, copy(q.expr))

function power(q::qExpr, n::Int)::qExpr
    if n == 1 
         return q 
    end
    return qExpr(QPower(q.statespace, n, q))
end


mutable struct QRoot <: QComposite
    statespace::StateSpace
    n::Int
    expr::qExpr
end
function root(q::qExpr, n::Int=2)::qExpr
    if n == 1 
        return q 
    end
    return qExpr(QRoot(q.statespace, n, q))
end
copy(q::QRoot) = QRoot(q.statespace, q.n, copy(q.expr))