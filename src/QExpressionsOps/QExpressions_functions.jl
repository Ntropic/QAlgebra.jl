import Base: exp, log

export QCommutator, QExp, QLog, QPower, power, QRoot, root

mutable struct QCommutator <: QMultiComposite
    statespace::StateSpace
    expr::Vector{QExpr}
end
function QCommutator(q1::QExpr, q2::QExpr)::QExpr
    return QExpr(QCommutator(q1.statespace, QExpr[q1, q2]))
end
copy(q::QCommutator) = QCommutator(q.statespace, [copy(x) for x in q.expr])

""" 
    QCompositeProduct

Represents a product of QComposites. 
"""
mutable struct QCompositeProduct <: QMultiComposite
    statespace::StateSpace         # State space of the product.
    expr::Vector{QExpr} 
    function QCompositeProduct(statespace::StateSpace, expr::Vector{QExpr} )
        new(statespace, expr)
    end
    function QCompositeProduct(expr::Vector{QExpr} )
        new(expr[1].statespace, expr)
    end
    function QCompositeProduct(expr::Vector{QComposite} )
        new(expr[1].statespace, [QExpr([x]) for x in expr])
    end
end
copy(q::QCompositeProduct) = QCompositeProduct(q.statespace, [copy(x) for x in q.expr])

### Non-simple qFunctions 
mutable struct QExp <: QComposite
    statespace::StateSpace
    expr::QExpr
end
function exp(q::QExpr)::QExpr
    return QExpr(QExp(q.statespace, q))
end
copy(q::QExp) = QExp(q.statespace, copy(q.expr))


mutable struct QLog <: QComposite
    statespace::StateSpace
    expr::QExpr
end
function log(q::QExpr)::QExpr
    return QExpr([QLog(q.statespace, q)])
end
copy(q::QLog) = QLog(q.statespace, copy(q.expr))


mutable struct QPower <: QComposite
    statespace::StateSpace
    n::Int
    expr::QExpr
end
copy(q::QPower) = QPower(q.statespace, q.n, copy(q.expr))

function power(q::QExpr, n::Int)::QExpr
    if n == 1 
         return q 
    end
    return QExpr(QPower(q.statespace, n, q))
end


mutable struct QRoot <: QComposite
    statespace::StateSpace
    n::Int
    expr::QExpr
end
function root(q::QExpr, n::Int=2)::QExpr
    if n == 1 
        return q 
    end
    return QExpr(QRoot(q.statespace, n, q))
end
copy(q::QRoot) = QRoot(q.statespace, q.n, copy(q.expr))
