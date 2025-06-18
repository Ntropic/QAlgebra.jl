import Base: exp, log

export qCommutator, qExp, qLog, qPower, power, qRoot

mutable struct qCommutator <: qMultiComposite
    statespace::StateSpace
    expr::Vector{qExpr}
end
function qCommutator(q1::qExpr, q2::qExpr)::qExpr
    return qExpr(qCommutator(q1.statespace, qExpr[q1, q2]))
end
copy(q::qCommutator) = qCommutator(q.statespace, [copy(x) for x in q.expr])


### Non-simple qFunctions 
mutable struct qExp <: qComposite
    statespace::StateSpace
    expr::qExpr
end
function exp(q::qExpr)::qExpr
    return qExpr(qExp(q.statespace, q))
end
copy(q::qExp) = qExp(q.statespace, copy(q.expr))


mutable struct qLog <: qComposite
    statespace::StateSpace
    expr::qExpr
end
function log(q::qExpr)::qExpr
    return qExpr(qLog(q.statespace, q))
end
copy(q::qLog) = qExp(q.statespace, copy(q.expr))


mutable struct qPower <: qComposite
    statespace::StateSpace
    n::Int
    expr::qExpr
end
copy(q::qPower) = qExp(q.statespace, n, copy(q.expr))

function power(q::qExpr, n::Int)::qExpr
    if n == 1 
         return q 
    end
    return qExpr(qPower(q.statespace, n, q))
end


mutable struct qRoot <: qComposite
    statespace::StateSpace
    n::Int
    expr::qExpr
end
function root(q::qExpr, n::Int=2)::qExpr
    if n == 1 
        return q 
    end
    return qExpr(qRoot(q.statespace, n, q))
end
copy(q::qRoot) = qExp(q.statespace, n, copy(q.expr))