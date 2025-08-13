import Base: exp, log, sqrt

export QCommutator, QExp, QLog, QPower, power, QRoot, root

struct QCommutator <: QMultiComposite
    statespace::StateSpace
    coeff_fun::CFunction
    expr::Vector{QExpr}
    function QCommutator(q1::QExpr, q2::QExpr)
        new(QCommutator(q1.statespace, copy(q1.statespace.fone), QExpr[q1, q2]))
    end
    function QCommutator(q1::QExpr, q2::QExpr, ::Val{:simp})
        new(QCommutator(q1.statespace, copy(q1.statespace.fone), QExpr[simplify_QExpr(q1), simplify_QExpr(q2)]))
    end
end
copy(q::QCommutator) = QCommutator(q.statespace, q.coeff_fun, q.expr)
modify_expr(q::QCommutator, expr::Vector{QExpr}) = QCommutator(q.statespace, q.coeff_fun, expr)
modify_coeff_expr(q::QCommutator, coeff_fun::CFunction, expr::Vector{QExpr}) = QCommutator(q.statespace, coeff_fun, expr)
modify_coeff(q::QCommutator, coeff_fun::CFunction) = QCommutator(q.statespace, coeff_fun, q.expr)

""" 
    QCompositeProduct

Represents a product of QComposites. 
"""
struct QCompositeProduct <: QMultiComposite
    statespace::StateSpace         # State space of the product.
    coeff_fun::CFunction
    expr::Vector{QExpr} 
    function QCompositeProduct(statespace::StateSpace, coeff_fun::CFunction, expr::Vector{QExpr} )
        new(statespace, copy(coeff_fun), expr)
    end
    function QCompositeProduct(coeff_fun::CFunction, expr::Vector{QExpr} )
        new(expr[1].statespace, copy(coeff_fun), copy.(expr))
    end
    function QCompositeProduct(coeff_fun::CFunction, expr::Vector{QComposite} )
        new(expr[1].statespace, copy(coeff_fun), [QExpr([x]) for x in expr])
    end
    function QCompositeProduct(statespace::StateSpace, expr::Vector{QExpr} )
        new(statespace, copy(statespace.fone), copy.(expr))
    end
    function QCompositeProduct(expr::Vector{QExpr} )
        new(expr[1].statespace, copy(expr[1].statespace.fone), copy.(expr))
    end
    function QCompositeProduct(expr::Vector{QComposite} )
        new(expr[1].statespace, copy(expr[1].statespace.fone), [QExpr([x]) for x in expr])
    end
    function QCompositeProduct(statespace::StateSpace, coeff_fun::CFunction, expr::Vector{QExpr}, ::Val{:simp})
        new(statespace, copy(coeff_fun), simplify_QExpr(expr))
    end
    function QCompositeProduct(coeff_fun::CFunction, expr::Vector{QExpr} , ::Val{:simp})
        new(expr[1].statespace, copy(coeff_fun), simplify_QExpr(expr))
    end
    function QCompositeProduct(statespace::StateSpace, expr::Vector{QExpr} , ::Val{:simp})
        new(statespace, copy(statespace.fone), simplify_QExpr(expr))
    end
    function QCompositeProduct(expr::Vector{QExpr} , ::Val{:simp})
        new(expr[1].statespace, copy(expr[1].statespace.fone), simplify_QExpr(expr))
    end
end
copy(q::QCompositeProduct) = QCompositeProduct(q.statespace, q.coeff_fun, q.expr)
modify_expr(q::QCompositeProduct, expr::Vector{QExpr}) = QCompositeProduct(q.statespace, q.coeff_fun, expr)
modify_coeff_expr(q::QCompositeProduct, coeff_fun::CFunction, expr::Vector{QExpr}) = QCompositeProduct(q.statespace, coeff_fun, expr)
modify_coeff(q::QCompositeProduct, coeff_fun::CFunction) = QCompositeProduct(q.statespace, coeff_fun, q.expr)

### Non-simple qFunctions 
struct QExp <: QComposite
    statespace::StateSpace
    coeff_fun::CFunction
    expr::QExpr
    function QExp(statespace::StateSpace, coeff_fun::CFunction, expr::QExpr)
        if length(expr) == 1 && isa(expr[1], QLog)
            return expr[1].expr
        end
        if is_numeric(expr) 
            sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
            return modify_coeff(expr[1], coeff_fun * exp(sum_of_coeff_funs)) 
        end
        new(statespace, copy(coeff_fun), simplify_QExpr(expr))
    end
end
function exp(q::QExpr)::QExpr
    return QExpr([QExp(q.statespace, q.statespace.fone, q)])
end
copy(q::QExp) = QExp(q.statespace, q.coeff_fun, q.expr)
modify_expr(q::QExp, expr::QExpr) = QExp(q.statespace, q.coeff_fun, expr)
modify_coeff_expr(q::QExp, coeff_fun::CFunction, expr::QExpr) = QExp(q.statespace, coeff_fun, expr)
modify_coeff(q::QExp, coeff_fun::CFunction) = QExp(q.statespace, coeff_fun, q.expr)
iszero(q::QExp) = iszero(q.coeff_fun) 

struct QLog <: QComposite
    statespace::StateSpace
    coeff_fun::CFunction
    expr::QExpr
    function QLog(statespace::StateSpace, coeff_fun::CFunction, expr::QExpr)
        if length(expr) == 1 && isa(expr[1], QExp)
            return expr[1].expr
        end
        if is_numeric(expr)
            sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
            return modify_coeff(expr[1], coeff_fun * log(sum_of_coeff_funs)) 
        end
        return new(statespace, copy(coeff_fun), simplify_QExpr(expr))
    end
end
function log(q::QExpr)::QExpr
    return QExpr([QLog(q.statespace, q.statespace.fone, q)])
end
copy(q::QLog) = QLog(q.statespace, q.coeff_fun, q.expr)
modify_expr(q::QLog, expr::QExpr) = QLog(q.statespace, q.coeff_fun, expr)
modify_coeff_expr(q::QLog, coeff_fun::CFunction, expr::QExpr) = QLog(q.statespace, coeff_fun, expr)
modify_coeff(q::QLog, coeff_fun::CFunction) = QLog(q.statespace, coeff_fun, q.expr)
iszero(q::QLog) = iszero(q.coeff_fun) #|| isone(q.expr) => that should be autosimplified

struct QPower <: QCompositeN
    statespace::StateSpace
    coeff_fun::CFunction
    n::Int
    expr::QExpr
    function QPower(statespace::StateSpace, coeff_fun::CFunction, n::Int, expr::QExpr)
        if is_numeric(expr)
            sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
            return modify_coeff(expr[1], coeff_fun * power(sum_of_coeff_funs, n))
        end
        return new(statespace, copy(coeff_fun), n, simplify_QExpr(expr))
    end
end
copy(q::QPower) = QPower(q.statespace, q.coeff_fun, q.n, q.expr)
modify_expr(q::QPower, expr::QExpr) = QPower(q.statespace, q.coeff_fun, q.n, expr)
modify_coeff_expr(q::QPower, coeff_fun::CFunction, expr::QExpr) = QPower(q.statespace, coeff_fun, q.n, expr)
modify_coeff(q::QPower, coeff_fun::CFunction) = QPower(q.statespace, coeff_fun, q.n, q.expr)

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


struct QRoot <: QCompositeN
    statespace::StateSpace
    coeff_fun::CFunction
    n::Int
    expr::QExpr
    function QRoot(statespace::StateSpace, coeff_fun::CFunction, n::Int, expr::QExpr)
        if is_numeric(expr)
            sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
            return modify_coeff(expr[1], coeff_fun * root(sum_of_coeff_funs, n))
        end 
        new(statespace, copy(coeff_fun), n, simplify_QExpr(expr))
    end
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
copy(q::QRoot) = QRoot(q.statespace, q.coeff_fun, q.n, q.expr)
modify_expr(q::QRoot, expr::QExpr) = QRoot(q.statespace, q.coeff_fun, q.n, expr)
modify_coeff_expr(q::QRoot, coeff_fun::CFunction, expr::QExpr) = QRoot(q.statespace, coeff_fun, q.n, expr)
modify_coeff(q::QRoot, coeff_fun::CFunction) = QRoot(q.statespace, coeff_fun, q.n, q.expr)
