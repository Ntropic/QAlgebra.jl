import Base: exp, log, sqrt

export QAtomProduct, QSum, Sum, ∑, QCommutator, QCompositeProduct, QExp, QLog, QPower, power, QRoot, root


""" 
    QAtomProduct

A product of QAtom expressions, i.e. qTerms or QAbstract.
It contains:
    - `statespace`: The state space in which the product is defined.
    - `coeff_fun`: The function of parameters for the Operator product
    - `expr`: A vector of qAtoms (qTerms or QAbstract) that are multiplied together.
"""
struct QAtomProduct <: QComposite
    statespace::StateSpace         # State space of the product.
    coeff_fun::CFunction            # function of scalar parameters => has +,-,*,/,^ defined 
    expr::Vector{QAtom}             # Vector of qAtoms (qTerms or QAbstract).
    separate_expectation_values::Bool 
    function QAtomProduct(statespace::StateSpace, coeff::T, expr::AbstractVector{<:QAtom}= QAtom[], separate_expectation_values::Bool=false) where T <: CFunction
        new(statespace, coeff, expr, separate_expectation_values)
    end
    function QAtomProduct(statespace::StateSpace, coeff::T, expr::S, separate_expectation_values::Bool=false) where {T <: CFunction, S <: QAtom}
        new(statespace, coeff, [expr], separate_expectation_values)
    end
    function QAtomProduct(statespace::StateSpace, expr::AbstractVector{<:QAtom}= QAtom[], separate_expectation_values::Bool=false) 
        new(statespace, statespace.c_one, expr, separate_expectation_values)
    end
    function QAtomProduct(statespace::StateSpace, expr::S, separate_expectation_values::Bool=false) where {S <: QAtom}
        new(statespace, statespace.c_one, [expr], separate_expectation_values)
    end
end
modify_expr(q::QAtomProduct, expr::Vector{QAtom})::QAtomProduct = QAtomProduct(q.statespace, q.coeff, expr, q.separate_expectation_values)
modify_coeff_expr(q::QAtomProduct, coeff::CFunction, expr::Vector{QAtom})::QAtomProduct = QAtomProduct(q.statespace, coeff, expr, q.separate_expectation_values)
modify_coeff(q::QAtomProduct, coeff::CFunction)::QAtomProduct = QAtomProduct(q.statespace, coeff, q.expr, q.separate_expectation_values)
each_term(q::QAtomProduct) = q.expr
each_coeff(q::QAtomProduct)::Vector{CFunction} = [q.coeff]
get_coeff(q::QComposite) = q.coeff_fun


"""
    QSum

A `QSum` represents the summation of a quantum Equation over indexes in a quantum expression.
It contains:
    - `expr`: The expression being summed over, which is a `QExpr` object.
    - `indexes`: A vector of strings representing the summation indexes (e.g., "i").
    - `subsystem_index`: The index of the subspace in which the indexes live. 
    - `element_indexes`: A vector of integers representing the position of the indexes in that subspace.
    - `neq`: A boolean indicating whether different indexes in the sum can refer to the same element in the subspace. 
            For example, the indexes i,j,k can refer to different elements in a much larger bath of elements. 
"""
struct QSum <: QComposite
    statespace::StateSpace
    expr::QExpr       # The expression being summed over.    # use expr in other QComposites except for QAtomProduct
    indexes::Vector{SubSpaceIndex}
    neq::Bool
end
function QSum(expr::QExpr, indexes::Vector{SubSpaceIndex}, neq::Base.Bool)
    statespace = expr.statespace
    if length(indexes) == 0
        return expr
    end
    return new(statespace,  expr, copy(indexes), subsystem_index, copy(element_indexes), neq)
end
copy(q::QSum)::QSum = QSum(q.expr, q.indexes, q.neq)
modify_expr(q::QSum, expr::QExpr) = QSum(expr, q.indexes, q.neq)
modify_expr_indexes(q::QSum, expr::QExpr, indexes::Vector{SubSpaceIndex}) = QSum(expr, indexes, q.neq)
each_term(q::QSum) = q.expr
each_coeff(q::QSum)::Vector{CFunction} = flatmap_to(each_coeff, each_term(q), CFunction)
modify_coeff(q::QSum, coeff::CFunction)::QSum = QSum(q.expr*coeff, q.indexes, q.neq)
get_coeff(q::QSum) = q.statespace.c_one

"""
    Sum(index::Union{String,Symbol,Vector{String},Vector{Symbol}}, expr::QExpr; neq::Bool=false) -> QSum

Constructor of a `QSum` struct. Defines the indexes to sum over, the expressions for which to apply the sum and optionally whether the sum is only over non equal indexes. 
"""
function Sum(indexes::Union{Vector{String},Vector{Symbol}}, expr::QExpr; neq::Bool=false)::QExpr
    statespace = expr.statespace
    subspace_indexes = SubSpaceIndex.(indexes, Ref(statespace.subspace_info))
    return QExpr(statespace, [QSum(statespace, expr, subspace_indexes, neq)])
end
function Sum(index::Union{String,Symbol}, expr::QExpr; neq::Bool=false)::QExpr
    return Sum([index], expr, neq=neq)
end
""" 
    ∑(index::Union{String,Symbol}, expr::QExpr; neq::Bool=false) -> QSum

Alternative way to call the `Sum` constructor. Sum(index, expr; neq) = ∑(index, expr; neq).
"""
∑(index::Union{String,Symbol}, expr::QExpr; neq::Bool=false) = Sum(index, expr, neq=neq)
∑(indexes::Union{Vector{String},Vector{Symbol}}, expr::QExpr; neq::Bool=false) = Sum(indexes, expr, neq=neq)
 


""" 
    QCompositeProduct

Represents a product of QComposites. 
"""
struct QCompositeProduct <: QMultiComposite
    statespace::StateSpace         # State space of the product.
    coeff_fun::CFunction
    expr::Vector{QComposite} 
end
function QCompositeProductCleanup(ss::StateSpace, coeff_fun::CFunction, expr::Vector{QComposite})
    if length(expr) == 0
        return IdentityQAtomProduct(ss, coeff_fun)
    elseif length(expr) == 1
        return modify_coeff(expr[1], coeff_fun*get_coeff(expr[1]))
    else
        coeff_fun_mod, expr_mod = separate_coeff_qcomposites(expr, statespace) 
        return QCompositeProduct(expr[1].statespace, coeff_fun * coeff_fun_mod, expr_mod)
    end
end
function QCompositeProductCleanup(ss::StateSpace, coeff_fun::CFunction, expr::Vector{QComposite}, ::Val{:nosimp})
    if length(expr) == 0
        return IdentityQAtomProduct(ss, coeff_fun)
    elseif length(expr) == 1
        return modify_coeff(expr[1], coeff_fun*get_coeff(expr[1]))
    else
        return QCompositeProduct(ss, coeff_fun, expr)
    end
end
modify_expr(q::QCompositeProduct, expr::Vector{QComposite}) = QCompositeProductCleanup(q.statespace, q.coeff_fun, expr)
modify_expr(q::QCompositeProduct, expr::Vector{QComposite},::Val{:nosimp}) = QCompositeProductCleanup(q.statespace, q.coeff_fun, expr, Val(:nosimp))
modify_coeff_expr(q::QCompositeProduct, coeff_fun::CFunction, expr::Vector{QComposite}) = QCompositeProductCleanup(q.statespace, coeff_fun, expr)
modify_coeff_expr(q::QCompositeProduct, coeff_fun::CFunction, expr::Vector{QComposite}, ::Val{:nosimp}) = QCompositeProductCleanup(q.statespace, coeff_fun, expr, Val(:nosimp))
modify_coeff(q::QCompositeProduct, coeff_fun::CFunction)::QCompositeProduct = QCompositeProductCleanup(q.statespace, coeff_fun, q.expr, Val(:nosimp))

struct QCommutator <: QMultiComposite
    statespace::StateSpace
    coeff_fun::CFunction
    expr::Vector{QExpr}
    function QCommutator(q1::QExpr, q2::QExpr)
        new(QCommutator(q1.statespace, copy(q1.statespace.c_one), QExpr[q1, q2]))
    end
    function QCommutator(q1::QExpr, q2::QExpr, ::Val{:simp})
        new(QCommutator(q1.statespace, copy(q1.statespace.c_one), QExpr[simplify_QExpr(q1), simplify_QExpr(q2)]))
    end
end
modify_expr(q::QCommutator, expr::Vector{QExpr}) = QCommutator(q.statespace, q.coeff_fun, expr)
modify_coeff_expr(q::QCommutator, coeff_fun::CFunction, expr::Vector{QExpr}) = QCommutator(q.statespace, coeff_fun, expr)
modify_coeff(q::QCommutator, coeff_fun::CFunction)::QCommutator = QCommutator(q.statespace, coeff_fun, q.expr)
each_term(q::QMultiComposite) = q.expr
each_coeff(q::QMultiComposite)::Vector{CFunction} = CFunction[q.coeff_fun; flatmap_to(each_coeff, each_term(q), CFunction)]
each_term(q::QComposite) = [q.expr]
each_coeff(q::QComposite)::Vector{CFunction} = CFunction[q.coeff_fun; each_coeff(q.expr)]

### Non-simple qFunctions 
struct QExp <: QComposite
    statespace::StateSpace
    coeff_fun::CFunction
    expr::QExpr
end
function QExp(coeff_fun::CFunction, expr::QExpr)
    statespace = expr.statespace
    if length(expr) == 1 && isa(expr[1], QLog)
        return expr[1].expr
    elseif is_numeric(expr) 
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
        return modify_coeff(expr[1], coeff_fun * exp(sum_of_coeff_funs)) 
    end
    return QExp(statespace, coeff_fun, simplify_QExpr(expr))
end
function exp(q::QExpr)::QExpr
    return QExpr([QExp(q.statespace.c_one, q)])
end
modify_expr(q::QExp, expr::QExpr) = QExp(q.coeff_fun, expr)
modify_coeff_expr(q::QExp, coeff_fun::CFunction, expr::QExpr) = QExp(coeff_fun, expr)
modify_coeff(q::QExp, coeff_fun::CFunction) = QExp( coeff_fun, q.expr)
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
        return new(statespace, coeff_fun, simplify_QExpr(expr))
    end
end
function log(q::QExpr)::QExpr
    return QExpr([QLog(q.statespace, q.statespace.c_one, q)])
end
modify_expr(q::QLog, expr::QExpr) = QLog(q.statespace, q.coeff_fun, expr)
modify_coeff_expr(q::QLog, coeff_fun::CFunction, expr::QExpr) = QLog(q.statespace, coeff_fun, expr)
modify_coeff(q::QLog, coeff_fun::CFunction)::QLog = QLog(q.statespace, coeff_fun, q.expr)
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
        return new(statespace, coeff_fun, n, simplify_QExpr(expr))
    end
end
modify_expr(q::QPower, expr::QExpr) = QPower(q.statespace, q.coeff_fun, q.n, expr)
modify_coeff_expr(q::QPower, coeff_fun::CFunction, expr::QExpr)::QPower = QPower(q.statespace, coeff_fun, q.n, expr)
modify_coeff(q::QPower, coeff_fun::CFunction) = QPower(q.statespace, coeff_fun, q.n, q.expr)

""" 
    power(q::QExpr, n::Int)::QExpr

Returns the q^n of a quantum expression, without expanding the expression.
"""
function power(q::QExpr, n::Int)::QExpr
    if n == 1 
         return q 
    end
    return QExpr([QPower(q.statespace, q.statespace.c_one, n, q)])
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
        new(statespace, coeff_fun, n, simplify_QExpr(expr))
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
    return QExpr([QRoot(q.statespace, q.statespace.c_one, n, q)])
end
function sqrt(q::QExpr)::QExpr
    return QExpr([QRoot(q.statespace, q.statespace.c_one, 2, q)])
end
modify_expr(q::QRoot, expr::QExpr) = QRoot(q.statespace, q.coeff_fun, q.n, expr)
modify_coeff_expr(q::QRoot, coeff_fun::CFunction, expr::QExpr)::QRoot = QRoot(q.statespace, coeff_fun, q.n, expr)
modify_coeff(q::QRoot, coeff_fun::CFunction) = QRoot(q.statespace, coeff_fun, q.n, q.expr)
