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
        return new(statespace, coeff, expr, separate_expectation_values)
    end
    function QAtomProduct(statespace::StateSpace, coeff::T, expr::S, separate_expectation_values::Bool=false) where {T <: CFunction, S <: QAtom}
        return new(statespace, coeff, [expr], separate_expectation_values)
    end
    function QAtomProduct(statespace::StateSpace, expr::AbstractVector{<:QAtom}= QAtom[], separate_expectation_values::Bool=false) 
        return new(statespace, statespace.c_one, expr, separate_expectation_values)
    end
    function QAtomProduct(statespace::StateSpace, expr::S, separate_expectation_values::Bool=false) where {S <: QAtom}
         new(statespace, statespace.c_one, [expr], separate_expectation_values)
    end
end
modify_expr(q::QAtomProduct, expr::Vector{QAtom})::QAtomProduct = QAtomProduct(q.statespace, q.coeff_fun, expr, q.separate_expectation_values)
modify_coeff_expr(q::QAtomProduct, coeff::CFunction, expr::Vector{QAtom})::QAtomProduct = QAtomProduct(q.statespace, coeff, expr, q.separate_expectation_values)
modify_coeff(q::QAtomProduct, coeff::CFunction)::QAtomProduct = QAtomProduct(q.statespace, coeff, q.expr, q.separate_expectation_values)
each_term(q::QAtomProduct) = q.expr
each_coeff(q::QAtomProduct)::Vector{CFunction} = [q.coeff_fun]
get_coeff(q::QComposite) = q.coeff_fun
multiply_coeff(q::QComposite, coeff::CFunction) = modify_coeff(q, get_coeff(q)*coeff)


"""
    QSum

A `QSum` represents the summation of a quantum expression over one or more
index blocks. Each block is a set of summation indexes, with an associated
flag indicating whether the indexes in that block must be all distinct.

It contains:
  - `statespace`: (`StateSpace`) the statespace of the sum.
  - `expr`      : (`QExpr`) the expression being summed over.
  - `blocks`    : (`Vector{Vector{SubSpaceIndex}}`) disjoint, nonempty blocks of summation indexes. Indexes in each block are sorted.
  - `neq_blocks`: (`BitVector`) same length as `blocks`, with `true` meaning all indexes in that block must be different, and `false` meaning they may coincide.

For example:
  ∑_{[i,j]}^{≠} ∑_[k]^{=}} f(i,j,k)

represents a sum where `i` and `j` must be different, while `k` is unrestricted.
"""
struct QSum <: QComposite
    statespace::StateSpace
    expr::QExpr
    eq_indexes::Vector{SubSpaceIndex}            # may be empty
    neq_blocks::Vector{Vector{SubSpaceIndex}}    # each block nonempty and sorted
end
function _QSum(statespace::StateSpace, expr::QExpr, eq_indexes::Vector{SubSpaceIndex}, neq_blocks::Vector{Vector{SubSpaceIndex}})::Vector{QComposite}
    qsum = QSum(statespace, expr, eq_indexes, neq_blocks)
    # run decollision + flatten
    comps = decollision_QSum(qsum)

    if length(comps) == 1 && comps[1] isa QSum
        return QComposite[comps[1]]
    else
        # wrap multiple composites back into an expr
        return comps
    end
end
function _QSum(statespace::StateSpace, expr::QExpr, indexes::Vector{SubSpaceIndex}; neq::Bool=false)::Vector{QComposite}
    isempty(indexes) && return QComposite[expr]
    idxs_sorted = sort(indexes, by=expanded)
    if neq && length(idxs_sorted) > 1
        return _QSum(statespace, expr, SubSpaceIndex[], Vector{Vector{SubSpaceIndex}}([idxs_sorted]))
    else
        return _QSum(statespace, expr, idxs_sorted, Vector{Vector{SubSpaceIndex}}())
    end
end
modify_expr(q::QSum, expr::QExpr) = _QSum(q.statespace, expr, q.eq_indexes, q.neq_blocks)
modify_expr(q::QSum, expr::Vector{QComposite}) = _QSum(q.statespace, QExpr(q.statespace, expr), q.eq_indexes, q.neq_blocks)
modify_expr_indexing(q::QSum, expr::QExpr, eq_indexes::Vector{SubSpaceIndex}, neq_blocks::Vector{Vector{SubSpaceIndex}}) = _QSum(q.statespace, expr, eq_indexes, neq_blocks)
each_term(q::QSum) = q.expr
each_coeff(q::QSum)::Vector{CFunction} = flatmap_to(each_coeff, each_term(q), CFunction)
multiply_coeff(q::QSum, coeff::CFunction)::QSum = modify_expr(q, multiply_coeff(q.expr, coeff) )
get_coeff(q::QSum) = q.statespace.c_one
all_indexes(q::QSum) = vcat(q.eq_indexes, q.neq_blocks...)
iter_all_indexes(q::QSum) = Iterators.flatten((q.eq_indexes, Iterators.flatten(q.neq_blocks)))
iter_all_indexes_with_refs(q::QSum) = Iterators.flatten(( ((q.eq_indexes, i, q.eq_indexes[i]) for i in eachindex(q.eq_indexes)), ((blk, j, blk[j]) for blk in q.neq_blocks for j in eachindex(blk))))   # also returns the current vector and index 
length_all_indexes(q::QSum) = length(q.eq_indexes) + sum(length(block) for block in q.neq_blocks)
function is_single_neq(q::QSum)::Bool 
    (length(q.eq_indexes) == 1 && isempty(q.neq_blocks)) || (isempty(q.eq_indexes) && length(q.neq_blocks) == 1)
end

"""
    Sum(index::Union{String,Symbol,Vector{String},Vector{Symbol}}, expr::QExpr; neq::Bool=false) -> QSum

Constructor of a `QSum` struct. Defines the indexes to sum over, the expressions for which to apply the sum and optionally whether the sum is only over non equal indexes. 
"""
function Sum(statespace::StateSpace, expr::QExpr, eq_indexes::Vector{SubSpaceIndex}, blocks::Vector{Vector{SubSpaceIndex}})::QExpr
    if length(blocks) > 1 
        sort!(blocks)
    end
    all_inds = vcat(eq_indexes, blocks...)
    if length(all_inds) == 0
        return expr 
    end
    new_blocks = Vector{Vector{SubSpaceIndex}}() 
    sizehint!(new_blocks, length(blocks))
    for block in blocks 
        if length(block) <= 1
            append!(eq_indexes, block)
            # remove block from blocks 
        else
            push!(new_blocks, block)
        end
    end
    all = copy(eq_indexes)
    append!(all, reduce(vcat, new_blocks; init=SubSpaceIndex[]))
    if length(sort_unique!(all)) != length(all) 
        error("Sum: indexes must be unique across blocks.")
    end
    return QExpr(statespace, QSum(statespace, expr, eq_indexes, new_blocks))
end
function Sum(statespace::StateSpace, indexes::Union{Vector{String},Vector{Symbol}}, expr::QExpr; neq::Bool=false)::QExpr
    subspace_indexes = SubSpaceIndex.(indexes, Ref(statespace.subspace_info))
    for (i, ind) in enumerate(subspace_indexes)
        for ind2 in subspace_indexes[i+1:end]
            if ind == ind2 
                error("Cannot sum twice over $(indexes[i]).")
            end
        end
    end
    for (sub_ind, index) in zip(subspace_indexes, indexes)
        curr_ensemble = statespace.subspace_info.ensemble_index_by_outer_index[outer(sub_ind)]
        if iszero(curr_ensemble)
            error("Subsystem $index not among ensemble indexes.")
        end
        how_many_non_sum = statespace.subspace_info.how_many_non_sum_by_ensemble[curr_ensemble]
        if how_many_non_sum > inner(sub_ind)
            curr_subspace = statespace.subspaces[outer(sub_ind)]
            possible_keys = curr_subspace.keys[how_many_non_sum+1:end]
            if length(possible_keys) > 0 
                error("Please use a summation index of the ensemble. You used $index, the available summation indexes are $possible_keys.")
            else
                error("No summation indexes defined for this ensemble subspace. 
                        Define in call to SubSpaceDefinitions via Tuple specifying non summation indexes, summation indexes and finally the OperatorSpace. ")
            end
        end
    end
    subspace_indexes = sort(subspace_indexes, by = expanded)
    return QExpr(statespace, _QSum(statespace, expr, subspace_indexes, neq=neq))
end
function Sum(statespace::StateSpace, index::Union{String,Symbol}, expr::QExpr; neq::Bool=false)::QExpr
    return Sum(statespace, [index], expr, neq=neq)
end
""" 
      ∑(index::Union{String,Symbol}, expr::QExpr; neq::Bool=false) -> QSum (use \\sum + Enter)
    sum(index::Union{String,Symbol}, expr::QExpr; neq::Bool=false) -> QSum
Alternative way to call the `Sum` constructor. Sum(index, expr; neq) = ∑(index, expr; neq).
"""
∑(index::Union{String,Symbol}, expr::QExpr; neq::Bool=false) = Sum(expr.statespace, index, expr, neq=neq)
∑(indexes::Union{Vector{String},Vector{Symbol}}, expr::QExpr; neq::Bool=false) = Sum(expr.statespace, indexes, expr, neq=neq)
Base.sum(index::Union{String,Symbol}, expr::QExpr; neq::Bool=false) = Sum(expr.statespace, index, expr, neq=neq)
Base.sum(indexes::Union{Vector{String},Vector{Symbol}}, expr::QExpr; neq::Bool=false) = Sum(expr.statespace, indexes, expr, neq=neq)
 


""" 
    QCompositeProduct

Represents a product of QComposites. 
"""
struct QCompositeProduct <: QMultiComposite
    statespace::StateSpace         # State space of the product.
    coeff_fun::CFunction
    expr::Vector{QComposite} 
end
function _QCompositeProduct(statespace::StateSpace, coeff_fun::CFunction, expr::Vector{QComposite})::Vector{QComposite}
    if length(expr) == 0
        return [IdentityQAtomProduct(statespace, coeff_fun)]
    elseif length(expr) == 1
        return [multiply_coeff(expr[1], coeff_fun)]
    else
        coeff_fun_mod, expr_mod = separate_coeff_qcomposites(expr, statespace) 
        return [QCompositeProduct(statespace, coeff_fun * coeff_fun_mod, expr_mod)]
    end
end
function _QCompositeProduct(statespace::StateSpace, coeff_fun::CFunction, expr::Vector{QComposite}, ::Val{:nosimp})::Vector{QComposite}
    if length(expr) == 0
        return QComposite[IdentityQAtomProduct(statespace, coeff_fun)]
    elseif length(expr) == 1
        return QComposite[multiply_coeff(expr[1], coeff_fun)]
    else
        return QComposite[_QCompositeProduct(statespace, coeff_fun, expr)]
    end
end
modify_expr(q::QCompositeProduct, expr::Vector{QComposite}) = _QCompositeProduct(q.statespace, q.coeff_fun, expr)
modify_expr(q::QCompositeProduct, expr::Vector{QComposite},::Val{:nosimp}) = _QCompositeProduct(q.statespace, q.coeff_fun, expr, Val(:nosimp))
modify_coeff_expr(q::QCompositeProduct, coeff_fun::CFunction, expr::Vector{QComposite}) = _QCompositeProduct(q.statespace, coeff_fun, expr)
modify_coeff_expr(q::QCompositeProduct, coeff_fun::CFunction, expr::Vector{QComposite}, ::Val{:nosimp}) = _QCompositeProduct(q.statespace, coeff_fun, expr, Val(:nosimp))
modify_coeff(q::QCompositeProduct, coeff_fun::CFunction)::QCompositeProduct = QCompositeProduct(q.statespace, coeff_fun, q.expr)

struct QCommutator <: QMultiComposite
    statespace::StateSpace
    coeff_fun::CFunction
    expr::Vector{QExpr}
    function QCommutator(statespace::StateSpace, coeff_fun::CFunction, expr::Vector{QExpr})
        @assert length(expr)==2 "QCommutators require two QExpr, got $(length(expr)) instead."
        new(statespace, coeff_fun, expr)
    end
end
function QCommutator(statespace::StateSpace, q1::QExpr, q2::QExpr, coeff_fun::CFunction)::Vector{QComposite}
    return QComposite[QCommutator(statespace, coeff_fun, QComposite[q1, q2])]
end
function QCommutator(statespace::StateSpace, q1::QExpr, q2::QExpr)::Vector{QComposite}
    return QComposite[QCommutator(statespace, statespace.c_one, QExpr(statespace, QComposite[q1, q2]))]
end
modify_expr(q::QCommutator, expr::Vector{QExpr}) = QComposite[QCommutator(q.statespace, q.coeff_fun, expr)]
modify_coeff_expr(q::QCommutator, coeff_fun::CFunction, expr::Vector{QExpr}) = QComposite[QCommutator(q.statespace, coeff_fun, expr)]
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
function _QExp(statespace::StateSpace, coeff_fun::CFunction, expr::QExpr)::Vector{QComposite}
    if length(expr) == 1 && isa(expr[1], QLog)
        return QComposite[expr[1].expr]
    elseif is_numeric(expr) 
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
        return QComposite[modify_coeff(expr[1], coeff_fun * exp(sum_of_coeff_funs)) ]
    end
    return QComposite[QExp(statespace, coeff_fun, simplify_QExpr(expr))]
end
function exp(q::QExpr)::QExpr
    return QExpr(q.statespace, _QExp(q.statespace, q.statespace.c_one, q))
end
modify_expr(q::QExp, expr::QExpr) = _QExp(q.statespace, q.coeff_fun, expr)
modify_coeff_expr(q::QExp, coeff_fun::CFunction, expr::QExpr) = _QExp(q.statespace, coeff_fun, expr)
modify_coeff(q::QExp, coeff_fun::CFunction) = QExp( q.statespace, coeff_fun, q.expr)
iszero(q::QExp) = iszero(q.coeff_fun) 

struct QLog <: QComposite
    statespace::StateSpace
    coeff_fun::CFunction
    expr::QExpr
end
function _QLog(statespace::StateSpace, coeff_fun::CFunction, expr::QExpr)::Vector{QComposite}
    if length(expr) == 1 && isa(expr[1], QExp)
        return QComposite[expr[1].expr * coeff_fun]
    end
    if is_numeric(expr)
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
        return QComposite[modify_coeff(expr[1], coeff_fun * log(sum_of_coeff_funs)) ]
    end
    return QComposite[QLog(statespace, coeff_fun, simplify_QExpr(expr))]
end
function log(q::QExpr)::QExpr
    return QExpr(q.statespace, _QLog(q.statespace, q.statespace.c_one, q))
end
modify_expr(q::QLog, expr::QExpr) = _QLog(q.statespace, q.coeff_fun, expr)
modify_coeff_expr(q::QLog, coeff_fun::CFunction, expr::QExpr) = _QLog(q.statespace, coeff_fun, expr)
modify_coeff(q::QLog, coeff_fun::CFunction)::QLog = QLog(q.statespace, coeff_fun, q.expr)
iszero(q::QLog) = iszero(q.coeff_fun) #|| isone(q.expr) => that should be autosimplified

struct QPower <: QCompositeN
    statespace::StateSpace
    coeff_fun::CFunction
    n::Int
    expr::QExpr
end
function _QPower(statespace::StateSpace, coeff_fun::CFunction, n::Int, expr::QExpr)::Vector{QComposite}
    if is_numeric(expr)
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
        return QComposite[modify_coeff(expr[1], coeff_fun * power(sum_of_coeff_funs, n))]
    end
    return QComposite[QPower(statespace, coeff_fun, n, expr)]
end
modify_expr(q::QPower, expr::QExpr) = _QPower(q.statespace, q.coeff_fun, q.n, expr)
modify_coeff_expr(q::QPower, coeff_fun::CFunction, expr::QExpr)::QPower = _QPower(q.statespace, coeff_fun, q.n, expr)
modify_coeff(q::QPower, coeff_fun::CFunction) = QPower(q.statespace, coeff_fun, q.n, q.expr)

""" 
    power(q::QExpr, n::Int)::QExpr

Returns the q^n of a quantum expression, without expanding the expression.
"""
function power(q::QExpr, n::Int)::QExpr
    if n == 1 
         return q 
    end
    return QExpr(q.statespace, _QPower(q.statespace.c_one, n, q))
end


struct QRoot <: QCompositeN
    statespace::StateSpace
    coeff_fun::CFunction
    n::Int
    expr::QExpr
end
function _QRoot(statespace::StateSpace, coeff_fun::CFunction, n::Int, expr::QExpr)::Vector{QComposite}
    if is_numeric(expr)
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
        return QComposite[modify_coeff(expr[1], coeff_fun * root(sum_of_coeff_funs, n))]
    end 
    return QComposite[QRoot(statespace, coeff_fun, n, expr)]
end
""" 
    root(q::QExpr, n::Int)::QExpr

Returns the n'th root of a quantum expression q^{1/n}, without expanding the expression.
"""
function root(q::QExpr, n::Int=2)::QExpr
    if n == 1 
        return q 
    end
    return QExpr(q.statespace, _QRoot(q.statespace, q.statespace.c_one, n, q))
end
function sqrt(q::QExpr)::QExpr
    return QExpr(q.statespace, _QRoot(q.statespace, q.statespace.c_one, 2, q))
end
modify_expr(q::QRoot, expr::QExpr) = _QRoot(q.statespace, q.coeff_fun, q.n, expr)
modify_coeff_expr(q::QRoot, coeff_fun::CFunction, expr::QExpr)::QRoot = QRoot(coeff_fun, q.n, expr)
modify_coeff(q::QRoot, coeff_fun::CFunction) = QRoot(q.statespace, coeff_fun, q.n, q.expr)
