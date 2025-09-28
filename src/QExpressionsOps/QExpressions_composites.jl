import Base: exp, log, sqrt
import ..QAlgebra: sort_unique!


export QAtomProduct, QAtomOrdered, OrderedQAtomProduct, permutation, QSum, ∑, NeqConstraint, neq, QCommutator, QCompositeProduct, QExp, QLog, QPower, power, QRoot, root

const SumIndexInput = Union{Symbol, String, SubSpaceIndex}

"""
    QAtomProduct

A product of QAtom expressions, i.e. qTerms or QAbstract.
It contains:
    - `qspace`      : The quantum space in which the product is defined.
    - `coeff_fun`   : The function of parameters for the Operator product
    - `expr`        : A vector of qAtoms (qTerms or QAbstract) that are multiplied together.
    - `separate_expectation_values` : Flag indicating whether expectation values factorise within cumulants.
    - `braket`       : Indicates whether the product should be rendered inside expectation brackets ⟨⋯⟩.
"""
struct QAtomProduct <: QComposite
    qspace::QSpace         # State space of the product.
    coeff_fun::CFunction            # function of scalar parameters => has +,-,*,/,^ defined 
    expr::Vector{QAtom}             # Vector of qAtoms (qTerms or QAbstract).
    separate_expectation_values::Bool 
    braket::Bool

    function QAtomProduct(qspace::QSpace, coeff::T, expr::AbstractVector{<:QAtom}=QAtom[], separate_expectation_values::Bool=false, braket::Bool=false) where T <: CFunction
        return new(qspace, coeff, Vector{QAtom}(expr), separate_expectation_values, braket)
    end
    function QAtomProduct(qspace::QSpace, coeff::T, expr::S, separate_expectation_values::Bool=false) where {T <: CFunction, S <: QAtom}
        return QAtomProduct(qspace, coeff, QAtom[expr], separate_expectation_values)
    end
    function QAtomProduct(qspace::QSpace, expr::AbstractVector{<:QAtom}= QAtom[], separate_expectation_values::Bool=false) 
        return QAtomProduct(qspace, qspace.c_one, expr, separate_expectation_values)
    end
    function QAtomProduct(qspace::QSpace, expr::S, separate_expectation_values::Bool=false) where {S <: QAtom}
        return QAtomProduct(qspace, qspace.c_one, QAtom[expr], separate_expectation_values)
    end
end

modify_expr(q::QAtomProduct, expr::Vector{QAtom})::Vector{QComposite} =
    QComposite[QAtomProduct(q.qspace, q.coeff_fun, expr, q.separate_expectation_values, q.braket)]
modify_expr(q::QAtomProduct, expr::Vector{QAtom}, ::Val{:nosimp})::Vector{QComposite} =
    QComposite[QAtomProduct(q.qspace, q.coeff_fun, expr, q.separate_expectation_values, q.braket)]
modify_coeff_expr(q::QAtomProduct, coeff::CFunction, expr::Vector{QAtom})::QAtomProduct =
    QAtomProduct(q.qspace, coeff, expr, q.separate_expectation_values, q.braket)
modify_coeff(q::QAtomProduct, coeff::CFunction)::QAtomProduct =
    QAtomProduct(q.qspace, coeff, q.expr, q.separate_expectation_values, q.braket)
each_term(q::QAtomProduct) = q.expr
each_coeff(q::QAtomProduct)::Vector{CFunction} = [q.coeff_fun]
get_coeff(q::QComposite) = q.coeff_fun
multiply_coeff(q::QComposite, coeff::CFunction) = modify_coeff(q, get_coeff(q)*coeff)
set_braket(q::QAtomProduct, val::Bool=true) = QAtomProduct(q.qspace, q.coeff_fun, q.expr, q.separate_expectation_values, val)


include("ConstrainedIndexes.jl")
"""
    QSum

Concrete representation of a quantum sum. Each instance stores one
`ConstrainedIndexBlock` per ensemble of the ambient `QSpace`, so all bound summation
indexes and their constraint matrices live inside `blocks`.
It contains:
  - `qspace::QSpace`: ambient space that supplies ensemble metadata.
  - `expr::QExpr`: body of the summation.
  - `blocks::Vector{ConstrainedIndexBlock}`: per-ensemble containers of indexes,
    cached ensemble slot numbers, and their constraint bit vectors. Empty
    blocks denote ensembles with no bound indexes.
"""
struct QSum <: QComposite
    qspace::QSpace
    expr::QExpr
    blocks::Vector{ConstrainedIndexBlock}
    function QSum(qspace::QSpace, expr::QExpr, blocks::Vector{ConstrainedIndexBlock})
        #info = qspace.subspace_info
        #length(blocks) == length(info.where_ensembles) || error("One block per ensemble required.")
        return new(qspace, expr, blocks)
    end
end

function _QSum(qspace::QSpace, expr::QExpr, blocks::Vector{ConstrainedIndexBlock})::Vector{QComposite}
    if isempty(blocks) || all(isempty(block.indexes) for block in blocks)
        return QComposite[expr]
    end
    qsum = QSum(qspace, expr, blocks)
    # TODO: reinstate decollision_QSum once it supports block-based construction
    # return decollision_QSum(qsum)
    return QComposite[qsum]
end
modify_expr(q::QSum, expr::QExpr, ::Val{:nodecollision}) = QComposite[QSum(q.qspace, expr, q.blocks)]
modify_expr(q::QSum, expr::Vector{QComposite}, ::Val{:nodecollision}) = QComposite[QSum(q.qspace, QExpr(q.qspace, expr), q.blocks)]
modify_expr(q::QSum, expr::QExpr, ::Val{:nosimp}) = QComposite[QSum(q.qspace, expr, q.blocks)]
modify_expr(q::QSum, expr::Vector{QComposite}, ::Val{:nosimp}) = QComposite[QSum(q.qspace, QExpr(q.qspace, expr, Val(:nosimp)), q.blocks)]
modify_expr(q::QSum, expr::QExpr) = _QSum(q.qspace, expr, q.blocks)
modify_expr(q::QSum, expr::Vector{QComposite}) = _QSum(q.qspace, QExpr(q.qspace, expr), q.blocks)
each_term(q::QSum) = q.expr
each_coeff(q::QSum)::Vector{CFunction} = flatmap_to(each_coeff, each_term(q), CFunction)
multiply_coeff(q::QSum, coeff::CFunction)::QSum = only(modify_expr(q, multiply_coeff(q.expr, coeff)))
get_coeff(q::QSum) = q.qspace.c_one

all_indexes(q::QSum)::Vector{SubSpaceIndex} = _flatten_indexes(q.blocks)
iter_all_indexes(q::QSum) = Base.Iterators.flatten((block.indexes for block in q.blocks))
iter_all_constraints(q::QSum) = Base.Iterators.flatten((block.constraints for block in q.blocks))
container_iter_all_indexes_with_refs(q::QSum) = Base.Iterators.flatten( ((block.indexes, i, block.indexes[i]) for i in eachindex(block.indexes)) for block in q.blocks)
length_all_indexes(q::QSum) = sum(length(block) for block in q.blocks)

_to_subspace_index(qspace::QSpace, idx::SubSpaceIndex) = idx
_to_subspace_index(qspace::QSpace, idx::Union{Symbol,String}) = SubSpaceIndex(idx, qspace.subspace_info)

indexes2subspaceindexes(qspace::QSpace, idx::SumIndexInput)::Vector{SubSpaceIndex} = SubSpaceIndex[_to_subspace_index(qspace, idx)]
indexes2subspaceindexes(qspace::QSpace, idxs::AbstractVector{<:SumIndexInput})::Vector{SubSpaceIndex} = SubSpaceIndex[_to_subspace_index(qspace, idx) for idx in idxs]

function _make_qsum_expr(qspace::QSpace, raw_indexes::Union{SumIndexInput, AbstractVector{<:SumIndexInput}}, expr::QExpr, constraints::Vararg{NeqConstraint}; neq::Bool=false)
    indexes = indexes2subspaceindexes(qspace, raw_indexes)
    isempty(indexes) && return expr
    normed_constraints::Vector{NeqConstraint{SubSpaceIndex}} = neqconstraint_of_SubSpaceIndex.(Ref(qspace), collect(constraints))
    blocks = sort_indexes_and_constraints_into_ensemble_blocks(qspace, indexes, neq, normed_constraints)
    comps = _QSum(qspace, expr, blocks)
    return QExpr(qspace, comps)
end

"""
    ∑(indexes, expr, constraints...; neq=false)
    Base.sum(indexes, expr, constraints...; neq=false)

Construct a `QSum` by binding summation indexes in `expr`. Indexes may be
supplied as a single entry or vector, where each entry is a `String`, `Symbol`,
or `SubSpaceIndex`. Additional `NeqConstraint`s restrict coinciding ensemble
slots, while `neq=true` forces all tracked summation indexes within each
ensemble to be distinct. Constraints are supplied positionally as extra
arguments (e.g. `∑(:i, expr, neq(:i, :j))`). All provided summation indexes
must be unique members of the given `qspace`.
"""
∑(index::SumIndexInput, expr::QExpr, constraints::NeqConstraint...; neq::Bool=false) = _make_qsum_expr(expr.qspace, index, expr, constraints...; neq=neq)
∑(indexes::AbstractVector{<:SumIndexInput}, expr::QExpr, constraints::NeqConstraint...; neq::Bool=false) = _make_qsum_expr(expr.qspace, indexes, expr, constraints...; neq=neq)

Base.sum(index::SumIndexInput, expr::QExpr, constraints::NeqConstraint...; neq::Bool=false) = _make_qsum_expr(expr.qspace, index, expr, constraints...; neq=neq)
Base.sum(indexes::AbstractVector{<:SumIndexInput}, expr::QExpr, constraints::NeqConstraint...; neq::Bool=false) = _make_qsum_expr(expr.qspace, indexes, expr, constraints...; neq=neq)



"""
    QCompositeProduct(qspace, coeff_fun, exprs)

Product of quantum composites collected as a single node. Simplifies scalar
coefficients while leaving the nested structure intact.
"""
struct QCompositeProduct <: QMultiComposite
    qspace::QSpace         # State space of the product.
    coeff_fun::CFunction
    expr::Vector{QComposite} 
end
function _QCompositeProduct(qspace::QSpace, coeff_fun::CFunction, expr::Vector{QComposite})::Vector{QComposite}
    if length(expr) == 0
        return [IdentityQAtomProduct(qspace, coeff_fun)]
    elseif length(expr) == 1
        return [multiply_coeff(expr[1], coeff_fun)]
    else
        coeff_fun_mod, expr_mod = separate_coeff_qcomposites(expr, qspace) 
        return [QCompositeProduct(qspace, coeff_fun * coeff_fun_mod, expr_mod)]
    end
end
function _QCompositeProduct(qspace::QSpace, coeff_fun::CFunction, expr::Vector{QComposite}, ::Val{:nosimp})::Vector{QComposite}
    if length(expr) == 0
        return QComposite[IdentityQAtomProduct(qspace, coeff_fun)]
    elseif length(expr) == 1
        return QComposite[multiply_coeff(expr[1], coeff_fun)]
    else
        return _QCompositeProduct(qspace, coeff_fun, expr)  # already returns a vector
    end
end
modify_expr(q::QCompositeProduct, expr::Vector{QComposite}) = _QCompositeProduct(q.qspace, q.coeff_fun, expr)
modify_expr(q::QCompositeProduct, expr::Vector{QComposite},::Val{:nosimp}) = _QCompositeProduct(q.qspace, q.coeff_fun, expr, Val(:nosimp))
modify_coeff_expr(q::QCompositeProduct, coeff_fun::CFunction, expr::Vector{QComposite}) = _QCompositeProduct(q.qspace, coeff_fun, expr)
modify_coeff_expr(q::QCompositeProduct, coeff_fun::CFunction, expr::Vector{QComposite}, ::Val{:nosimp}) = _QCompositeProduct(q.qspace, coeff_fun, expr, Val(:nosimp))
modify_coeff(q::QCompositeProduct, coeff_fun::CFunction)::QCompositeProduct = QCompositeProduct(q.qspace, coeff_fun, q.expr)

"""
    QCommutator(qspace, exprs; coeff_fun=qspace.c_one)

Composite holding two `QExpr` factors representing a commutator. The actual
commutator algebra is implemented in [`Commutator`](@ref); this struct simply
stores the symbolic tuple `[A, B]` with an optional scalar prefactor.
"""
struct QCommutator <: QMultiComposite
    qspace::QSpace
    coeff_fun::CFunction
    expr::Vector{QExpr}
    function QCommutator(qspace::QSpace, coeff_fun::CFunction, expr::Vector{QExpr})
        @assert length(expr)==2 "QCommutators require two QExpr, got $(length(expr)) instead."
        new(qspace, coeff_fun, expr)
    end
end
function QCommutator(qspace::QSpace, q1::QExpr, q2::QExpr, coeff_fun::CFunction)::Vector{QComposite}
    return QComposite[QCommutator(qspace, coeff_fun, QComposite[q1, q2])]
end
function QCommutator(qspace::QSpace, q1::QExpr, q2::QExpr)::Vector{QComposite}
    return QComposite[QCommutator(qspace, qspace.c_one, QExpr(qspace, QComposite[q1, q2]))]
end
modify_expr(q::QCommutator, expr::Vector{QExpr}) = QComposite[QCommutator(q.qspace, q.coeff_fun, expr)]
modify_expr(q::QCommutator, expr::Vector{QExpr}, ::Val{:nosimp}) =
    QComposite[QCommutator(q.qspace, q.coeff_fun, expr)]
modify_coeff_expr(q::QCommutator, coeff_fun::CFunction, expr::Vector{QExpr}) = QComposite[QCommutator(q.qspace, coeff_fun, expr)]
modify_coeff(q::QCommutator, coeff_fun::CFunction)::QCommutator = QCommutator(q.qspace, coeff_fun, q.expr)
each_term(q::QMultiComposite) = q.expr
each_coeff(q::QMultiComposite)::Vector{CFunction} = CFunction[q.coeff_fun; flatmap_to(each_coeff, each_term(q), CFunction)]
each_term(q::QComposite) = [q.expr]
each_coeff(q::QComposite)::Vector{CFunction} = CFunction[q.coeff_fun; each_coeff(q.expr)]

### Non-simple qFunctions 
"""
    QExp(qspace, coeff_fun, expr)

Composite for the symbolic exponential `coeff_fun * exp(expr)` used inside
`QExpr` terms without expanding the series.
"""
struct QExp <: QComposite
    qspace::QSpace
    coeff_fun::CFunction
    expr::QExpr
end
function _QExp(qspace::QSpace, coeff_fun::CFunction, expr::QExpr)::Vector{QComposite}
    if length(expr) == 1 && isa(expr[1], QLog)
        return QComposite[expr[1].expr]
    elseif isnumeric(expr) 
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
        return QComposite[modify_coeff(expr[1], coeff_fun * exp(sum_of_coeff_funs)) ]
    end
    return QComposite[QExp(qspace, coeff_fun, simplify_QExpr(expr))]
end
function _QExp(qspace::QSpace, coeff_fun::CFunction, expr::QExpr, ::Val{:nosimp})::Vector{QComposite}
    if length(expr) == 1 && isa(expr[1], QLog)
        return QComposite[expr[1].expr]
    elseif isnumeric(expr)
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
        return QComposite[modify_coeff(expr[1], coeff_fun * exp(sum_of_coeff_funs))]
    end
    return QComposite[QExp(qspace, coeff_fun, expr)]
end
function exp(q::QExpr)::QExpr
    return QExpr(q.qspace, _QExp(q.qspace, q.qspace.c_one, q))
end
modify_expr(q::QExp, expr::QExpr) = _QExp(q.qspace, q.coeff_fun, expr)
modify_expr(q::QExp, expr::QExpr, ::Val{:nosimp}) = _QExp(q.qspace, q.coeff_fun, expr, Val(:nosimp))
modify_coeff_expr(q::QExp, coeff_fun::CFunction, expr::QExpr) = _QExp(q.qspace, coeff_fun, expr)
modify_coeff(q::QExp, coeff_fun::CFunction) = QExp( q.qspace, coeff_fun, q.expr)
iszero(q::QExp) = iszero(q.coeff_fun) 

"""
    QLog(qspace, coeff_fun, expr)

Composite representing `coeff_fun * log(expr)`.
"""
struct QLog <: QComposite
    qspace::QSpace
    coeff_fun::CFunction
    expr::QExpr
end
function _QLog(qspace::QSpace, coeff_fun::CFunction, expr::QExpr)::Vector{QComposite}
    if length(expr) == 1 && isa(expr[1], QExp)
        return QComposite[expr[1].expr * coeff_fun]
    end
    if isnumeric(expr)
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
        return QComposite[modify_coeff(expr[1], coeff_fun * log(sum_of_coeff_funs)) ]
    end
    return QComposite[QLog(qspace, coeff_fun, simplify_QExpr(expr))]
end
function _QLog(qspace::QSpace, coeff_fun::CFunction, expr::QExpr, ::Val{:nosimp})::Vector{QComposite}
    if length(expr) == 1 && isa(expr[1], QExp)
        return QComposite[expr[1].expr * coeff_fun]
    end
    if isnumeric(expr)
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
        return QComposite[modify_coeff(expr[1], coeff_fun * log(sum_of_coeff_funs))]
    end
    return QComposite[QLog(qspace, coeff_fun, expr)]
end
function log(q::QExpr)::QExpr
    return QExpr(q.qspace, _QLog(q.qspace, q.qspace.c_one, q))
end
modify_expr(q::QLog, expr::QExpr) = _QLog(q.qspace, q.coeff_fun, expr)
modify_expr(q::QLog, expr::QExpr, ::Val{:nosimp}) = _QLog(q.qspace, q.coeff_fun, expr, Val(:nosimp))
modify_coeff_expr(q::QLog, coeff_fun::CFunction, expr::QExpr) = _QLog(q.qspace, coeff_fun, expr)
modify_coeff(q::QLog, coeff_fun::CFunction)::QLog = QLog(q.qspace, coeff_fun, q.expr)
iszero(q::QLog) = iszero(q.coeff_fun) #|| isone(q.expr) => that should be autosimplified

"""
    QPower(qspace, coeff_fun, n, expr)

Symbolic integer power `coeff_fun * expr^n` kept in unevaluated form.
"""
struct QPower <: QCompositeN
    qspace::QSpace
    coeff_fun::CFunction
    n::Int
    expr::QExpr
end
function _QPower(qspace::QSpace, coeff_fun::CFunction, n::Int, expr::QExpr)::Vector{QComposite}
    if isnumeric(expr)
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
        return QComposite[modify_coeff(expr[1], coeff_fun * power(sum_of_coeff_funs, n))]
    end
    return QComposite[QPower(qspace, coeff_fun, n, expr)]
end
function _QPower(qspace::QSpace, coeff_fun::CFunction, n::Int, expr::QExpr, ::Val{:nosimp})::Vector{QComposite}
    if isnumeric(expr)
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
        return QComposite[modify_coeff(expr[1], coeff_fun * power(sum_of_coeff_funs, n))]
    end
    return QComposite[QPower(qspace, coeff_fun, n, expr)]
end
modify_expr(q::QPower, expr::QExpr) = _QPower(q.qspace, q.coeff_fun, q.n, expr)
modify_expr(q::QPower, expr::QExpr, ::Val{:nosimp}) = _QPower(q.qspace, q.coeff_fun, q.n, expr, Val(:nosimp))
modify_coeff_expr(q::QPower, coeff_fun::CFunction, expr::QExpr)::QPower = _QPower(q.qspace, coeff_fun, q.n, expr)
modify_coeff(q::QPower, coeff_fun::CFunction) = QPower(q.qspace, coeff_fun, q.n, q.expr)

""" 
    power(q::QExpr, n::Int; force_symbolic::Bool=false)::QExpr

Returns `q^n` without expanding the expression. When `force_symbolic=true`, the
result is wrapped in a `QPower` composite; otherwise, if `q` consists of a
single `QAtomProduct`, the exponentiation falls back to the general `^`
operator (which multiplies the expression directly).
"""
function power(q::QExpr, n::Int; force_symbolic::Bool=false)::QExpr
    if n == 1
        return q
    end
    if !force_symbolic && length(q.terms) == 1 && isa(q.terms[1], QAtomProduct)
        return q ^ n
    end
    return QExpr(q.qspace, _QPower(q.qspace, q.qspace.c_one, n, q))
end


"""
    QRoot(qspace, coeff_fun, n, expr)

Symbolic n-th root `coeff_fun * expr^(1/n)` left unevaluated.
"""
struct QRoot <: QCompositeN
    qspace::QSpace
    coeff_fun::CFunction
    n::Int
    expr::QExpr
end
function _QRoot(qspace::QSpace, coeff_fun::CFunction, n::Int, expr::QExpr)::Vector{QComposite}
    if isnumeric(expr)
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
        return QComposite[modify_coeff(expr[1], coeff_fun * root(sum_of_coeff_funs, n))]
    end 
    return QComposite[QRoot(qspace, coeff_fun, n, expr)]
end
function _QRoot(qspace::QSpace, coeff_fun::CFunction, n::Int, expr::QExpr, ::Val{:nosimp})::Vector{QComposite}
    if isnumeric(expr)
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in expr)
        return QComposite[modify_coeff(expr[1], coeff_fun * root(sum_of_coeff_funs, n))]
    end
    return QComposite[QRoot(qspace, coeff_fun, n, expr)]
end
""" 
    root(q::QExpr, n::Int)::QExpr

Returns the n'th root of a quantum expression q^{1/n}, without expanding the expression.
"""
function root(q::QExpr, n::Int=2)::QExpr
    if n == 1 
        return q 
    end
    return QExpr(q.qspace, _QRoot(q.qspace, q.qspace.c_one, n, q))
end
function sqrt(q::QExpr)::QExpr
    return QExpr(q.qspace, _QRoot(q.qspace, q.qspace.c_one, 2, q))
end
modify_expr(q::QRoot, expr::QExpr) = _QRoot(q.qspace, q.coeff_fun, q.n, expr)
modify_expr(q::QRoot, expr::QExpr, ::Val{:nosimp}) = _QRoot(q.qspace, q.coeff_fun, q.n, expr, Val(:nosimp))
modify_coeff_expr(q::QRoot, coeff_fun::CFunction, expr::QExpr)::QRoot = QRoot(coeff_fun, q.n, expr)
modify_coeff(q::QRoot, coeff_fun::CFunction) = QRoot(q.qspace, coeff_fun, q.n, q.expr)

"""
    ExpectedValues(obj)

Return a copy of `obj` where every `QAtomProduct` is wrapped in an expectation value.
Works on individual products, composite expressions, differential equations and
collections thereof. Non-composite quantum objects are returned unchanged.
"""
Expectation(q::QAtomProduct) = q.braket ? q : set_braket(q, true)
Expectation(q::QExpr) = QExpr(q.qspace, [Expectation(term) for term in q.terms])
Expectation(q::QSum) = modify_expr(q, Expectation(q.expr), Val(:nodecollision))[1]
Expectation(q::QComposite) = modify_expr(q, Expectation(q.expr), Val(:nosimp))[1]
Expectation(q::QMultiComposite) = modify_expr(q, Expectation.(q.expr), Val(:nosimp))[1]
Expectation(q::QAtom) = q
