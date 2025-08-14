module QExpressions
using ..QSpace
using ..FFunctions
using ..StringUtils
using ComplexRationals
import Base: show, adjoint, conj, iterate, getindex, length, eltype, +, -, sort, *, ^, product, iszero, copy
using ..QAlgebra: get_default 
using ..FFunctions: isnumeric
export QObj, QAtom, QAbstract, QComposite, QCompositeProduct, QMultiComposite, QTerm, QAtomProduct, qExpr, QSum, Sum, ∑, diff_QEq, QEq, base_operators, simplify, simplifyqAtomProduct, flatten, neq, d_dt

# ==========================================================================================================================================================
# --------> Base Types and Their Constructors <---------------------------------------------------------------------------------------------------------
# ==========================================================================================================================================================
# We are constructing terms and equations as an Abstract Syntax Tree 
Is = Union{Int,Vector{Int}}

""" 
    QObj

The abstract type `QObj` is the base type for all quantum expressions in this module.
"""
abstract type QObj end  # most general
""" 
    QAtom

The abstract type `QAtom` is a subtype of `QObj` and represents elementary operator definitions, such as QTerm and QAbstract. 
"""
abstract type QAtom <: QObj end # elementary operator definitions 
""" 
    QComposite

The abstract type `QComposite` is a subtype of `QObj` and represents composite expressions, 
such as QSum and QAtomProduct which consist of QAtom, QAbstract or QComposite objects themselves. 
"""
abstract type QComposite <: QObj end  # products and sums of operator definitions

""" 
    QMultiComposite 

Abstract type for composite expressions that contain a Vector of qExpr objects.
"""
abstract type QMultiComposite <: QComposite end

"""
    QTerm

A `QTerm` represents a single term in a quantum expression. It contains:
    - `op_indices`: A vector of indices representing the operators in the term, which are also defined in a StateSpace.
"""
struct QTerm <: QAtom
    statespace::StateSpace
    op_indices::Vector{Is}
    function QTerm(statespace::StateSpace, op_indices::Vector{Is})
        new(statespace, op_indices)
    end
end
function copy(q::QTerm)::QTerm
    return QTerm(q.statespace, copy(q.op_indices))
end

"""
    QAbstract(indices::Vector{Int})

A purely‐symbolic abstract operator
    - key_index: The index of the abstract_key in the state space.
    - sub_index: The index of the suboperator in the state_space 
    - exponent: The exponent of the operator.
    - dag: A boolean indicating whether the operator is daggered (default = `false`)
    - operator_type: A reference to the operator of which it is a type 
    - index_map: Keeps track of indexes, that are equal (for neq transformations)
Is an instance of an OperatorType 
"""
mutable struct QAbstract <: QAtom
    statespace::StateSpace
    key_index::Int
    sub_index::Int
    exponent::Int
    dag::Bool
    operator_type::OperatorType
    index_map::Vector{Tuple{Int,Int}}
end
function QAbstract(statespace::StateSpace, operator_type::OperatorType, key_index::Int, sub_index::Int=-1, exponent::Int=1, dag::Bool=false; index_map::Vector{Tuple{Int,Int}}=Tuple{Int,Int}[])
    return QAbstract(statespace, key_index, sub_index, exponent, dag, operator_type, index_map)
end
function copy(q::QAbstract)::QAbstract
    return QAbstract(q.statespace, q.operator_type, q.key_index, q.sub_index, q.exponent, q.dag, index_map=copy(q.index_map))
end


""" 
    QAtomProduct

A product of QAtom expressions, i.e. qTerms or QAbstract.
It contains:
    - `statespace`: The state space in which the product is defined.
    - `coeff_fun`: The function of parameters for the Operator product
    - `expr`: A vector of qAtoms (qTerms or QAbstract) that are multiplied together.
"""
mutable struct QAtomProduct <: QComposite
    statespace::StateSpace         # State space of the product.
    coeff_fun::FFunction            # function of scalar parameters => has +,-,*,/,^ defined 
    expr::Vector{QAtom}             # Vector of qAtoms (qTerms or QAbstract).
    separate_expectation_values::Bool 
    function QAtomProduct(statespace::StateSpace, coeff::FFunction, expr::AbstractVector{<:QAtom}= QAtom[], separate_expectation_values::Bool=false)
        new(statespace, coeff, expr, separate_expectation_values)
    end
    function QAtomProduct(statespace::StateSpace, coeff::Number, var_exponents::Vector{Int}, expr::QAtom, separate_expectation_values::Bool=false)
        f_fun = FAtom(coeff, var_exponents)
        return new(statespace, f_fun, [expr], separate_expectation_values)
    end
    function QAtomProduct(statespace::StateSpace, coeff::Number, var_exponents::Vector{Int}, expr::AbstractVector{<:QAtom}, separate_expectation_values::Bool=false)
        f_fun = FAtom(coeff, var_exponents)
        return new(statespace, f_fun, expr, separate_expectation_values)
    end
    function QAtomProduct(statespace::StateSpace, coeff::Number, expr::AbstractVector{<:QAtom}, separate_expectation_values::Bool=false)
        f_fun = coeff* statespace.fone
        return new(statespace, f_fun, expr, separate_expectation_values)
    end
    function QAtomProduct(statespace::StateSpace, coeff::Number, expr::QAtom, separate_expectation_values::Bool=false)
        f_fun = coeff * statespace.fone
        return new(statespace, f_fun, [expr], separate_expectation_values)
    end
end
function copy(q::QAtomProduct)::QAtomProduct
    return QAtomProduct(q.statespace, copy(q.coeff_fun), [copy(s) for s in q.expr], q.separate_expectation_values)
end

"""
    qExpr

A `qExpr` represents a quantum equation, consisting of a Vector of quantum Expressions representing the additive terms of the equation.
It also contains a reference to the state space in which the equation is defined.
"""
mutable struct qExpr <: QObj
    statespace::StateSpace
    terms::Vector{QObj}              #AbstractVector{<:QObj}
    function qExpr(statespace::StateSpace, terms::AbstractVector{<:QObj})
        if isempty(terms) 
            # add neotral zero term
            zero_term = QAtomProduct(statespace, statespace.fone*0, QAtom[])
            terms = [zero_term]
        end
        return new(statespace, terms)
    end
end
function qExpr(term::T) where T<:QObj
    return qExpr(term.statespace, [term])
end
function qExpr(statespace::StateSpace, prod::T) where T<:QObj
    return qExpr(statespace, [prod])
end
function qExpr(statespace::StateSpace, terms::QAtom)
    return qExpr(statespace, QAtomProduct(statespace, statespace.fone, [terms]))
end
function qExpr(terms::AbstractVector{<:QObj})
    return qExpr(terms[1].statespace, terms)
end
function copy(q::qExpr)::qExpr
    return qExpr(q.statespace, [copy(s) for s in q.terms])
end

"""
    QSum

A `QSum` represents the summation of a quantum Equation over indexes in a quantum expression.
It contains:
    - `expr`: The expression being summed over, which is a `qExpr` object.
    - `indexes`: A vector of strings representing the summation indexes (e.g., "i").
    - `subsystem_index`: The index of the subspace in which the indexes live. 
    - `element_indexes`: A vector of integers representing the position of the indexes in that subspace.
    - `neq`: A boolean indicating whether different indexes in the sum can refer to the same element in the subspace. 
            For example, the indexes i,j,k can refer to different elements in a much larger bath of elements. 
"""
mutable struct QSum <: QComposite
    statespace::StateSpace
    expr::qExpr       # The expression being summed over.    # use expr in other qComposites except for QAtomProduct
    indexes::Vector{String}   # The summation index (e.g. "i").
    subsystem_index::Int  # The subspace index where the summation index was found.
    element_indexes::Vector{Int}    # The position in that subspace.
    neq::Bool
end
function copy(q::QSum)::QSum
    return QSum(q.statespace, copy(q.expr), copy(q.indexes), q.subsystem_index, copy(q.element_indexes), q.neq)
end

"""
    diff_QEq

A `diff_QEq` represents a differential equation of the form:

    d/dt ⟨Op⟩ = RHS

It represents time evolution of operator expectation values, and wraps the symbolic structure of such an equation.

# Fields
- `left_hand_side::QTerm`: The LHS operator being differentiated.
- `expr::qExpr`: The RHS symbolic expression.
- `statespace::StateSpace`: The StateSpace in which the equation is defined.
- `braket::Bool`: Whether to use braket notation ⟨⋯⟩ (default = `true`).
"""
abstract type QEq <: QObj end

struct diff_QEq <: QEq
    statespace::StateSpace
    left_hand_side::QAtomProduct
    expr::qExpr
    braket::Bool
end
function copy(q::diff_QEq)::diff_QEq
    return diff_QEq(q.statespace, copy(q.left_hand_side), copy(q.expr), q.braket)
end

"""
    diff_QEq(lhs::QTerm, rhs::qExpr, statespace::StateSpace; braket=true)

Construct a [`diff_QEq`](@ref) that represents the time derivative of ⟨lhs⟩ = rhs.

Automatically applies `neq()` to the RHS to expand sums over distinct indices.
"""
function diff_QEq(statespace::StateSpace, left_hand_side::QAtomProduct, expr::qExpr; braket::Bool=true)
    new_rhs = neq(expr)
    return diff_QEq(statespace, left_hand_side, new_rhs, braket)
end

"""
    Sum(index::Union{String,Symbol,Vector{String},Vector{Symbol}}, expr::qExpr; neq::Bool=false) -> QSum

Constructor of a `QSum` struct. Defines the indexes to sum over, the expressions for which to apply the sum and optionally whether the sum is only over non equal indexes. 
"""
function Sum(indexes::Union{Vector{String},Vector{Symbol}}, expr::qExpr; neq::Bool=false)::qExpr
    index_strs = [string(index) for index in indexes]
    ss = expr.statespace
    the_s_ind::Int = -1
    e_inds::Vector{Int} = []
    for index_str in index_strs
        found = false
        for (s_ind, sub) in enumerate(ss.subspaces)
            for (e_ind, key) in enumerate(sub.keys)
                if key == index_str
                    if the_s_ind == -1
                        the_s_ind = s_ind
                    else
                        if s_ind != the_s_ind
                            error("Index $index_str found in multiple subspaces. Please specify a single subspace.")
                        end
                    end
                    push!(e_inds, e_ind)
                    found = true
                    break
                end
            end
            if found
                break
            end
        end
        if !found
            error("Index $index_str not found in any subspace keys in the state space.")
        end
    end
    if length(index_strs) > 0
        if length(e_inds) != length(unique(e_inds))
            error("Duplicate indexes found in the input.")
        end
        e_inds = sort(e_inds)
        return qExpr(ss, [QSum(ss, expr, index_strs, the_s_ind, e_inds, neq)])
    else
        return expr
    end
end
function Sum(index::Union{String,Symbol}, expr::qExpr; neq::Bool=false)::qExpr
    return Sum([index], expr, neq=neq)
end
""" 
    ∑(index::Union{String,Symbol}, expr::qExpr; neq::Bool=false) -> QSum

Alternative way to call the `Sum` constructor. Sum(index, expr; neq) = ∑(index, expr; neq).
"""
∑(index::Union{String,Symbol}, expr::qExpr; neq::Bool=false) = Sum(index, expr, neq=neq)
∑(indexes::Union{Vector{String},Vector{Symbol}}, expr::qExpr; neq::Bool=false) = Sum(indexes, expr, neq=neq)
 
#### Helper Functions #######################################################################################
# Define iteration for qExpr so that iterating over it yields its QTerm's.
function iterate(q::qExpr, state::Int=1)
    state > length(q.terms) && return nothing
    return q.terms[state], state + 1
end
function iterate(q::T, state::Int=1) where T <: QComposite
    error("Cannot iterate over a QComposite of type $(T).")
end

function getindex(q::qExpr, i::Int)
    q.terms[i]
end
function getindex(q::T, i::Int) where T <: QComposite
    q.expr[i]
end

# Optionally, define length and eltype.
length(q::qExpr) = length(q.terms)

iszero(q::qExpr) = length(q.terms) == 0 || all(iszero, q.terms)
iszero(q::QAtomProduct) = iszero(q.coeff_fun)
iszero(q::T) where T<:QComposite = iszero(q.expr)
iszero(q::T) where T<:QMultiComposite = any(iszero, q.expr)



include("QExpressionsOps/QExpressions_functions.jl")
include("QExpressionsOps/QExpressions_helper.jl") # Helper functions for QAtomProduct simplify

include("QExpressionsOps/QExpressions_string2term.jl")
include("QExpressionsOps/QExpressions_base_operators.jl")
include("QExpressionsOps/QExpressions_sorting.jl")
include("QExpressionsOps/QExpressions_simplify.jl")



include("QExpressionsOps/QExpressionsAlgebra.jl")
include("QExpressionsOps/QExpressionsPrint.jl")

include("QExpressionsOps/QSum_modify.jl")

include("QExpressionsOps/QExpressions_welldefined.jl")
include("QExpressionsOps/QExpressions_substitute.jl")
include("QExpressionsOps/QExpressions_reorder.jl")

include("QExpressionsOps/QExpressions_cumulants.jl")


"""
    d_dt(statespace::StateSpace, expr)

Evaluate the time derivative of an expression `expr` in the context of the given state space `ss`.

This function expects that `expr` is an equation (i.e. an Expr with an equal sign as its head),
of the form

    LHS = RHS
The function then returns a `diff_QEq` constructed from the left-hand side QTerm and the right-hand side qExpr.
"""
function d_dt(left_hand::Union{QAtomProduct,qExpr}, right_hand::qExpr)::diff_QEq
    # Check if expr is an equality.
    qstate = right_hand.statespace

    if left_hand isa qExpr
        if left_hand.statespace != qstate
            error("Left and right sides of the equation must be in the same state space.")
        end
        if length(left_hand.terms) != 1
            error("Left-hand side of the equation must consist of a single QTerm.")
        end
        left_hand = left_hand.terms[1]
        if !isa(left_hand, QAtomProduct)
            error("Left-hand side of the equation must be a QAtomProduct. Or a QAtomProduct wrapped in a qExpr. ")
        end
    end
    if !isnumeric(left_hand.coeff_fun) && abs(left_hand.coeff_fun - 1) != 0
        error("Left-hand side of the equation must be a QTerm with a purely numeric coefficient of 1.")
    end
    if !iszero(left_hand.coeff_fun.var_exponents)
        error("Left-hand side of the equation must be a QTerm with no variable exponents.")
    end
    # Return a diff_QEq constructed from these sides.
    return diff_QEq(qstate, left_hand, right_hand)
end
end