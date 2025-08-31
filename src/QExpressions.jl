module QExpressions
using ..QSpace
using ..CFunctions
using ..StringUtils
using ComplexRationals
import Base: show, adjoint, conj, iterate, getindex, length, eltype, +, -, sort, *, ^, product, iszero, copy
using ..QAlgebra: FLIP_IF_FIRST_TERM_NEGATIVE, DO_BRACED
using ..CFunctions: isnumeric
export QEq, QObj, QAtom, QAbstract, QComposite, QCompositeN, QMultiComposite, QTerm, QExpr, diff_QEq, base_operators, d_dt #simplify

# ==========================================================================================================================================================
# --------> Base Types and Their Constructors <---------------------------------------------------------------------------------------------------------
# ==========================================================================================================================================================
# We are constructing terms and equations as an Abstract Syntax Tree 
Is = Vector{Int}

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
    QCompositeN

The abstract type `QCompositeN` is a subtype of `QComposite` and represents composite expressions, 
such as QPower and QRoot which have an additional element `n` with integer value.
"""
abstract type QCompositeN <: QComposite end  # QComposite with additional argument n

""" 
    QMultiComposite 

Abstract type for composite expressions that contain a Vector of QExpr objects.
"""
abstract type QMultiComposite <: QComposite end

""" 
    QEq

The abstract type `QEq` is the base type for all quantum expressions in this module.
"""
abstract type QEq end  # most general

"""
    QTerm

A `QTerm` represents a single term in a quantum expression. It contains:
    - `op_indices`: A vector of indices representing the operators in the term, which are also defined in a StateSpace.
"""
struct QTerm <: QAtom
    op_indices::Vector{Vector{Int}}
    function QTerm(op_indices::Vector{Vector{Int}})
        return new(copy.(op_indices))
    end
    function QTerm(op_indices::Vector{Vector{Int}}, ::Val{:nocopy})
        return new(op_indices)
    end
end 
@inline function Base.getindex(qterm::QTerm, i::Int)
    return qterm.op_indices[i]
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
struct QAbstract <: QAtom
    key_index::Int
    sub_index::Int
    exponent::Int
    dag::Bool
    operator_type::OperatorType
    index_map::Vector{Tuple{SubSpaceIndex,SubSpaceIndex}}
    function QAbstract(operator_type::OperatorType, key_index::Int, sub_index::Int=-1, exponent::Int=1, dag::Bool=false, index_map::Vector{Tuple{SubSpaceIndex,SubSpaceIndex}}=Tuple{SubSpaceIndex,SubSpaceIndex}[])
        return new(key_index, sub_index, exponent, dag, operator_type, index_map)
    end
end
dag_copy(q::QAbstract)::QAbstract = QAbstract(q.operator_type, q.key_index, q.sub_index, q.exponent, !q.dag, q.index_map)
add_to_index_map(q::QAbstract, added_index_pair::Tuple{SubSpaceIndex,SubSpaceIndex}) = QAbstract(q.operator_type, q.key_index, q.sub_index, q.exponent, q.dag, vcat(q.index_map, added_index_pair))
change_exp_dag(q::QAbstract, new_exp::Int, new_dag::Bool) = QAbstract(q.operator_type, q.key_index, q.sub_index, new_exp, new_dag, q.index_map)


"""
    QExpr

A `QExpr` represents a quantum equation, consisting of a Vector of quantum Expressions representing the additive terms of the equation.
It also contains a reference to the state space in which the equation is defined.
"""
struct QExpr <: QObj
    statespace::StateSpace
    terms::Vector{QComposite}              #AbstractVector{<:QComposite}    
    function QExpr(statespace::StateSpace, terms::AbstractVector{<:QComposite})
        if isempty(terms) 
            # add neotral zero term
            zero_term = QAtomProduct(statespace, statespace.c_zero, QAtom[])
            terms = [zero_term]
        end
        return new(statespace, terms)
    end
    function QExpr(statespace::StateSpace, terms::AbstractVector{<:QComposite}, ::Val{:simp})
        if isempty(terms) 
            # add neotral zero term
            zero_term = QAtomProduct(statespace, statespace.c_zero, QAtom[])
            terms = [zero_term]
        end
        return new(statespace, simplify_QExpr(Vector{QComposite}(terms)))
    end
    function QExpr(terms::AbstractVector{<:QComposite})
        return new(terms[1].statespace, copy(terms))
    end
    function QExpr(statespace::StateSpace, prod::T) where T<:QComposite
        return new(statespace, QComposite[prod])
    end
    function QExpr(statespace::StateSpace, terms::QAtom)
        return new(statespace, QComposite[QAtomProduct(statespace,terms)])
    end
    function QExpr(terms::AbstractVector{<:QComposite}, ::Val{:simp})
        return new(terms[1].statespace, simplify_QExpr(Vector{QComposite}(terms)))
    end
end
copy(q::QExpr)::QExpr = QExpr(q.statespace, q.terms)
length(q::QExpr) = length(q.terms)
each_term(q::QExpr) = q.terms
each_coeff(q::QExpr)::Vector{CFunction} = flatmap_to(each_coeff, each_term(q), CFunction)
multiply_coeff(q::QExpr, coeff::CFunction) = QExpr(q.statespace, [multiply_coeff(s, coeff) for s in q.terms])

include("QExpressionsOps/QExpressions_composites.jl")
include("QExpressionsOps/QExpressions_helper.jl") 


"""
    diff_QEq

A `diff_QEq` represents a differential equation of the form:

    d/dt ⟨Op⟩ = RHS

It represents time derivative of an operator expectation value, and wraps the symbolic structure of such an equation.

# Fields
- `left_hand_side::QTerm`: The LHS operator being differentiated.
- `expr::QExpr`: The RHS symbolic expression.
- `statespace::StateSpace`: The StateSpace in which the equation is defined.
- `braket::Bool`: Whether to use braket notation ⟨⋯⟩ (default = `true`).
"""
struct diff_QEq <: QEq
    statespace::StateSpace
    left_hand_side::QAtomProduct
    expr::QExpr 
    braket::Bool
end
copy(q::diff_QEq)::diff_QEq = diff_QEq(q.statespace, copy(q.left_hand_side), copy(q.expr), q.braket)

"""
    diff_QEq(lhs::QTerm, rhs::QExpr, statespace::StateSpace; braket=true)

Construct a [`diff_QEq`](@ref) that represents the time derivative of ⟨lhs⟩ = rhs.

Automatically applies `neq()` to the RHS to expand sums over distinct indices.
"""
function diff_QEq(statespace::StateSpace, left_hand_side::QAtomProduct, expr::QExpr; braket::Bool=true)
    if !contains_abstract(left_hand_side)
        where_acting = which_ensemble_acting(left_hand_side)
        new_rhs = neq(expr, where_acting)
        return diff_QEq(statespace, left_hand_side, new_rhs, braket)
    else
        return diff_QEq(statespace, left_hand_side, expr, braket)
    end
end


#### Helper Functions #######################################################################################
# Define iteration for QExpr so that iterating over it yields its QTerm's.
function iterate(q::QExpr, state::Int=1)
    state > length(q.terms) && return nothing
    return q.terms[state], state + 1
end
function iterate(q::T, state::Int=1) where T <: QComposite
    error("Cannot iterate over a QComposite of type $(T).")
end

function getindex(q::QExpr, i::Int)
    q.terms[i]
end
function getindex(q::T, i::Int) where T <: QComposite
    q.expr[i]
end

# Optionally, define length and eltype.
iszero(q::QExpr) = length(q.terms) == 0 || all(iszero, q.terms)
iszero(q::QAtomProduct) = iszero(q.coeff_fun)
iszero(q::QSum) = iszero(q.expr)
iszero(q::T) where T<:QComposite = iszero(q.coeff_fun) || iszero(q.expr)
iszero(q::T) where T<:QMultiComposite = iszero(q.coeff_fun) || any(iszero, q.expr) 

include("QExpressionsOps/QExpressions_base_operators.jl")
include("QExpressionsOps/QExpressions_sort.jl")
include("QExpressionsOps/QExpressions_simplify.jl")

include("QExpressionsOps/QExpressions_algebra.jl")
include("QExpressionsOps/QExpressions_print.jl")

include("QExpressionsOps/QSum_modify.jl")

include("QExpressionsOps/QExpressions_welldefined.jl")
include("QExpressionsOps/QExpressions_substitute.jl")
include("QExpressionsOps/QExpressions_repartition.jl")

include("QExpressionsOps/QExpressions_cumulants.jl")


"""
    d_dt(statespace::StateSpace, expr)

Evaluate the time derivative of an expression `expr` in the context of the given state space `ss`.

This function expects that `expr` is an equation (i.e. an Expr with an equal sign as its head),
of the form

    LHS = RHS
The function then returns a `diff_QEq` constructed from the left-hand side QTerm and the right-hand side QExpr.
"""
function d_dt(left_hand::Union{QAtomProduct,QExpr}, right_hand::QExpr)::diff_QEq
    # Check if expr is an equality.
    qstate = right_hand.statespace

    if left_hand isa QExpr
        if left_hand.statespace != qstate
            error("Left and right sides of the equation must be in the same state space.")
        end
        if length(left_hand.terms) != 1
            error("Left-hand side of the equation must consist of a single QTerm.")
        end
        left_hand = left_hand.terms[1]
        if !isa(left_hand, QAtomProduct)
            error("Left-hand side of the equation must be a QAtomProduct. Or a QAtomProduct wrapped in a QExpr. ")
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