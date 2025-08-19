module QExpressions
using ..QSpace
using ..CFunctions
using ..StringUtils
using ComplexRationals
import Base: show, adjoint, conj, iterate, getindex, length, eltype, +, -, sort, *, ^, product, iszero, copy
using ..QAlgebra: FLIP_IF_FIRST_TERM_NEGATIVE, DO_BRACED
using ..CFunctions: isnumeric
export QEq, QObj, QAtom, QAbstract, QComposite, QCompositeN, QCompositeProduct, QMultiComposite, QTerm, QAtomProduct, QExpr, QSum, Sum, ∑, diff_QEq, base_operators, simplify, simplify_QAtomProduct, flatten, neq, d_dt

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
    index_map::Vector{Tuple{Int,Int}}
    function QAbstract(operator_type::OperatorType, key_index::Int, sub_index::Int=-1, exponent::Int=1, dag::Bool=false; index_map::Vector{Tuple{Int,Int}}=Tuple{Int,Int}[])
        return new(key_index, sub_index, exponent, dag, operator_type, copy.(index_map))
    end
end
dag_copy(q::QAbstract)::QAbstract = QAbstract(q.operator_type, q.key_index, q.sub_index, q.exponent, !q.dag, index_map=q.index_map)
add_to_index_map(q::QAbstract, added_index_pair::Tuple{Int,Int}) = QAbstract(q.operator_type, q.key_index, q.sub_index, q.exponent, q.dag, vcat(q.index_map, added_index_pair))

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
modify_expr(q::QAtomProduct, expr::Vector{QAtom})::QAtomProduct = QAtomProduct(q.statespace, q.coeff, expr, q.separate_expectation_values, Val(:dont_check_time))
modify_coeff_expr(q::QAtomProduct, coeff::CFunction, expr::Vector{QAtom})::QAtomProduct = QAtomProduct(q.statespace, coeff, expr, q.separate_expectation_values)
modify_coeff(q::QAtomProduct, coeff::CFunction)::QAtomProduct = QAtomProduct(q.statespace, coeff, q.expr, q.separate_expectation_values)
each_term(q::QAtomProduct) = q.expr
each_coeff(q::QAtomProduct)::Vector{CFunction} = [q.coeff]


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
            zero_term = QAtomProduct(statespace, statespace.c_one*0, QAtom[])
            terms = [zero_term]
        end
        return new(statespace, terms)
    end
    function QExpr(statespace::StateSpace, terms::AbstractVector{<:QComposite}, ::Val{:simp})
        if isempty(terms) 
            # add neotral zero term
            zero_term = QAtomProduct(statespace, statespace.c_one*0, QAtom[])
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
    indexes::Vector{String}   # The summation index (e.g. "i").
    subsystem_index::Int  # The subspace index where the summation index was found.
    element_indexes::Vector{Int}    # The position in that subspace.
    neq::Bool
    function QSum(statespace::StateSpace, expr::QExpr, indexes::Vector{String}, subsystem_index::Int, element_indexes::Vector{Int}, neq::Bool)
        return new(statespace,  expr, copy(indexes), subsystem_index, copy(element_indexes), neq)
    end
    function QSum(statespace::StateSpace, expr::QExpr, indexes::Vector{String}, subsystem_index::Int, element_indexes::Vector{Int}, neq::Bool, ::Val{:simp})
        return new(statespace,  expr, copy(indexes), subsystem_index, copy(element_indexes), neq)
    end
end
copy(q::QSum)::QSum = QSum(q.statespace,q.expr, q.indexes, q.subsystem_index, q.element_indexes, q.neq)
modify_expr(q::QSum, expr::QExpr) = QSum(q.statespace, expr, q.indexes, q.subsystem_index, q.element_indexes, q.neq)
modify_expr_indexes(q::QSum, expr::QExpr, indexes::Vector{String}, subsystem_index::Int, element_indexes::Vector{Int}) = QSum(q.statespace, expr, indexes, subsystem_index, element_indexes, q.neq)
each_term(q::QSum) = q.expr
each_coeff(q::QExpr)::Vector{CFunction} = flatmap_to(each_coeff, each_term(q), CFunction)

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
    new_rhs = neq(expr)
    return diff_QEq(statespace, left_hand_side, new_rhs, braket)
end

"""
    Sum(index::Union{String,Symbol,Vector{String},Vector{Symbol}}, expr::QExpr; neq::Bool=false) -> QSum

Constructor of a `QSum` struct. Defines the indexes to sum over, the expressions for which to apply the sum and optionally whether the sum is only over non equal indexes. 
"""
function Sum(indexes::Union{Vector{String},Vector{Symbol}}, expr::QExpr; neq::Bool=false)::QExpr
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
        sort!(e_inds)
        return QExpr(ss, [QSum(ss, expr, index_strs, the_s_ind, e_inds, neq, Val(:simp))])
    else
        return expr
    end
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

include("QExpressionsOps/QExpressions_functions.jl")
include("QExpressionsOps/QExpressions_helper.jl") # Helper functions for QAtomProduct simplify

include("QExpressionsOps/QExpressions_base_operators.jl")
include("QExpressionsOps/QExpressions_sort.jl")
include("QExpressionsOps/QExpressions_simplify.jl")

include("QExpressionsOps/QExpressions_string2term.jl")
include("QExpressionsOps/QExpressions_algebra.jl")
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