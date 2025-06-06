module qExpressions
using ..qSpace
using ..FFunctions
using ..StringUtils
using ComplexRationals
import Base: show, adjoint, conj, iterate, length, eltype, +, -, sort, *, ^, product, iszero, copy

export qObj, qAtom, qComposite, qMultiComp, qTerm, qAtomProduct, qExpr, qSum, Sum, ∑, diff_qEQ, term, base_operators,simplify, flatten, neq, d_dt


# ==========================================================================================================================================================
# --------> Base Types and Their Constructors <---------------------------------------------------------------------------------------------------------
# ==========================================================================================================================================================

""" 
    qObj

The abstract type `qObj` is the base type for all quantum expressions in this module.
"""
abstract type qObj end  # most general
""" 
    qAtom

The abstract type `qAtom` is a subtype of `qObj` and represents elementary operator definitions, such as qTerm and qAbstract. 
"""
abstract type qAtom <: qObj end # elementary operator definitions 
""" 
    qComposite

The abstract type `qComposite` is a subtype of `qObj` and represents composite expressions, 
such as qSum and qAtomProduct which consist of qAtom, qAbstract or qComposite objects themselves. 
"""
abstract type qComposite <: qObj end  # products and sums of operator definitions

"""
    qMultiComp

The abstract type `qMultiComposite` is a subtype of `qComposite` and represents composite expressions that contain multiple `qExpr` objects, such as `qCompositeProd` and `qCommutator`. 
"""
abstract type qMultiComposite <: qComposite end  

Is = Union{Int,Vector{Int}}

"""
    qTerm

A `qTerm` represents a single term in a quantum expression. It contains:
    - `op_indices`: A vector of indices representing the operators in the term, which are also defined in a StateSpace.
"""
mutable struct qTerm <: qAtom
    op_indices::Vector{Is}
end

"""
    qAbstract(indices::Vector{Int})

A purely‐symbolic abstract operator
    - key_index: The index of the abstract_key in the state space.
    - sub_index: The index of the suboperator in the state_space 
    - exponent: The exponent of the operator.
    - dag: A boolean indicating whether the operator is daggered (default = `false`)
    - operator_type: A reference to the operator of which it is a type 
    - index_map: Keeps track of indexes, that are equal (for neq transformations)
Is an instance of an OperatorType 
"""
mutable struct qAbstract <: qAtom
    key_index::Int
    sub_index::Int
    exponent::Int
    dag::Bool
    operator_type::OperatorType
    index_map::Vector{Tuple{Int,Int}}
end
function qAbstract(operator_type::OperatorType, key_index::Int, sub_index::Int=-1, exponent::Int=1, dag::Bool=false; index_map::Vector{Tuple{Int,Int}}=Tuple{Int,Int}[])
    return qAbstract(key_index, sub_index, exponent, dag, operator_type, index_map)
end

""" 
    qAtomProduct

A product of qAtom expressions, i.e. qTerms or qAbstract.
It contains:
    - `statespace`: The state space in which the product is defined.
    - `coeff_fun`: The function of parameters for the Operator product
    - `expr`: A vector of qAtoms (qTerms or qAbstract) that are multiplied together.
"""
mutable struct qAtomProduct <: qComposite
    statespace::StateSpace         # State space of the product.
    coeff_fun::FFunction            # function of scalar parameters => has +,-,*,/,^ defined 
    expr::Vector{qAtom}             # Vector of qAtoms (qTerms or qAbstract).
    function qAtomProduct(statespace::StateSpace, coeff::FFunction, expr::Vector{<:qAtom}= [])
        if isempty(expr)
           # add neutral operator
           expr = [statespace.neutral_op]
        end
        new(statespace, coeff, expr)
    end
    function qAtomProduct(statespace::StateSpace, coeff::Number, var_exponents::Vector{Int}, expr::qAtom)
        f_fun = FAtom(coeff, var_exponents) 
        return new(statespace, f_fun, [expr]) 
    end
    function qAtomProduct(statespace::StateSpace, coeff::Number, var_exponents::Vector{Int}, expr::AbstractVector{<:qAtom})
        f_fun = FAtom(coeff, var_exponents)
        return new(statespace, f_fun, expr)
    end
    function qAtomProduct(statespace::StateSpace, coeff::Number, expr::AbstractVector{<:qAtom})
        f_fun = coeff* statespace.fone
        return new(statespace, f_fun, expr)
    end
    function qAtomProduct(statespace::StateSpace, coeff::Number, expr::qAtom)
        f_fun = coeff * statespace.fone
        return new(statespace, f_fun, [expr]) 
    end
end

"""
    qExpr

A `qExpr` represents a quantum equation, consisting of a Vector of quantum Expressions representing the additive terms of the equation.
It also contains a reference to the state space in which the equation is defined.
"""
mutable struct qExpr
    statespace::StateSpace
    terms::Vector{qComposite}         # Vector of terms
    function qExpr(statespace::StateSpace, terms::Vector{<:qComposite})
        if isempty(terms) 
            # add neotral zero term
            zero_term = qAtomProduct(statespace.fone*0, qTerm(statespace.neutral_op))
            push!(terms, zero_term)
        end
        return new(statespace, terms)
    end
end
function qExpr(statespace::StateSpace, prod::qComposite)
    return qExpr(statespace, [prod])
end
function qExpr(statespace::StateSpace, terms::qAtom)
    return qExpr(statespace, qAtomProduct(statespace, statespace.fone, [terms]))
end

""" 
    qSum

A `qSum` represents the summation of a quantum Equation over indexes in a quantum expression.
It contains:
    - `expr`: The expression being summed over, which is a `qExpr` object.
    - `indexes`: A vector of strings representing the summation indexes (e.g., "i").
    - `subsystem_index`: The index of the subspace in which the indexes live. 
    - `element_indexes`: A vector of integers representing the position of the indexes in that subspace.
    - `neq`: A boolean indicating whether different indexes in the sum can refer to the same element in the subspace. 
            For example, the indexes i,j,k can refer to different elements in a much larger bath of elements. 
"""
mutable struct qSum <: qComposite
    statespace::StateSpace
    expr::qExpr       # The expression being summed over.    # use expr in other qComposites except for qAtomProduct
    indexes::Vector{String}   # The summation index (e.g. "i").
    subsystem_index::Int  # The subspace index where the summation index was found.
    element_indexes::Vector{Int}    # The position in that subspace.
    neq::Bool
end

"""
    diff_qEQ

A `diff_qEQ` represents a differential equation of the form:

    d/dt ⟨Op⟩ = RHS

It represents time evolution of operator expectation values, and wraps the symbolic structure of such an equation.

# Fields
- `left_hand_side::qTerm`: The LHS operator being differentiated.
- `expr::qExpr`: The RHS symbolic expression.
- `statespace::StateSpace`: The StateSpace in which the equation is defined.
- `braket::Bool`: Whether to use braket notation ⟨⋯⟩ (default = `true`).
- `do_sigma::Bool`: Whether to display Pauli operators as `σₓ`, etc. (default = `true`).
"""
mutable struct diff_qEQ <: qComposite
    left_hand_side::qAtomProduct
    expr::qExpr
    statespace::StateSpace
    braket::Bool
    do_sigma::Bool
end

"""
    diff_qEQ(lhs::qTerm, rhs::qExpr, statespace::StateSpace; braket=true, do_sigma=true)

Construct a [`diff_qEQ`](@ref) that represents the time derivative of ⟨lhs⟩ = rhs.

Automatically applies `neq()` to the RHS to expand sums over distinct indices.
"""
function diff_qEQ(left_hand_side::qAtomProduct, expr::qExpr, statespace::StateSpace; braket::Bool=true, do_sigma::Bool=false)
    new_rhs = neq(expr)
    return diff_qEQ(left_hand_side, new_rhs, statespace, braket, do_sigma)
end

"""
    Sum(index::Union{String,Symbol,Vector{String},Vector{Symbol}}, expr::qExpr; neq::Bool=false) -> qSum

Constructor of a `qSum` struct. Defines the indexes to sum over, the expressions for which to apply the sum and optionally whether the sum is only over non equal indexes. 
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
        return qExpr(ss, [qSum(ss, expr, index_strs, the_s_ind, e_inds, neq)])
    else
        return expr
    end
end
function Sum(index::Union{String,Symbol}, expr::qExpr; neq::Bool=false)::qExpr
    return Sum([index], expr, neq=neq)
end
""" 
    ∑(index::Union{String,Symbol}, expr::qExpr; neq::Bool=false) -> qSum

Alternative way to call the `Sum` constructor. Sum(index, expr; neq) = ∑(index, expr; neq).
"""
∑(index::Union{String,Symbol}, expr::qExpr; neq::Bool=false) = Sum(index, expr, neq=neq)
∑(indexes::Union{Vector{String},Vector{Symbol}}, expr::qExpr; neq::Bool=false) = Sum(indexes, expr, neq=neq)
 
#### Helper Functions #######################################################################################
# Define iteration for qExpr so that iterating over it yields its qTerm's.
function iterate(q::qExpr, state::Int=1)
    state > length(q.terms) && return nothing
    return (q.terms[state], state + 1)
end

# Optionally, define length and eltype.
length(q::qExpr) = length(q.terms)
length(q::qSum) = length(q.expr.terms)
length(q::qAtomProduct) = length(q.expr)
iszero(q::qAtomProduct) = iszero(q.coeff_fun)
iszero(q::qExpr) = length(q.terms) == 0 || all(iszero, q.terms)
iszero(q::qSum) = iszero(q.expr)

function copy(q::qTerm)::qTerm
    return qTerm(copy(q.op_indices))
end
function copy(q::qAbstract)::qAbstract
    return qAbstract(copy(q.key_index), copy(q.sub_index), copy(q.exponent), copy(q.dag), copy(q.operator_type), copy(q.index_map))
end
function copy(q::qAtomProduct)::qAtomProduct
    return qAtomProduct(q.statespace, copy(q.coeff_fun), copy(q.expr))
end
function copy(q::qSum)::qSum
    return qSum(copy(q.statespace), copy(q.expr), copy(q.indexes), copy(q.subsystem_index), copy(q.element_indexes), copy(q.neq))
end
function copy(q::qExpr)::qExpr
    return qExpr(q.statespace, copy(q.terms))
end

function copy(q::diff_qEQ)::diff_qEQ
    return diff_qEQ(copy(q.left_hand_side), copy(q.expr), copy(q.statespace), copy(q.braket), copy(q.do_sigma))
end


include("qExpressionsOps/qExpressions_Helper.jl") # Helper functions for qAtomProduct simplify

include("qExpressionsOps/qExpressions_string2term.jl")
include("qExpressionsOps/qExpressions_base_operators.jl")
include("qExpressionsOps/qExpressions_sorting.jl")
include("qExpressionsOps/qExpressions_simplify.jl")



include("qExpressionsOps/qExpressionsAlgebra.jl")
include("qExpressionsOps/qExpressionsPrint.jl")

include("qExpressionsOps/qSum_modify.jl")

#


"""
    d_dt(statespace::StateSpace, expr)

Evaluate the time derivative of an expression `expr` in the context of the given state space `ss`.

This function expects that `expr` is an equation (i.e. an Expr with an equal sign as its head),
of the form

    LHS = RHS
The function then returns a `diff_qEQ` constructed from the left-hand side qTerm and the right-hand side qExpr.
"""
function d_dt(left_hand::Union{qAtomProduct,qExpr}, right_hand::qExpr)::diff_qEQ
    # Check if expr is an equality.
    qstate = right_hand.statespace

    if left_hand isa qExpr
        if left_hand.statespace != qstate
            error("Left and right sides of the equation must be in the same state space.")
        end
        if length(left_hand.terms) != 1
            error("Left-hand side of the equation must consist of a single qTerm.")
        end
        left_hand = left_hand.terms[1]
    end
    if abs(left_hand.coeff - 1) > 1e-10
        error("Left-hand side of the equation must be a qTerm with coeff 1.")
    end
    if !iszero(left_hand.coeff_fun.var_exponents)
        error("Left-hand side of the equation must be a qTerm with no variable exponents.")
    end
    # Return a diff_qEQ constructed from these sides.
    return diff_qEQ(left_hand, right_hand, qstate)
end
end