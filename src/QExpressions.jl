module QExpressions
using ..QSpaces
using ..CFunctions
import ..CFunctions: expand
using ..StringUtils
using ComplexRationals
using SparseArrays
import Base: show, adjoint, conj, iterate, getindex, length, eltype, +, -, sort, *, /, ^, product, iszero, copy
using ..QAlgebra: FLIP_IF_FIRST_TERM_NEGATIVE, DO_BRACED, EXPAND_CUMULANTS, vecvec_or, vecvec_or!, findfirstfreeafterbefore, sorted_unique_push!, bubble_insert_unique!
using ..CFunctions: isnumeric
export QObj, QAtom, QAbstract, QComposite, QCompositeN, QMultiComposite, QTerm, QExpr, QCumulant, diffQEq, Expectation, base_operators, d_dt

export @define, @define_basics, QExpr2CFunction

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
such as `AbstractQSum` (and its aliases `QSum`/`QIntegral`) and `QAtomProduct`, which consist of `QAtom`,
`QAbstract` or `QComposite` objects themselves.
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

""" QParent

An abstract type to store parents of other quantum types, such as `QEq`s and `diffQEq`s.
""" 
abstract type QParent <: QObj end 

"""
    QTerm

A `QTerm` represents a single term in a quantum expression. It contains:
    - `op_indices`: A vector of indices representing the operators in the term, which are also defined in a QSpace.
"""
struct QTerm <: QAtom
    op_indices::Vector{Is}
    time_index::Int 
    function QTerm(op_indices::Vector{Is}, time_index::Int=-1)
        return new(copy.(op_indices), time_index)
    end
end 
@inline function Base.getindex(qterm::QTerm, i::Int)
    return qterm.op_indices[i]
end
modify_expr(q::QTerm, new_op_indices::Vector{Is}) = QTerm(new_op_indices, q.time_index)
function modify_time_index(q::QTerm, new_time_index::Int)::QTerm
    @assert q.time_index != -1 "Cannot change time_index of non time dependent QTerm."
    QTerm(q.op_indices, new_time_index)
end
of_time(q::QTerm) = q.time_index != -1


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
    time_index::Int
    operator_type::OperatorType
    index_map::Vector{Tuple{SubSpaceIndex,SubSpaceIndex}}
    function QAbstract(::Val{:nocheck}, operator_type::OperatorType, key_index::Int, sub_index::Int=-1, exponent::Int=1, dag::Bool=false, time_index::Int=-1, index_map::Vector{Tuple{SubSpaceIndex,SubSpaceIndex}}=Tuple{SubSpaceIndex,SubSpaceIndex}[])
        return new(key_index, sub_index, exponent, dag, time_index, operator_type, index_map)
    end
    function QAbstract(operator_type::OperatorType, key_index::Int, sub_index::Int=-1, exponent::Int=1, dag::Bool=false, time_index::Int=-1, index_map::Vector{Tuple{SubSpaceIndex,SubSpaceIndex}}=Tuple{SubSpaceIndex,SubSpaceIndex}[])
        if operator_type.of_time && time_index < 0
            error("Constructing an operator of time but no time index >= 0 has been provided!")
        end
        return new(key_index, sub_index, exponent, dag, time_index, operator_type, index_map)
    end
end
dag_copy(q::QAbstract)::QAbstract = QAbstract(q.operator_type, q.key_index, q.sub_index, q.exponent, !q.dag, q.time_index, q.index_map)
add_to_index_map(q::QAbstract, added_index_pair::Tuple{SubSpaceIndex,SubSpaceIndex}) = QAbstract(q.operator_type, q.key_index, q.sub_index, q.exponent, q.dag, q.time_index, vcat(q.index_map, added_index_pair))
add_to_index_map(q::QAbstract, added_index_pair::Vector{Tuple{SubSpaceIndex,SubSpaceIndex}}) = QAbstract(q.operator_type, q.key_index, q.sub_index, q.exponent, q.dag, q.time_index, vcat(q.index_map, added_index_pair))
modify_exp_dag(q::QAbstract, new_exp::Int, new_dag::Bool) = QAbstract(q.operator_type, q.key_index, q.sub_index, new_exp, new_dag, q.time_index, q.index_map)
modify_time_index(q::QAbstract, new_time_index::Int) = QAbstract(q.operator_type, q.key_index, q.sub_index, q.exponent, q.dag, new_time_index, q.index_map)
of_time(q::QAbstract) = q.operator_type.of_time

#### Reference the AST upwards aswell #################################
const ParentRef  = Union{Nothing, WeakRef}
const ParentCell = Base.RefValue{ParentRef}
"""
    QExpr

A `QExpr` represents a quantum equation, consisting of a Vector of quantum Expressions representing the additive terms of the equation.
It also contains a reference to the state space in which the equation is defined.
"""
struct QExpr <: QParent
    qspace::QSpace
    terms::Vector{QComposite}
    function QExpr(qspace::QSpace, terms::AbstractVector{<:QComposite}, ::Val{:nosimp})
        vec_terms = Vector{QComposite}(terms)
        if isempty(vec_terms)
            zero_term = QAtomProduct(qspace, qspace.c_zero, QAtom[])
            vec_terms = QComposite[zero_term]
        end
        new(qspace, vec_terms)
    end
    function QExpr(qspace::QSpace, terms::AbstractVector{<:QComposite})
        simplified = simplify_QExpr(Vector{QComposite}(terms))
        return QExpr(qspace, simplified, Val(:nosimp))
    end
    function QExpr(qspace::QSpace, prod::T) where T<:QComposite
        return QExpr(qspace, QComposite[prod], Val(:nosimp))
    end
    function QExpr(qspace::QSpace, term::QAtom)
        return QExpr(qspace, QComposite[QAtomProduct(qspace, term)], Val(:nosimp))
    end
end
length(q::QExpr) = length(q.terms)
each_term(q::QExpr) = q.terms
each_coeff(q::QExpr)::Vector{CFunction} = flatmap_to(each_coeff, each_term(q), CFunction)
multiply_coeff(q::QExpr, coeff::CFunction) = QExpr(q.qspace, [multiply_coeff(s, coeff) for s in q.terms])

include("QExpressionsOps/QExpressions_composites.jl")
include("QExpressionsOps/QExpressions_helper.jl") 
include("QExpressionsOps/QExpressions_expand.jl")
include("QExpressionsOps/QExpressions_cumulants.jl")

"""
    diffQEq

Represents a differential equation of the form `d/dt ⟨Op⟩ = RHS`, storing the left-hand
side operator product and right-hand side expression in a common `QSpace`.
"""
struct diffQEq <: QParent
    qspace::QSpace
    left_hand_side::QAtomProduct
    expr::QExpr 
    function diffQEq(qspace::QSpace, left_hand_side::QAtomProduct, expr::QExpr, ::Val{:raw})
        new(qspace, left_hand_side, expr)
    end
end


"""
    diffQEq(lhs::QTerm, rhs::QExpr, qspace::QSpace)

Construct a [`diffQEq`](@ref) that represents the time derivative of ⟨lhs⟩ = rhs.

Automatically applies `neq()` to the RHS to expand sums over distinct indices.
"""
function _normalize_diff_time(qspace::QSpace, lhs::QAtomProduct, rhs::QExpr)
    lhs_usage = contains_which_t_indexes(lhs)
    rhs_usage = contains_which_t_indexes(rhs)
    combined = lhs_usage .| rhs_usage
    active_positions = findall(combined)
    if length(active_positions) > 1
        labels = [pos == 1 ? :t : Symbol("t$(pos - 1)") for pos in active_positions]
        error("diffQEq requires at most one time index; found $(join(string.(labels), ", ")).")
    elseif length(active_positions) == 1
        idx = active_positions[1] - 1
        if idx != 0
            sp = Substitution_t(idx, 0)
            lhs_expr = substitute(QExpr(qspace, lhs), sp)
            length(lhs_expr.terms) == 1 || error("Time substitution changed LHS arity in diffQEq.")
            lhs_term = lhs_expr.terms[1]
            lhs_term isa QAtomProduct || error("Time substitution changed LHS type to $(typeof(lhs_term)).")
            rhs = substitute(rhs, sp)
            return lhs_term, rhs
        end
    end
    return lhs, rhs
end

function diffQEq(qspace::QSpace, left_hand_side::QAtomProduct, expr::QExpr)
    left_hand_side, expr = _normalize_diff_time(qspace, left_hand_side, expr)
    left_hand_side = Expectation(left_hand_side)
    expr = Expectation(expr)
    @assert !(contains_non_simple_QObj(expr)) "Differential requires simple AbstractQSums, i.e. no aggregators in QComposites (such as QExp, QLog...) and no nested AbstractQSums (multiple and complex indexing at the same level is possible, and immediate nesting is automatically simplified to composite indexes)."
    if !contains_abstract(left_hand_side) && !contains_abstract(expr)
        return reorder(neq(diffQEq(qspace, left_hand_side, expr, Val(:nosimp))))
    else
        return diffQEq(qspace, left_hand_side, expr, Val(:nosimp))
    end
end
function diffQEq(qspace::QSpace, left_hand_side::QAtomProduct, expr::QExpr, ::Val{:nosimp})
    left_hand_side = Expectation(left_hand_side)
    expr = Expectation(expr)
    return diffQEq(qspace, left_hand_side, expr, Val(:raw))
end


"""
    d_dt(qspace::QSpace, expr)

Evaluate the time derivative of an expression `expr` in the context of the given state space `ss`.

This function expects that `expr` is an equation (i.e. an Expr with an equal sign as its head),
of the form

    LHS = RHS
The function then returns a `diffQEq` constructed from the left-hand side QTerm and the right-hand side QExpr.
"""
function d_dt(left_hand::Union{QAtomProduct,QExpr}, right_hand::QExpr)::diffQEq
    # Check if expr is an equality.
    qstate = right_hand.qspace

    if left_hand isa QExpr
        if left_hand.qspace != qstate
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
    # Return a diffQEq constructed from these sides.
    return diffQEq(qstate, left_hand, right_hand)
end
function Expectation(eq::diffQEq)
    lhs = Expectation(eq.left_hand_side)
    rhs = Expectation(eq.expr)
    return diffQEq(eq.qspace, lhs, rhs, Val(:raw))
end

#### Helper Functions #######################################################################################
QNotAtom = Union{QComposite, QExpr, diffQEq}
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


include("QExpressionsOps/QExpressions_base_operators.jl")
include("QExpressionsOps/QExpressions_sort.jl")
include("QExpressionsOps/QExpressions_simplify.jl")

include("QExpressionsOps/QExpressions_properties.jl")
include("QExpressionsOps/QExpressions_algebra.jl")
include("QExpressionsOps/QExpressions_print.jl")

include("QExpressionsOps/QSum_modify.jl")

include("QExpressionsOps/QExpressions_welldefined.jl")
include("QExpressionsOps/QExpressions_substitute.jl")
include("QExpressionsOps/QExpressions_reorder.jl")
include("QExpressionsOps/QSum_decollision.jl") 

include("QExpressionsOps/OrderedQExpressions.jl")

import ..CFunctions: define_cabstract, define_ctype, list_cabstracts, list_ctypes, c_abstract_exists

function QExpr2CFunction(q::QExpr)::CFunction 
    if isnumeric(q) 
        if length(q.terms) < 1 
            error("Cannot extract CFunction from empty QExpr.") 
        end 
        coeff = get_coeff(q.terms[1]) 
        for r in q.terms[2:end] 
            coeff += get_coeff(r) 
        end 
    return coeff 
    else 
        error("Requires numeric QExpr, no quantum operators present.") 
    end 
end 

"""
    @define qspace, name
    @define qspace, name, fun

Two forms:

1) `@define ss name`
   - Calls `define_cabstract(ss, name)` to register a new abstract and returns a `QExpr`.
   - Binds a global `const name::QExpr` in the caller's module.

2) `@define ss Name fun`
   - Registers a custom ctype named `String(name)` using `fun` (a *numeric* `QExpr`)
     converted via `QExpr2CFunction`.
   - Defines constructors:
       Name(args::CFunction...)
       Name(coeff::ComplexRational, args::CFunction...)
       Name(qs::QExpr...)
       Name(coeff::ComplexRational, qs::QExpr...)
"""
macro define(qspace, name, fun=nothing)
    # Normalize the binding name to a Symbol (for variables & method names)
    n_sym = name isa Symbol ? name : Symbol(name)
    n_str = String(n_sym)   # for registration APIs
    ctype_sym = gensym(:ctype)  # internal const for ctype, not user-facing
    CR1 = ComplexRational(1,0,1)

    if fun === nothing
        # 1) @define ss name
        return esc(quote
            # Register and build a QExpr for the new abstract
            const $(n_sym)::QExpr = begin
                define_cabstract($qspace.param_info, $n_str)
                # Build the QExpr representing this abstract
                let __ab__ = $qspace.param_info.abstract_definitions[end]
                    QExpr($qspace, [
                        QAtomProduct($qspace,
                                     CAbstract($qspace.param_info,
                                               $CR1,
                                               __ab__.index))
                    ])
                end
            end
            $(n_sym)  # make the macro call evaluate to the QExpr
        end)
    else
        # 2) @define ss Name fun
        # fun is provided by the caller; don’t eval it in the macro. Convert at runtime.
        return esc(quote
            const $(ctype_sym) = define_ctype($qspace.param_info, $n_str, QExpr2CFunction($fun))
            if $(ctype_sym).has_abstract
                # Constructor 3: Name(qs::QExpr...)
                function $(n_sym)(qs::QExpr...)
                    c_exprs::Vector{CFunction} = QExpr2CFunction.(collect(qs))
                    return QExpr($qspace, QAtomProduct($qspace, CCustomType($qspace.param_info, $CR1, c_exprs, $(ctype_sym)), QTerm[]))
                end
            else
                $(n_sym)::QExpr = QExpr($qspace, QAtomProduct($qspace, CCustomType($qspace.param_info, $CR1, CFunction[$(ctype_sym).fun], $(ctype_sym)), QTerm[]))
            end
            $(n_sym)
        end)
    end
end
macro define_basics(qspace, var=:q)
    return esc(quote
        var0 = @define($qspace, var)
        @define($qspace, sin, 1//(2*1im) * (exp(1im*var0) - exp(-1im*var0)))
        @define($qspace, cos, 1//2 * (exp(1im*var0) + exp(-1im*var0)))
    end)
end

list_cabstracts(qspace::QSpace) = list_cabstracts(qspace.param_info)
list_ctypes(qspace::QSpace) = list_ctypes(qspace.param_info)
end
