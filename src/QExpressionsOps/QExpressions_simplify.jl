function same_term_type(t1::QAtom, t2::QAtom) # generic case 
    return false
end
function same_term_type(t1::QTerm, t2::QTerm)
    return t1.op_indices == t2.op_indices && t1.time_index == t2.time_index
end
function same_term_type(t1::QAbstract, t2::QAbstract) 
    return t1.key_index == t2.key_index && t1.sub_index == t2.sub_index && t1.index_map == t2.index_map && t1.time_index == t2.time_index
end
function same_term_type(t1::QAtomProduct, t2::QAtomProduct)::Bool  # different coefficients but same operator content
    # check operators individually -> must all be the same!
    if length(t1.expr) != length(t2.expr) || t1.separate_expectation_values != t2.separate_expectation_values 
        return false
    end
    for (a,b) in zip(t1.expr, t2.expr) 
        if !same_term_type(a,b)
            return false 
        end
    end
    return true
end

function same_term_type(s1::QSum, s2::QSum)::Bool
    return s1.eq_indexes == s2.eq_indexes && s1.neq_blocks == s2.neq_blocks
end
function same_term_type(e1::QExpr, e2::QExpr)::Bool
    return all(same_term_type.(e1.expr, e2.expr))
end
function same_term_type(q1::S, q2::S)::Bool   where S<:QMultiComposite
    return all(same_term_type.(q1.expr, q2.expr))
end
function same_term_type(q1::S, q2::S)::Bool   where S<:QComposite
    return same_term_type(q1.expr, q2.expr)
end
function same_term_type(q1::S, q2::S)::Bool   where S<:QCompositeN
    return same_term_type(q1.expr, q2.expr) && q1.n == q2.n
end
function same_term_type(t1::S, t2::T)::Bool where {S<:QComposite, T<:QComposite}   # term grouping either not implemented for this object or objects aren'T the same 
    return false 
end

function combine_term_sum(t1::QAtomProduct, t2::QAtomProduct)::QAtomProduct
    return QAtomProduct(t1.qspace, t1.coeff_fun + t2.coeff_fun, t1.expr)
end
function combine_term_sum(s1::QSum, s2::QSum)::QSum
    return only(modify_expr(s1, s1.expr + s2.expr))
end
function combine_term_sum(q1::S, q2::S)::QComposite where S<:QComposite 
    return modify_coeff(q1, q1.coeff_fun + q2.coeff_fun)
end

function simplify_QExpr(terms::Vector{QComposite})::Vector{QComposite}
    # If there are no terms, return an empty QExpr.
    if isempty(terms)
        return QComposite[]
    end
    
    # First, sort QExpr without modifying the original.
    sort!(terms)

    combined_terms = QComposite[]
    i = 1
    curr_term = terms[1]
    while i < length(terms)
        # Combine adjacent like terms.
        next_term = terms[i+1]
        if same_term_type(curr_term, next_term)
            curr_term = combine_term_sum(curr_term, next_term)
        else
            if !iszero(curr_term)
                push!(combined_terms, curr_term)
            end
            curr_term = next_term
        end
        i += 1
    end
    if !iszero(curr_term)
        push!(combined_terms, curr_term)
    end
    return combined_terms
end
function simplify_QExpr(q::QExpr)::QExpr 
    return QExpr(q.qspace, simplify_QExpr(q.terms))
end

#####################################################################################################################

include("QExpressions_Simplify/QExpressions_simplify_atoms.jl")
include("QExpressions_Simplify/QExpressions_simplify_composites.jl")

#####################################################################################################################
function simplify_QExp(q::QExp)::Vector{QComposite}
    if length(q.expr) == 1 && isa(q.expr[1], QLog)
        return q.expr[1].expr
    end
    if isnumeric(q.expr) 
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in q.expr)
        return [ modify_coeff(q.expr[1], q.coeff_fun * exp(sum_of_coeff_funs)) ]
    end
    return [ q ]
end

function simplify_QLog(q::QLog)::Vector{QComposite}
    if length(q.expr) == 1 && isa(q.expr[1], QExp)
        return q.expr[1].expr
    end
    if isnumeric(q.expr)
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in q.expr)
        return [ modify_coeff(q.expr[1], q.coeff_fun * log(sum_of_coeff_funs)) ]
    end
    return [ q ]
end

function simplify_QPower(q::QPower)::Vector{QComposite}
    if isnumeric(q.expr)
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in q.expr)
        return [ modify_coeff(q.expr[1], q.coeff_fun * power(sum_of_coeff_funs, q.n)) ]
    end
    return [ q ]
end
function simplify_QRoot(q::QRoot)::Vector{QComposite}
    if isnumeric(q.expr)
        sum_of_coeff_funs = sum(qi.coeff_fun for qi in q.expr)
        return [ modify_coeff(q.expr[1], q.coeff_fun * root(sum_of_coeff_funs, q.n)) ]
    end    
    return [ q ]
end


"""
    simplify(q) -> simplified

Simplifies `QExpr`-based symbolic quantum expressions by recursively reducing internal structures:

- `QAtomProduct`: Applies pairwise simplifications repeatedly and merges results into canonical `QAtomProduct`s.
- `QMultiComposite`: Simplifies each expression element-wise.
- `QComposite`: Simplifies its internal expression and returns a new `QComposite`.
- `QExpr`: Flattens and simplifies terms, then combines like terms where possible.
- `QSum`: Simplifies its expression array and returns a new `QSum`.
- `diff_QEq`: Replaces its right-hand side with a simplified version.

Returns either a single simplified object or a list of canonical components depending on input type.

function simplify(q::QAtomProduct)::Vector{QComposite}
    return simplify_QAtomProduct(q)
end

function simplify(q::T)::Vector{QComposite} where T <: QMultiComposite
    simplified = simplify.(q.expr)
    return modify_expr(q, map(only, simplified))
end

function simplify(q::T)::Vector{QComposite} where T <: QComposite
    return modify_expr(q, simplify(q.expr))
end 

function simplify(q::QLog)::Vector{QComposite} 
    return simplify_QLog(only(modify_expr(q, simplify(q.expr))))
end 
function simplify(q::QExp)::Vector{QComposite} 
    return simplify_QExp(only(modify_expr(q, simplify(q.expr))))
end 

function simplify(q::QExpr)::QExpr
    return QExpr(q.qspace, simplify) #simplify_QExpr(mapreduce(simplify, vcat, q.terms))) # don't also simplify its elements. they should already be simplified! 
end
function simplify(q::diff_QEq)::diff_QEq
    simp_rhs = simplify(q.expr)
    return diff_QEq(q.qspace, q.left_hand_side, simp_rhs, q.do_braket)
end
"""
