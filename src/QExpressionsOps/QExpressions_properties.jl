export is_numeric, contains_abstract, contains_time

""" 
    is_numeric(t::QTerm, qspace::StateSpace) -> Bool
    is_numeric(t::QAbstract, qspace::StateSpace) -> Bool
    is_numeric(p::QAtomProduct) -> Bool
    is_numeric(s::QSum) -> Bool
    is_numeric(expr::QExpr) -> Bool

Returns true either if it is zero or it has only neutral elements for operators.
"""
function is_numeric(op_indices::Vector{Vector{Int}}, statespace::StateSpace)::Bool
    return statespace.I_op == op_indices
end
function is_numeric(t::QTerm, statespace::StateSpace)::Bool
    return statespace.I_op == t.op_indices
end
function is_numeric(t::QAbstract, statespace::StateSpace)::Bool
    return false
end
function is_numeric(e::QAtomProduct)
    if length(e.expr) > 1
        return false 
    elseif length(e.expr) == 0 
        return true
    else 
        return iszero(e.coeff_fun) || all([is_numeric(e.expr[1], e.statespace) for x in e.expr]) 
    end
end
function is_numeric(e::T) where T<:QComposite
    return iszero(e.coeff_fun)
end
function is_numeric(e::T) where T<:QMultiComposite
    return iszero(e.coeff_fun) 
end

function is_numeric(s::QSum)::Bool
    return false # is_numeric(s.expr)
end
function is_numeric(expr::QExpr)::Bool
    terms = expr.terms
    isempty(terms) && return true  # No terms = numeric 0
    return all(is_numeric, terms)
end


""" 
    contains_abstract(q::QObj) -> Bool

Checks if the quantum object contains an abstract operator among its leaves.
"""
function contains_abstract(term::QExpr)::Bool
    return any(contains_abstract, term.terms)
end

function contains_abstract(term::T)::Bool where T<:QComposite
    return contains_abstract(term.expr)
end

function contains_abstract(term::T)::Bool where T<:QMultiComposite
    return any(contains_abstract, term.expr)
end

function contains_abstract(term::QAtomProduct)::Bool
    return any(t -> isa(t, QAbstract), term.expr)
end

function contains_abstract(term::diff_QEq)::Bool 
    return contains_abstract(term.expr) && contains_abstract(term.left_hand_side)
end

import ..CFunctions: contains_c_indexes
""" 
    contains_c_indexes(f::Union{CFunction, QObj}, indexes::Vector{Int})::Bool

Checks if any of the CFunctions (in a QObj) depend on the indexes. 
"""
contains_c_indexes(q::QExpr, indexes::Vector{Int})::Bool = any(q -> contains_c_indexes(q, indexes), q.terms) 
contains_c_indexes(q::QAtom, indexes::Vector{Int}) = error("Cannot be applied to QAtom")
function contains_c_indexes(q::T, indexes::Vector{Int})::Bool where T <: QComposite 
    return contains_c_indexes(q.coeff_fun, indexes) || contains_c_indexes(q.expr, indexes)
end
function contains_c_indexes(q::M, indexes::Vector{Int})::Bool where M <: QMultiComposite
    return contains_c_indexes(q.coeff_fun) || any(t -> contains_c_indexes(x, indexes), q.expr)
end
contains_c_indexes(q::QSum, indexes::Vector{Int})::Bool = any(q -> contains_c_indexes(q, indexes), q.expr) 
contains_c_indexes(q::diff_QEq, indexes::Vector{Int})::Bool = contains_c_indexes(q.expr, indexes)
function get_t_indexes(param_info::ParameterInfo, t_ind::Int=-1)::Vector{Int} 
    if t_ind == -1 
        return param_info.indexes_of_t
    else
        return param_info.indexes_by_t_index[t_ind+1]
    end
end
""" 
    contains_time(q::QObj) -> Bool 

Checks is the quantum object depends on time. Doesn't work for QAtoms!
"""
contains_time(q::T, t_ind=-1) where T<: QAtom = error("Cannot get time indexes from QAtom. Try QComposites, QExpr, of diff_QEq instead. ")
function contains_time(q::T; t_ind=-1)::Bool where T <: QObj
    indexes = get_t_indexes(q.param_info, t_ind)
    return contains_c_indexes(q, indexes)
end

###################

function simple_isa(q::QExpr, type::Type)::Bool
    return length(q) == 1 && isa(q.terms[1], type)
end
function isaQAtomProduct(q::QExpr)::Bool
    if length(q) > 1
        return false
    else
        # length(q) == 1 
        if isa(q.terms[1], QAtomProduct)
            return true
        end
    end
end
import Base: isone
function isone(q::QAtomProduct)::Bool
    if is_numeric(q) && isnumeric(q.coeff_fun)
        return isone(c)
    end
    return false
end
function isone(q::QExpr)::Bool
    return simple_isa(q, QAtomProduct) && isone(q.terms[1])
end

##################

function where_neutral(q::QTerm, statespace::StateSpace)::Vector{Bool}
    return [op == neut for (op, neut) in zip(q.op_indices, statespace.I_op)]
end
function where_neutral(q::QAbstract, statespace::StateSpace)::Vector{Bool}
    return q.operator_type.expanded_ss_acting   # should never be modified! copy would be safer, but slower
end
function where_acting(q::QTerm, statespace::StateSpace)::Vector{Bool}
    return [op != neut for (op, neut) in zip(q.op_indices, statespace.I_op)]
end
function where_acting(q::QAbstract, statespace::StateSpace)::Vector{Bool}
    return .!q.operator_type.expanded_ss_acting  # should never be modified! copy would be safer, but slower
end
function where_acting(q::QAtomProduct)
    # combine the action of all of its constituents via OR 
    if length(q.expr) == 0
        return zeros(Bool, length(statespace.I_op))
    else
        return reduce(.|, [where_acting(expr, statespace) for expr in q.expr])
    end
end

function commutes_QAtom(q1::QAbstract, q2::QAbstract, statespace::StateSpace)::Bool   # for QAtom can check 
    # check if all elements of where neutral are NAND
    return statespace.operatortype_info.commute_fun(q1.key_index, q1.sub_index, q1.dag, q2.key_index, q2.sub_index, q2.dag)
end
@inline commutes_QAtom_inds(inds::Vector{Int}, q1::QAbstract, q2::QAbstract, statespace) = commutes_QAtom(q1, q2, statespace) 

function commutes_QAtom(q1::QTerm, q2::QTerm, statespace::StateSpace)::Bool
    a_q1 = where_acting(q1, statespace)
    a_q2 = where_acting(q2, statespace)
    inds = findall(a_q1 .& a_q2)
    isempty(inds) && return true
    return commutes_QAtom_inds(inds, q1, q2, statespace)
end
@inline function commutes_QAtom_inds(inds::Vector{Int}, q1::QTerm, q2::QTerm, statespace::StateSpace)::Bool
    @inbounds for ind in inds
        if !statespace.subspaces[statespace.subspace_info.outer_ss_of_expanded[ind]].op_set.commutes(q1[ind], q2[ind])
            return false
        end
    end
    return true
end

# Add the mixed method once:
function commutes_QAtom(qt::QTerm, qa::QAbstract, statespace::StateSpace)::Bool
    a_t = where_acting(qt, statespace)
    a_a = where_acting(qa, statespace)
    return !any(a_t .& a_a) 
end
@inline commutes_QAtom(qa::QAbstract, qt::QTerm, statespace::StateSpace) = commutes_QAtom(qt, qa, statespace::StateSpace)
@inline commutes_QAtom_inds(inds::Vector{Int}, q1::QTerm, q2::QAbstract, statespace::StateSpace) = length(inds) == 0
@inline commutes_QAtom_inds(inds::Vector{Int}, q1::QAbstract, q2::QTerm, statespace::StateSpace) = length(inds) == 0

function any_overlaps(multi_where_acting::Vector{Vector{Bool}})
    n = length(multi_where_acting)
    if n ≤ 1
        return false, multi_where_acting[1]
    end

    added = multi_where_acting[1]
    for i in 2:n
        if any(added .& multi_where_acting[i])
            return true, added
        end
        added = added .| multi_where_acting[i]
    end
    return false, added
end
@inline function commutes(q1::QAtomProduct, q2::QAtomProduct)::Bool
    statespace = q1.statespace
    acts1 = where_acting.(q1.expr, Ref(statespace))  # cache acting masks for q1 atoms
    acts2 = where_acting.(q2.expr, Ref(statespace))  # cache acting masks for q2 atoms
    @inbounds for (ai, where_a1) in zip(q1.expr, acts1)
        for (aj, where_a2) in zip(q2.expr, acts2)
            inds = findall(where_a1 .& where_a2)          # overlap indices for (ai, bj)
            if !commutes_QAtom_inds(inds, ai, aj, statespace)
                return false
            end
        end
    end
    return true
end
function commutes(Q1::QExpr, Q2::QExpr)::Bool
    statespace = Q1.statespace
    # collect non-commuting pairs
    noncomm_pairs = Tuple{Int,Int}[]
    for (i, x1) in enumerate(Q1.terms)
        for (j, x2) in enumerate(Q2.terms)
            if !commutes(x1, x2)#, statespace)
                push!(noncomm_pairs, (i, j))
            end
        end
    end

    isempty(noncomm_pairs) && return true      # all commute
    length(noncomm_pairs) == 1 && return false # exactly one conflict → cannot cancel

    # build commutator sum for all non-commuting pairs
    comms = QExpr([])
    for (i, j) in noncomm_pairs
        push!(comms.terms, Commutator(Q1.terms[i], Q2.terms[j]))
    end

    simplified = simplify_QExpr(comms)
    return isempty(simplified.terms) || all(iszero, simplified.terms)
end
function commutes(Q1::S, Q2::T) where {S<:QComposite,T<:QComposite}
    return commutes(Q1.expr, Q2.expr)
end

# define internal commutes function for QMultiComposite 
# do the internal degrees of freedom commute? 
function QCommutator_commutes(Q::QCommutator)::Bool
    return commutes(Q.expr[1], Q.expr[2])
end
each_commutes(exprs::Vector{QExpr}, Q2::QExpr)::Bool = all(commutes(Q1, Q2) for Q1 in exprs)
each_commutes(Q1::QExpr, exprs::Vector{QExpr})::Bool = all(commutes(Q1, Q2) for Q2 in exprs)
each_commutes(exprs1::Vector{QExpr}, exprs2::Vector{QExpr}) = all(commutes(Q1, Q2) for Q1 in exprs1, Q2 in exprs2)

commutes(Q1::S, Q2::T) where {S<:QMultiComposite,T<:QMultiComposite} = each_commutes(Q1.expr, Q2.expr)
commutes(Q1::S, Q2::T) where {S<:QComposite,T<:QMultiComposite} = each_commutes(Q1.expr, Q2.expr)
commutes(Q2::S, Q1::T) where {S<:QMultiComposite,T<:QComposite} = each_commutes(Q1.expr, Q2.expr)


###########################################################
import Base: ==

function ==(a::QTerm, b::QTerm)
    return a.op_indices == b.op_indices
end
function ==(a::QAbstract, b::QAbstract)
    return a.key_index == b.key_index && a.sub_index == b.sub_index && a.exponent == b.exponent && a.dag == b.dag && a.index_map == b.index_map
end
function ==(a::QAtomProduct, b::QAtomProduct)
    return a.coeff_fun == b.coeff_fun && all([ai == bi for (ai, bi) in zip(a.expr, b.expr)])
end

function ==(a::QExpr, b::QExpr)
    if length(a) != length(b)
        return false
    end
    if a.statespace != b.statespace
        return false
    end
    return all([ai == bi for (ai, bi) in zip(a, b)])
end
function ==(a::QSum, b::QSum)
    if a.element_indexes != b.element_indexes
        return false
    end
    if a.subsystem_index != b.subsystem_index
        return false
    end
    if a.neq != b.neq
        return false
    end
    return a.expr == b.expr
end
function ==(expr::QExpr, n::Number)
    if is_numeric(expr)
        if length(simple_expr.terms) == 0
            return isapprox(0, n)
        else
            return isapprox(simple_expr.terms[1].coeff, n)
        end
    end
    return false
end
function ==(n::Number, expr::QExpr)
    return expr == n  # Symmetric
end
