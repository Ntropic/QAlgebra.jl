export is_t_var, is_t, is_local, contains_non_simple_QObj, contains_non_simple, contains_abstract, contains_time, contains_which_t_indexes, max_moment_of_terms, where_acting
export is_unitary, is_hermitian, substitution_properties_fulfilled, same_qspace, qspace_check
import ..CFunctions: isnumeric, CFunction, CAtom
""" 
    isnumeric(t::QObj) -> Bool

Returns true either if it is zero or it has only neutral elements for operators.
"""
function isnumeric(op_indices::Vector{Vector{Int}}, qspace::QSpace)::Bool
    return qspace.I_op == op_indices
end
function isnumeric(t::QTerm, qspace::QSpace)::Bool
    return qspace.I_op == t.op_indices
end
function isnumeric(t::QTerm, index::Int, qspace::QSpace)::Bool
    return qspace.I_op[index] == t.op_indices[index]
end
function isnumeric(t::QAbstract, qspace::QSpace)::Bool
    return false
end
function isnumeric(e::QAtomProduct)
    if length(e.expr) > 1
        return false 
    elseif length(e.expr) == 0 
        return true
    else 
        return iszero(e.coeff_fun) || all([isnumeric(e.expr[1], e.qspace) for x in e.expr]) 
    end
end
function isnumeric(e::T) where T<:QComposite
    return iszero(e.coeff_fun)
end
function isnumeric(e::T) where T<:QMultiComposite
    return iszero(e.coeff_fun) 
end

function isnumeric(s::QSum)::Bool
    return false # isnumeric(s.expr)
end
function isnumeric(expr::QExpr)::Bool
    terms = expr.terms
    isempty(terms) && return true  # No terms = numeric 0
    return all(isnumeric, terms)
end

function is_t_var(t::QExpr)::Bool
    # must be a single QAtomProduct term that is numeric
    if length(t) != 1 || !(t.terms[1] isa QAtomProduct) || !isnumeric(t.terms[1])
        return false
    end

    coeff_fun = t.terms[1].coeff_fun
    # coeff_fun must be a CAtom with unit coefficient
    if !isa(coeff_fun, CAtom) || !isone(coeff_fun.coeff)
        return false
    end

    inds = findall(!=(0), coeff_fun.var_exponents)
    # must be exactly one variable with exponent 1
    return length(inds) == 1 &&
           coeff_fun.var_exponents[inds[1]] == 1 &&
           t.qspace.params[inds[1]].is_t
end

is_t(t::QExpr)::Bool = is_t_var(t)

function is_t(prod::QAtomProduct)::Bool
    isnumeric(prod) || return false
    return _is_t_cfunction(prod.coeff_fun, prod.qspace)
end

is_t(::QObj)::Bool = false

@inline function _is_t_cfunction(f::CAtom, qspace::QSpace)::Bool
    isone(f.coeff) || return false
    exps = f.var_exponents
    idxs = findall(!iszero, exps)
    length(idxs) == 1 || return false
    idx = idxs[1]
    exps[idx] == 1 || return false
    return qspace.params[idx].is_t
end

@inline _is_t_cfunction(::CFunction, ::QSpace)::Bool = false

"""
    contains_non_simple_QObj(q::QObj) -> Bool 

Does a QObj contain non Basic Quantum Objects, such as QExp, QLog or QMultiComposites?
"""
contains_non_simple_QObj(q::QExpr, sum_is_simple::Bool=true)::Bool = any(t -> contains_non_simple_QObj(t, sum_is_simple), q.terms)
contains_non_simple_QObj(q::T, sum_is_simple::Bool=true) where {T <: QComposite} = true 
contains_non_simple_QObj(q::QAtomProduct, sum_is_simple::Bool=true)::Bool = false 
function contains_non_simple_QObj(q::QSum, sum_is_simple::Bool=true)::Bool 
    if sum_is_simple
        return any(x -> contains_non_simple_QObj(x, sum_is_simple), q.expr)
    end 
    return true
end

""" 
    contains_non_simple(q::QObj) -> Bool

Checks if it contains any non-simple QObjects, QAbstracts or non-simple CFunctions.
"""
contains_non_simple(q::QExpr, sum_is_simple::Bool=true)::Bool = any(t -> contains_non_simple(t, sum_is_simple), q.terms)
contains_non_simple(q::T, sum_is_simple::Bool=true) where {T <: QComposite} = true 
contains_non_simple(q::QAtomProduct, sum_is_simple::Bool=true)::Bool = contains_non_simple_CFunction(q.coeff_fun) || contains_abstract(q)
function contains_non_simple(q::QSum, sum_is_simple::Bool=true)::Bool
    if sum_is_simple
        return any(t -> contains_non_simple(t, sum_is_simple), q.expr)
    end
    return true 
end

""" 
    is_local(q::QObj) -> Bool 

Does an operator contains only local operations, or also entangling operators?
""" 
is_local(q::QExpr)::Bool = all(is_local, q.terms)
is_local(q::T) where {T <: QComposite} = is_local(q.expr)
is_local(q::QCompositeProduct) = false
is_local(q::QCommutator)::Bool = is_local(q.expr[1]*q.expr[2] - q.expr[2]*q.expr[1])
is_local(q::QAtomProduct)::Bool = sum(where_acting(q))<= 1

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
function contains_c_indexes(q::QAtomProduct, indexes::Vector{Int})::Bool 
    return contains_c_indexes(q.coeff_fun, indexes) 
end
function contains_c_indexes(q::T, indexes::Vector{Int})::Bool where T <: QComposite 
    return contains_c_indexes(q.coeff_fun, indexes) || contains_c_indexes(q.expr, indexes)
end
function contains_c_indexes(q::M, indexes::Vector{Int})::Bool where M <: QMultiComposite
    return contains_c_indexes(q.coeff_fun) || any(t -> contains_c_indexes(x, indexes), q.expr)
end
contains_c_indexes(q::QSum, indexes::Vector{Int})::Bool = any(q -> contains_c_indexes(q, indexes), q.expr) 
contains_c_indexes(q::diff_QEq, indexes::Vector{Int})::Bool = contains_c_indexes(q.expr, indexes)


contains_t_indexes(q::QExpr, indexes::Vector{Int})::Bool = any(q -> contains_t_indexes(q, indexes), q.terms) 
contains_t_indexes(q::QAtom, indexes::Vector{Int}) = error("Cannot be applied to QAtom")
function contains_t_indexes(q::QAtomProduct, indexes::Vector{Int})::Bool 
    return contains_c_indexes(q.coeff_fun, indexes) || any(x -> x.time_index != -1, q.expr)
end
function contains_t_indexes(q::T, indexes::Vector{Int})::Bool where T <: QComposite 
    return contains_c_indexes(q.coeff_fun, indexes) || contains_t_indexes(q.expr, indexes)
end
function contains_t_indexes(q::M, indexes::Vector{Int})::Bool where M <: QMultiComposite
    return contains_c_indexes(q.coeff_fun) || any(t -> contains_t_indexes(x, indexes), q.expr)
end
contains_t_indexes(q::QSum, indexes::Vector{Int})::Bool = any(q -> contains_t_indexes(q, indexes), q.expr) 
contains_t_indexes(q::diff_QEq, indexes::Vector{Int})::Bool = contains_t_indexes(q.expr, indexes) || contains_t_indexes(q.left_hand_side)
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
@inline contains_time(q::T, t_ind=0) where T<:QAtom = error("Cannot get time indexes from QAtom. Try QComposites, QExpr, of diff_QEq instead. ")

@inline function contains_time(q::T, t_ind=0)::Bool where T <: QObj
    indexes = get_t_indexes(q.qspace.param_info, t_ind)
    return contains_t_indexes(q, indexes)
end
""" 
    contains_which_t_indexes(q::QObj) -> BitVector 

Returns a Boolean Vector of whether each time index is present in the QObj, 
with time indexes starting at `t_index=0` and ending at `t_index=max_t_ind`  
"""
contains_which_t_indexes(q::T) where T<:QAtom = error("Cannot get time indexes from QAtom. Try QComposites, QExpr, of diff_QEq instead. ")
function contains_which_t_indexes(q::T)::BitVector where T <: QObj
    max_t_index = q.qspace.max_t_ind
    return [contains_t_indexes(q, get_t_indexes(q.qspace.param_info, t_ind)) for t_ind in 0:max_t_index] 
end

@inline function max_moment_of_terms(q::QAtomProduct)::Int
    @assert length(q.expr) == 1 && isa(q.expr[1], QTerm) "QAtomProduct can only contain a single QTerm to specify moment of operator."
    return sum(where_acting(q.expr[1], q.qspace))
end 
@inline function max_moment_of_terms(q::T)::Int where T <: QComposite
    return max_moment_of_terms(q.expr)
end
@inline function max_moment_of_terms(q::T)::Int where T <: QMultiComposite
    return maximum(max_moment_of_terms.(x) for x in q.expr)
end
@inline function max_moment_of_terms(q::QExpr)::Int 
    return maximum(max_moment_of_terms(t) for t in q.terms)
end
@inline function max_moment_of_terms(q::diff_QEq)::Int 
    return max(max_moment_of_terms(q.left_hand_side), max_moment_of_terms(q.expr))
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
import Base: isone, iszero
function isone(q::QAtomProduct)::Bool
    if isnumeric(q) && isnumeric(q.coeff_fun)
        return isone(q.coeff_fun)
    end
    return false
end
function isone(q::QExpr)::Bool
    return simple_isa(q, QAtomProduct) && isone(q.terms[1])
end

# Optionally, define length and eltype.
iszero(q::QExpr) = length(q.terms) == 0 || all(iszero, q.terms)
iszero(q::QAtomProduct) = iszero(q.coeff_fun)
iszero(q::QSum) = iszero(q.expr)
iszero(q::T) where T<:QComposite = iszero(q.coeff_fun) || iszero(q.expr)
iszero(q::T) where T<:QMultiComposite = iszero(q.coeff_fun) || any(iszero, q.expr) 

##################

function where_neutral(q::QTerm, qspace::QSpace)::BitVector
    return [op == neut for (op, neut) in zip(q.op_indices, qspace.I_op)]
end
function where_neutral(q::QAbstract, qspace::QSpace)::BitVector
    return q.operator_type.expanded_ss_acting   # should never be modified! copy would be safer, but slower
end
function where_acting(q::QTerm, qspace::QSpace)::BitVector
    return [op != neut for (op, neut) in zip(q.op_indices, qspace.I_op)]
end
function where_acting(q::QAbstract, qspace::QSpace)::BitVector
    return .!q.operator_type.expanded_ss_acting  # should never be modified! copy would be safer, but slower
end
function where_acting(q::QAtomProduct)::BitVector
    # combine the action of all of its constituents via OR
    qspace = q.qspace 
    if length(q.expr) == 0
        return falses(length(qspace.I_op))
    else
        return mapreduce(expr -> where_acting(expr, qspace), .|, q.expr)
    end
end
where_acting(q::QExpr)::BitVector = mapreduce(t -> where_acting(t), .|, q.terms)
function where_acting(q::T)::BitVector where {T<:QComposite}
    return where_acting(q.expr)
end
function where_acting(q::QCumulant)::BitVector
    acting = falses(length(q.qspace.I_op))
    for idx in q.where_acting
        acting[idx] = true
    end
    return acting .| where_acting(q.expr)
end
function where_acting(q::T)::BitVector where {T<:QMultiComposite}
    return mapreduce(expr -> where_acting(expr, qspace), .|, q.expr)
end
function where_acting(q::QSum)::BitVector
    acting = where_acting(q.expr)
    for ind in iter_all_indexes(q)
        acting[expanded(ind)] = true
    end
    return acting
end

function commutes_QAtom(q1::QAbstract, q2::QAbstract, qspace::QSpace)::Bool   # for QAtom can check 
    # check if all elements of where neutral are NAND
    if q1.time_index != q2.time_index
        return false 
    end
    return qspace.operatortype_info.commute_fun(q1.key_index, q1.sub_index, q1.dag, q2.key_index, q2.sub_index, q2.dag)
end
@inline commutes_QAtom_inds(inds::Vector{Int}, q1::QAbstract, q2::QAbstract, qspace) = commutes_QAtom(q1, q2, qspace) 

function commutes_QAtom(q1::QTerm, q2::QTerm, qspace::QSpace)::Bool
    if q1.time_index != q2.time_index
        return false 
    end
    a_q1 = where_acting(q1, qspace)
    a_q2 = where_acting(q2, qspace)
    inds = findall(a_q1 .& a_q2)
    isempty(inds) && return true
    return commutes_QAtom_inds(inds, q1, q2, qspace)
end
@inline function commutes_QAtom_inds(inds::Vector{Int}, q1::QTerm, q2::QTerm, qspace::QSpace)::Bool
    @inbounds for ind in inds
        if !qspace.subspaces[qspace.subspace_info.outer_ss_of_expanded[ind]].op_set.commutes(q1[ind], q2[ind])
            return false
        end
    end
    return true
end

# Add the mixed method once:
function commutes_QAtom(qt::QTerm, qa::QAbstract, qspace::QSpace)::Bool
    if qt.time_index != qa.time_index
        return false 
    end
    a_t = where_acting(qt, qspace)
    a_a = where_acting(qa, qspace)
    return !any(a_t .& a_a) 
end
@inline commutes_QAtom(qa::QAbstract, qt::QTerm, qspace::QSpace) = commutes_QAtom(qt, qa, qspace::QSpace)
@inline commutes_QAtom_inds(inds::Vector{Int}, q1::QTerm, q2::QAbstract, qspace::QSpace) = length(inds) == 0
@inline commutes_QAtom_inds(inds::Vector{Int}, q1::QAbstract, q2::QTerm, qspace::QSpace) = length(inds) == 0

function any_overlaps(multi_where_acting::Vector{BitVector})
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
    qspace = q1.qspace
    acts1 = where_acting.(q1.expr, Ref(qspace))  # cache acting masks for q1 atoms
    acts2 = where_acting.(q2.expr, Ref(qspace))  # cache acting masks for q2 atoms
    @inbounds for (ai, where_a1) in zip(q1.expr, acts1)
        for (aj, where_a2) in zip(q2.expr, acts2)
            inds = findall(where_a1 .& where_a2)          # overlap indices for (ai, bj)
            if !commutes_QAtom_inds(inds, ai, aj, qspace)
                return false
            end
        end
    end
    return true
end
function commutes(Q1::QExpr, Q2::QExpr)::Bool
    qspace = Q1.qspace
    # collect non-commuting pairs
    noncomm_pairs = Tuple{Int,Int}[]
    for (i, x1) in enumerate(Q1.terms)
        for (j, x2) in enumerate(Q2.terms)
            if !commutes(x1, x2)#, qspace)
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
function commutes(Q1::S, Q2::T) where {S<:QComposite,T<:QAtomProduct}
    return commutes(Q1.expr, QExpr(Q2.qspace, Q2))
end
function commutes(Q1::S, Q2::T) where {S<:QAtomProduct,T<:QComposite}
    return commutes(QExpr(Q1.qspace, Q1), Q2.expr)
end
# for QCompositeProduct we need to track this differently. 

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
    return a.op_indices == b.op_indices && a.time_index == b.time_index
end
function ==(a::QAbstract, b::QAbstract)
    return a.key_index == b.key_index && a.sub_index == b.sub_index && a.exponent == b.exponent && a.dag == b.dag && a.index_map == b.index_map && a.time_index == b.time_index
end
function ==(a::QAtomProduct, b::QAtomProduct)
    return a.coeff_fun == b.coeff_fun && all([ai == bi for (ai, bi) in zip(a.expr, b.expr)]) && a.qspace == b.qspace
end

function ==(a::QExpr, b::QExpr)
    if length(a) != length(b)
        return false
    end
    if a.qspace != b.qspace
        return false
    end
    return all([ai == bi for (ai, bi) in zip(a, b)])
end
function ==(a::QSum, b::QSum)
    a.qspace == b.qspace || return false
    a.eq_indexes == b.eq_indexes     || return false
    a.neq_blocks == b.neq_blocks || return false
    return a.expr == b.expr
end

function ==(expr::QExpr, n::Number)
    if isnumeric(expr)
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


#### Check substitution properties 
function is_unitary(q::QExpr)::Bool
    return isone(q*q')
end
function is_hermitian(q::QExpr)::Bool
    return q==q'
end
function substitution_properties_fulfilled(a::QAbstract, q::QExpr)::Bool 
    # check if all properties are fulfilled 
    if a.exponent != 1 || a.dag || length(a.index_map)>0
        error("Cannot substitute exponentiate, Daggered abstract operators or ones with index-maps defined. 
            Index maps are respected within substitution, but not in the substitution definition.")
    end
    op_type = a.operator_type
    of_time = op_type.of_time
    hermitian = op_type.hermitian 
    unitary = op_type.unitary
    expanded_ss_acting = op_type.expanded_ss_acting
    # check each of these:
    if !of_time && contains_time(q)
        error("Abstract operator is't time dependent but QExpr $q is. ")
    end
    acting = where_acting(q) 
    # check if any true element of acting is not true in expanded_ss_acting
    if any(acting .& .!expanded_ss_acting)
        error("Abstract operator is defined on expanded subspaces: $expanded_ss_acting, but QExpr acts on $acting.")
    end
    if hermitian && !is_hermitian(q) 
        error("Abstract operator expected to be hermitian, but QExpr is not: $q ≠ $(q').")
    end
    if unitary && !is_unitary(q)
        error("Abstract operator expected to be hermitian, but QExpr is not: $q ⋅ $(q') = $(q*q') instead of 1.")
    end
    return true 
end

########## Statespace check infra ##############################################
function same_qspace(a::S, b::T)::Bool where {S<:QNotAtom,T<:QNotAtom}
    return a.qspace === b.qspace
end
# Only check at outermost call sites. Internal calls use _NOCHK.
@inline qspace_check_if(::Val{true}, a, b) =
    (a.qspace === b.qspace) || _qspace_throw(a, b)
@inline qspace_check_if(::Val{false}, a, b) = nothing

@noinline function _qspace_throw(a, b)
    throw(AssertionError("Objects must share the same qspace; got $(summary(a)) vs $(summary(b))"))
end
