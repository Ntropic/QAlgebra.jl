export Dag, Commutator, is_numeric, same_statespace

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
    return all([is_numeric(x, e.statespace) for x in e.expr]) || iszero(e.coeff_fun)
end
function is_numeric(e::T) where T<:QComposite
    return is_numeric(e.expr) || iszero(e.coeff_fun)
    #error("is_numeric (without given statespace) not implemented for QComposite subtype $(typeof(e))")
end
function is_numeric(e::T) where T<:QMultiComposite
    return iszero(e.coeff_fun) || all([is_numeric(e1, statespace) for e1 in e.expr])
end

function is_numeric(s::QSum)::Bool
    return false # is_numeric(s.expr)
end
function is_numeric(expr::QExpr)::Bool
    terms = expr.terms
    if isempty(terms)
        return true  # No terms = numeric 0
    else
        return all(is_numeric(t) for t in terms)
    end
end

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


function where_neutral(q::QTerm, statespace::StateSpace)::Vector{Bool}
    return [op == neut for (op, neut) in zip(q.op_indices, statespace.I_op)]
end
function where_neutral(q::QAbstract, statespace::StateSpace)::Vector{Bool}
    return q.operator_type.subspaces   # should never be modified! copy would be safer, but slower
end
function where_acting(q::QTerm, statespace::StateSpace)::Vector{Bool}
    return [op != neut for (op, neut) in zip(q.op_indices, statespace.I_op)]
end
function where_acting(q::QAbstract, statespace::StateSpace)::Vector{Bool}
    return .!q.operator_type.subspaces  # should never be modified! copy would be safer, but slower
end
function where_acting(q::QAtomProduct)
    # combine the action of all of its constituents via OR 
    if length(q.expr) == 0
        return zeros(Bool, length(statespace.I_op))
    else
        return reduce(.|, [where_acting(expr, statespace) for expr in q.expr])
    end
end

function commutes_QAbstract(q1::QAbstract, q2::QAbstract, statespace::StateSpace)::Bool   # for QAtom can check 
    # check if all elements of where neutral are NAND
    return !any(.&(where_acting(q1, statespace), where_acting(q2, statespace)))
end
function commutes_QTerm(q1::QTerm, q2::QTerm, statespace::StateSpace)::Bool
    a_q1 = where_acting(q1, statespace)
    a_q2 = where_acting(q2, statespace)
    overlap = .&(a_q1, a_q2)
    inds = findall(a .& b)
    if length(inds) == 0 # commutes! 
        return true
    else # find overlap 
        for ind in inds
            if statespace.subspace[statespace.subspace_by_ind[ind]].commutes(q1[ind], q2[ind])
                return false
            end
        end
        return true
    end
end

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
function commutes(q1::QAtomProduct, q2::QAtomProduct)::Bool
    a_q1 = where_acting(q1, statespace)
    a_q2 = where_acting(q2, statespace)
    overlap = .&(a_q1, a_q2)
    inds = findall(a .& b)
    if length(inds) == 0 # commutes! 
        return true
    else
        term_q1 = isa.(q1.expr, Ref(QTerm))
        term_q2 = isa.(q2.expr, Ref(QTerm))
        # check if the QAtomProduct contains only QTerm objects or if not that the qAbstract terms overlap with none of the other terms in either q1 or q2 
        if all(term_q1) && all(term_q2)
            for ind in inds
                if statespace.subspace[statespace.subspace_by_ind[ind]].commutes(q1[ind], q2[ind])
                    return false
                end
            end
            return true
        else
            # Mixed QTerms and qAbstract → ensure no qAbstract in one overlaps with any term in the other
            # --- Final case: mixed QTerms and QAbstracts ---
            # Compare each QAtom in q1 with all in q2
            if sum(.!term_q1) > 1 || sum(.!term_q2) > 1
                return false
            end
            a_q1s = where_acting.(q1.expr[term_q1], Ref(statespace))
            a_q2s = where_acting.(q2.expr[term_q2], Ref(statespace))
            has_overlaps, added = any_overlaps(vcat(a_q1s, a_q2s))
            if has_overlaps
                return false
            else
                # check overlaps with qAtoms 
                t_acting = reduce(.&, vcat(q1.expr[.!term_q1], q2.expr[.!term_q2]))
                if any(.&(t_acting, added))
                    return false
                else
                    # check the qTerms amongst themselves. there should only be one qTerm in a QAtomProduct 
                    return commutes_QTerm(q1.expr[.!term_q1][1], q2.expr[.!term_q2][1])
                end
            end
        end
    end
end
function commutes(Q1::QExpr, Q2::QExpr)::Bool
    # check each term in Q1 with each term in Q2
    # for terms that don't commute with another, we multiply them with one another and check if they cancel out after all. 
    # element-wise test 
    elements_commute::Matrix{Bool} = zeros(Bool, length(Q1.expr), length(Q2.expr))
    for (i, t1) in enumerate(Q1.expr)
        for (j, t2) in enumerate(Q2.expr)
            elements_commute[i, j] = commutes(t1, t2)
        end
    end
    # find non commuters in matrix 
    non_commuter_inds = findall(.!elements_commute) # indices of non-commuting terms. 
    if length(non_commuter_inds) == 0
        return true # all commute, so the whole expression comm
    elseif length(non_commuter_inds) == 1
        return false # only one non-commuter, so the whole expression doesn't commute. 
    else
        # check if the non-commuting terms cancel out 
        qprods = QComposite[]
        for (i, j) in non_commuter_inds
            append!(qprods, Commutator(Q1.expr[i], Q2.expr[j]))
        end
        qprods_simplified = simplify_QExpr(qprods)
        if length(qprods_simplified) == 0
            return true # they cancel out, so the whole expression
        elseif length(qprods_simplified)
            return all(iszero(t) for t in qprods_simplified) # all zero
        end
    end
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

function same_statespace(a::S, b::T)::Bool where {S<:QComposite,T<:QComposite}
    return a.statespace == b.statespace
end

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


### Basic Operations: 

####### Unary Minus #############################################################
function -(t::QAtomProduct)::QAtomProduct
    return QAtomProduct(t.statespace, -t.coeff_fun, copy(t.expr))
end
function -(t::QExpr)::QExpr
    return QExpr(t.statespace, .-t.terms)
end
function -(t::T)::T where T<:QComposite
    return modify_expr(t, -t.expr)
end

#### Binary + ####################################################################
function +(Q1::QExpr, Q2::QExpr)::QExpr
    return QExpr(Q1.statespace, vcat(Q1.terms, Q2.terms))
end
function +(Q1::QExpr, Q2::T)::QExpr where {T<:QComposite}
    return QExpr(Q1.statespace, vcat(Q1.terms, Q2))
end
function +(Q1::QExpr, N::Number)::QExpr
    new_terms = vcat(Q1.terms, QAtomProduct(Q1.statespace, Q1.statespace.c_one * N, QTerm(Q1.statespace.I_op)))
    return QExpr(Q1.statespace, new_terms)
end
+(N::Number, Q1::QExpr)::QExpr = Q1 + N

#### Binary - ####################################################################
function -(Q1::QExpr, Q2::QExpr)::QExpr
    new_terms = simplify_QExpr(vcat(Q1.terms, .-Q2.terms))
    return QExpr(Q1.statespace, new_terms)
end
function -(Q1::QExpr, Q2::T)::QExpr where {T<:QComposite}
    new_terms = vcat(Q1.terms, -Q2)
    return QExpr(Q1.statespace, new_terms)
end
(-(Q2::T, Q1::QExpr)::QExpr) where {T<:QComposite} = -Q1 + Q2

function -(Q1::QExpr, N::Number)::QExpr
    new_terms = vcat(Q1.terms, QAtomProduct(Q1.statespace, -Q1.statespace.c_one * N, QTerm(Q1.statespace.I_op)))
    return QExpr(Q1.statespace, new_terms)
end
-(N::Number, Q1::QExpr)::QExpr = -Q1 + N

#### Multiply ####################################################################
function trivial_multiply(Q1::QAtomProduct, Q2::QAtomProduct)::QAtomProduct   
    e1, e2 = Q1.expr, Q2.expr
    n1, n2 = length(e1), length(e2)

    newexpr = Vector{QAtom}(undef, n1 + n2)
    if n1 > 0
        newexpr[1:n1] = e1
    end

    if n2 > 0
        newexpr[n1+1:end] = e2
    end
    return modify_coeff_expr(Q1, Q1.coeff_fun * Q2.coeff_fun, newexpr)
end
# Multiplies two QTerm’s from the same statespace. Returns a vector of QTerm’s that are the result of this multiplication and corresponding ComplexRational coefficients. 
function multiply_qterm(t1::Vector{Vector{Int}}, t2::Vector{Vector{Int}}, ss::StateSpace)::Tuple{Vector{Vector{Is}},Vector{ComplexRational}} 
    i = 0
    for s in ss.subspaces, _ in eachindex(s.ss_inner_ind)
        i += 1
        push!(results, s.op_set.op_product(t1[i], t2[i]))
    end

    new_inds = Vector{Is}[]
    new_coeffs = ComplexRational[]
    for combo in Base.Iterators.product(results...)
        push!(new_coeffs, foldl(*, (f for (f, _) in combo); init=ComplexRational(1, 0, 1)))
        push!(new_inds, [idx for (_, idx) in combo])
    end
    return new_inds, new_coeffs
end
function multiply_qterm(t1::QTerm, t2::QTerm, ss::StateSpace)::Tuple{Vector{QTerm},Vector{ComplexRational}}
    new_inds, new_coeffs = multiply_qterm(t1.op_indices, t2.op_indices, ss)
    return QTerm[QTerm(inds) for inds in new_inds], new_coeffs
end

function *(t1::QTerm, t2::QTerm, statespace::StateSpace)  # should never be called by the user!!
    return multiply_qterm(t1, t2, statespace)
end

#### Main Multiplication Functions ################################################
function *(p1::QAtomProduct, p2::QAtomProduct)::Vector{QComposite}
    p = trivial_multiply(p1, p2)  # append the terms of p1 and p2.
    return simplify_QAtomProduct(p)    # simplify the product. 
end

function *(num::Number, p1::QAtomProduct)::Vector{QComposite}
    if iszero(num) 
        return QComposite[] 
    else
        return [modify_coeff(p1, p1.coeff_fun * num)]
    end
end
function *(p1::QAtomProduct, num::Number)::Vector{QComposite}
    return num * p1
end

function *(p1::T1, p2::T2)::Vector{QComposite} where {T1<:QComposite,T2<:QComposite}
    return [QCompositeProduct(QComposite[p1, p2])]
end
function *(p1::QSum, p2::T2)::Vector{QSum} where T2<:QComposite
    new_expr::Vector{QComposite} = []
    for t in p1.expr
        append!(new_expr, t * p2)
    end
    return [QSum(p1.statespace, QExpr(new_expr, p1.indexes, p1.subsystem_index, p1.element_indexes, p1.neq), Val(:simp))]
end
function *(p1::QCompositeProduct, p2::T2)::Vector{QComposite} where T2<:QComposite
    curr_comp = p1.expr
    last_prod = p1.expr[end] * p2
    return_vec::Vector{QCompositeProduct} = []
    for last in last_prod
        if isa(last, QCompositeProduct)
            push!(return_vec, QCompositeProduct(vcat(curr_comp, QExpr([p2])), Val(:simp)))
        else
            push!(return_vec, QCompositeProduct(vcat(curr_comp[1:end-1], QExpr([last])), Val(:simp)))
        end
    end
    return return_vec
end
function *(p2::T2, p1::QCompositeProduct)::Vector{QComposite} where T2<:QComposite
    curr_comp = p1.expr
    first_prod = p2 * p1.expr[1]
    return_vec::Vector{QCompositeProduct} = []
    for first in first_prod
        if isa(first, QCompositeProduct)
            push!(return_vec, QCompositeProduct(vcat(p2, curr_comp)))
        else
            push!(return_vec, QCompositeProduct(vcat(first_prod, curr_comp[2:end])))
        end
    end
    return return_vec
end

function *(Q1::QExpr, Q2::QExpr)::QExpr
    if Q1.statespace != Q2.statespace
        error("Cannot multiply QExpr’s from different statespaces.")
    end
    new_terms::AbstractVector{QComposite} = []
    for t1 in Q1.terms
        for t2 in Q2.terms
            append!(new_terms, t1 * t2)
        end
    end
    return QExpr(Q1.statespace, new_terms)
end
function *(Q1::QExpr, Q2::T)::QExpr where T<:QComposite
    terms::Vector{QComposite} = QComposite[]
    for q in Q1.terms 
        append!(terms, q*Q2)
    end
    return QExpr(Q1.statespace, terms)
end
function *(Q1::T, Q2::QExpr)::QExpr where T<:QComposite
    terms::Vector{QComposite} = QComposite[]
    for q in Q2.terms
        append!(terms, Q1 * q) 
    end
    return QExpr(Q1.statespace, terms)
end

function *(Q1::QExpr, num::Number)::QExpr
    if num == 0
        return QExpr(Q1.statespace, QComposite[])  # or Vector{QComposite}() if preferred
    end
    terms::Vector{QComposite} = QComposite[]
    for q in Q1.terms 
        append!(terms, q*num )
    end
    return QExpr(Q1.statespace, terms)
end
*(num::Number, Q1::QExpr)::QExpr = Q1 * num


function *(Q1::T, num::Number)::Vector{QComposite} where T<:QComposite
    return [modify_expr(Q1, Q1.expr * num)]
end
(*(num::Number, Q2::T)::Vector{QComposite}) where T<:QComposite = Q2 * num

##### Exponentiation ###################################################
function ^(Q::Union{QExpr,T}, n::Integer) where T<:QComposite
    if n < 0
        error("Negative exponent not defined for subtype of QComposite or QExpr.")
    elseif n == 0
        return Identity(Q.statespace)
    end

    result = Identity(Q.statespace)
    base = Q
    exp = n
    while exp > 0
        if isodd(exp)
            result = result * base
        end
        base = base * base
        exp = exp ÷ 2
    end
    return result
end

##### Commutators ###################################################
# --- Define the commutator for QExpr ---
"""
    Commutator(Q1::Union{QExpr, s}, Q2::Union{QExpr, T}) where {S <: QExpr, T <: QExpr}

Computes the commutator [Q1, Q2] = Q1 * Q2 - Q2 * Q1.
Both QExpr's must share the same statespace.
"""
function Commutator(Q1::S, Q2::T)::Vector{QComposite} where {S<:QComposite,T<:QComposite}
    # QExpr multiplication is already defined.
    return vcat(Q1 * Q2, .-Q2 * Q1)
end
(Commutator(Q1::QExpr, Q2::T)::QExpr) where {T<:QComposite} = Q1 * Q2 - Q2 * Q1
(Commutator(Q1::T, Q2::QExpr)::QExpr) where {T<:QComposite} = Q1 * Q2 - Q2 * Q1
Commutator(Q1::QExpr, Q2::QExpr)::QExpr = Q1 * Q2 - Q2 * Q1

function +(Q::QExpr, exprs::Vector{QExpr})::QExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return c + Q
end
function +(exprs::Vector{QExpr}, Q::QExpr)::QExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return c + Q
end
function -(Q::QExpr, exprs::Vector{QExpr})::QExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return Q - c
end
function -(exprs::Vector{QExpr}, Q::QExpr)::QExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return c - Q
end
function *(Q::QExpr, exprs::Vector{QExpr})::QExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return Q * c
end
function *(exprs::Vector{QExpr}, Q::QExpr)::QExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return c * Q
end



function Identity(qspace::StateSpace)
    return QExpr(qspace, QAtomProduct(qspace, qspace.c_one, QAtom[QTerm(qspace.I_op)]))
end

"""
    Dag(qspace::StateSpace, t::QTerm) -> QTerm
    Dag(t::QExpr) -> QExpr
    Dag(t::T) -> T where T <: QComposite
    
Returns the Hermitian conjugate (dagger) of a QTerm, QExpr or QSum.
Overloads the adjoint function, which can be called via `t′`.
"""
function Dag(qspace::StateSpace, t::QTerm)::Vector{Tuple{QTerm,ComplexRational}}
    # Build a vector of the op_set for each operator factor, in the same order as t.op_indices.
    new_op_inds::Vector{Vector{Tuple{ComplexRational,Is}}} = []
    curr_op_inds = t.op_indices
    i = 1
    for subspace in qspace.subspaces
        for _ in 1:length(subspace.keys)
            op = curr_op_inds[i]
            push!(new_op_inds, subspace.op_set.op_dag(op))
            i += 1
        end
    end
    # Construct terms for each combination of new_op_inds
    terms::Vector{QTerm} = QTerm[]
    coeffs::Vector{ComplexRational} = ComplexRational[]
    for combo in product(new_op_inds...)
        curr_inds::Vector{Vector{Int}} = []
        curr_coeff::ComplexRational = ComplexRational(1, 0, 1)
        for i in 1:length(combo)
            curr_coeff *= combo[i][1]
            push!(curr_inds, combo[i][2])
        end
        push!(terms, QTerm(curr_inds))
        push!(coeffs, curr_coeff)
    end
    return collect(zip(terms, coeffs))
end
function Dag(qspace::StateSpace, t::QAbstract)::Vector{Tuple{QAbstract,ComplexRational}}
    # check if is hermitian
    if t.operator_type.hermitian
        return Tuple{QAbstract,ComplexRational}[(t, one(ComplexRational))]
    end
    new_t = dag_copy(t)
    return Tuple{QAbstract,ComplexRational}[(new_t, one(ComplexRational))]
end

function Dag(p::QAtomProduct)::Vector{QAtomProduct}
    # reverse order elementwise Dag 
    terms::Vector{Vector{Tuple{QAtom,ComplexRational}}} = []
    for t in reverse(p.expr)
        push!(terms, Dag(p.statespace, t))
    end
    coeff_fun = p.coeff_fun
    new_atom_products::Vector{QAtomProduct} = []
    for combo in Iterators.product(terms...)
        curr_atoms = [a[1] for a in combo]
        curr_coeff = [a[2] for a in combo]
        coeff = reduce(*, curr_coeff)
        if !iszero(coeff)
            push!(new_atom_products, QAtomProduct(p.statespace, coeff_fun * coeff, curr_atoms, p.separate_expectation_values))
        end
    end
    return new_atom_products
end

function Dag(Q::QExpr)::QExpr
    q_comps::Vector{QComposite} = []
    for q in Q.terms
        append!(q_comps, Dag(q))
    end
    return QExpr(q_comps)
end
function Dag(t::T)::Vector{QComposite} where T<:QComposite
    t_new = copy(t)
    t_new.expr = Dag(t.expr)
    return QComposite[t_new]
end
function Dag(t::QMultiComposite)::Vector{QMultiComposite}
    dag_exprs::Vector{QExpr} = []
    for expr in reverse(t.expr)
        append!(dag_exprs, Dag(expr))
    end
    return [modify_expr(t, dag_exprs)]
end

adjoint(Q::QExpr) = Dag(Q)
adjoint(Q::T) where {T<:QComposite} = Dag(Q)
