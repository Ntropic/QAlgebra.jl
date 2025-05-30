export Dag, Commutator, is_numeric, same_statespace

""" 
    is_numeric(t::qTerm, qspace::StateSpace) -> Bool
    is_numeric(t::qAbstract, qspace::StateSpace) -> Bool
    is_numeric(p::qProd) -> Bool
    is_numeric(s::qSum) -> Bool
    is_numeric(expr::qExpr) -> Bool

Returns true if only the coefficient of the term(s) is non-zero.
"""
function is_numeric(e::qAtom, statespace::StateSpace)
    error("is_numeric (with given statespace) not implemented for qAbstract subtype $(typeof(e))")
end
function is_numeric(e::qComposite)
    error("is_numeric (without given statespace) not implemented for qComposite subtype $(typeof(e))")
end

function is_numeric(t::qTerm, statespace::StateSpace)::Bool
    return statespace.neutral_op == t.op_indices
end
function is_numeric(t::qAbstract, statespace::StateSpace)::Bool
    return false 
end
function is_numeric(p::qProd)::Bool
    return all(is_numeric(t, p.statespace) for t in p.terms)
end
function is_numeric(s::qSum)::Bool
    return is_numeric(s.expr)
end
function is_numeric(expr::qExpr)::Bool
    expr_s = simplify(expr)
    terms = expr_s.terms

    if isempty(terms)
        return true  # No terms = numeric 0
    else
        return all(is_numeric(t, expr_s.statespace) for t in terms)
    end
end

function where_neutral(q::qAtom, statespace::StateSpace)::Vector{Bool}
    return [op==neut for (op, neut) in zip(q.op_indices, statespace.neutral_op)]
end
function where_neutral(q::qAbstract, statespace::StateSpace)::Vector{Bool}
    return copy(q.operator_type.subspaces)
end
function where_acting(q::qAtom, statespace::StateSpace)::Vector{Bool}
    return [op!=neut for (op, neut) in zip(q.op_indices, statespace.neutral_op)]
end
function where_acting(q::qAbstract, statespace::StateSpace)::Vector{Bool}
    return copy(q.operator_type.non_subspaces)
end

function qAtom_commute(q1::qAtom, q2::qAtom, statespace::StateSpace)::Bool
    # check if all elements of where neutral are NAND
    act_q1 = where_acting(q1, statespace)
    act_q2 = where_acting(q2, statespace)
    return all([nand(a,b) for (a,b) in zip(act_q1, act_q2)])
end

function same_statespace(a::qComposite, b::qComposite)::Bool
    return a.statespace == b.statespace
end

import Base: ==

function ==(a::qTerm, b::qTerm)
    return a.op_indices == b.op_indices
end
function ==(a::qAbstract, b::qAbstract)
    return a.key_index == b.key_index && a.sub_index == b.sub_index && a.exponent == b.exponent && a.dag == b.dag && a.index_map == b.index_map 
end
function ==(a::qProd, b::qProd)
    return a.coeff_fun == b.coeff_fun &&  all([ai == bi for (ai, bi) in zip(a.expr, b.expr)])
end

function ==(a::qExpr, b::qExpr)
    simple_a = simplify(a)
    simple_b = simplify(b)
    if length(simple_a) != length(simple_b)
        return false
    end
    if simple_a.statespace != simple_b.statespace
        return false
    end
    return all([ai == bi for (ai, bi) in zip(simple_a, simple_b)])
end
function ==(a::qSum, b::qSum)
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
function ==(expr::qExpr, n::Number)
    function isapprox_num(x, y; atol=1e-12)
        return isapprox(x, y, atol=atol)
    end
    simple_expr = simplify(expr)
    if is_numeric(simple_expr)
        if length(simple_expr.terms) == 0
            return isapprox_num(0, n)
        else
            return isapprox_num(simple_expr.terms[1].coeff, n)
        end
    end
    return false
end
function ==(n::Number, expr::qExpr)
    return expr == n  # Symmetric
end


### Basic Operations: 

####### Unary Minus #############################################################
function -(t::qProd)::qProd
    return qProd(t.statespace, -t.coeff_fun, copy(t.expr))
end
function -(t::qExpr)::qExpr
    return qExpr(.-t.terms, t.statespace)  
end
function -(t::qComposite)::qComposite
    t_new = copy(t)
    t_new.expr = -t_new.expr  # Negate the expression inside the composite.
    return t_new 
end

#### Binary + ####################################################################
function +(Q1::qExpr, Q2::qExpr)::qExpr
    new_terms = vcat(Q1.terms, Q2.terms)
    # Optionally: group like terms here.
    return qExpr(Q1.statespace, new_terms)
end
function +(Q1::qExpr, Q2::qComposite)::qExpr
    new_terms = copy(Q1.terms)
    push!(new_terms, Q2)
    return qExpr(Q1.statespace, new_terms)
end

#### Binary - ####################################################################
function -(Q1::qExpr, Q2::qExpr)::qExpr
    new_terms = vcat(Q1.terms, .-Q2.terms)
    # Optionally: group like terms here.
    return qExpr(Q1.statespace, new_terms)
end
function -(Q1::qExpr, Q2::qComposite)::qExpr
    new_terms = copy(Q1.terms)
    push!(new_terms, .-Q2)
    return qExpr(Q1.statespace, new_terms)
end

#### Multiply ####################################################################
function trivial_multiply(q1::qProd, q2::qProd)::qProd
    return qProd(Q1.statespace, Q1.coeff_fun*Q2.coeff_fun, vcat(Q1.terms, Q2.terms))
end
# Multiplies two qTerm’s from the same statespace. Returns a vector of qTerm’s that are the result of this multiplication and corresponding ComplexRational coefficients. 
function multiply_qterm(t1::qTerm, t2::qTerm, statespace::StateSpace)::Tuple{Vector{qTerm}, Vector{ComplexRational}}
    # Multiply coefficients and add variable exponents.
    coeff_base = ComplexRational(1,0,1)

    # For the operator indexes, iterate over the subspaces.
    results_per_subspace::Vector{Vector{Tuple{ComplexRational,Is}}} = Vector{Tuple{ComplexRational,Is}}[]
    # We assume the length of op_indices equals the number of subspaces.
    i = 0
    for subspace in statespace.subspaces
        op_set = subspace.op_set
        for (ind, key) in zip(subspace.statespace_inds, subspace.keys)
            i += 1
            op1 = t1.op_indices[i]
            op2 = t2.op_indices[i]
            # op_product should return a vector of tuples (factor, index).
            res = op_set.op_product(op1, op2)
            push!(results_per_subspace, res)
        end
    end

    new_terms = qTerm[]
    new_coeffs = ComplexRational[]
    # Iterate over the Cartesian product.
    for combo in product(results_per_subspace...)
        # combo is a tuple, one element per subspace.
        new_factor = coeff_base
        new_term::Vector{Is} = []
        for c in combo
            new_factor *= c[1]
            push!(new_term, c[2])
        end
        push!(new_terms, qTerm(new_term))
        push!(new_coeffs, new_factor)
    end
    return new_terms, new_coeffs
end
function multiply_qterms(terms::Vector{qTerm}, statespace::StateSpace)::Tuple{Vector{qTerm}, Vector{ComplexRational}}
    n = length(terms)
    if n == 1
        return terms, [one(ComplexRational)]
    elseif n == 2
        return multiply_qterms(terms[1], terms[2], statespace)
    elseif n > 2 
        curr_terms, curr_coeffs = multiply_qterms(terms[1:end-1], statespace)
        final_terms::Vector{qTerm} = qTerm[]
        final_coeffs::Vector{ComplexRational} = ComplexRational[]
        for i in 1:length(curr_terms)
            new_terms, new_coeffs = multiply_qterms([curr_terms[i]], terms[end], statespace)
            append!(final_terms, new_terms)
            append!(final_coeffs, new_coeffs .* curr_coeffs[i])
        end
        return final_terms, final_coeffs
    else
        error("Invalid number of terms in multiply_qterms.")
    end
end
function remove_subsequent_qTerms(p::qProd)::Vector{qProd}
    # first group terms into qTerm boxes and qAbstract boxes 
    ss = p.statespace
    terms::Vector{Vector{qAtom}} = []
    coeffs::Vector{ComplexRational} = []
    i = 1 
    curr_terms::Vector{qAtom} = []
    for term in p.expr
        if istype(term, qTerm) 
            push!(curr_terms, term)
        else 
            if length(curr_terms) > 0 
                multiplied_terms, multiplied_coeffs = multiply_qterms(curr_terms, ss)
                curr_terms = Vector{qAtom}[]
                new_coeffs = Vector{ComplexRational}[]
                for (t, c) in zip(terms, coeffs)
                    for (mt, mc) in zip(multiplied_terms, multiplied_coeffs)
                        push!(new_terms, vcat(t, mt, terms))
                        push!(new_coeffs, c * mc)
                    end
                end
                terms = new_terms 
                coeffs = new_coeffs 
            else
               for i in 1:length(terms)
                   push!(terms[i], term)
                end
            end
        end
    end
    return qProd[qProd(ss, p.coeff_fun*c, t) for (t, c) in zip(terms, coeffs)]
end
function commuting_sort(p::qProd)::qProd
    ss = p.statespace 
    terms = p.expr
    coeff_fun = p.coeff_fun

    n = length(terms) 
    # create sort keys for each term 
    sort_keys = [sort_key(t) for t in terms]
    is_qTerm::Vector{Bool} = [isa(t, qTerm) for t in terms]
    i = 1   
    while i < n # reduce n when removing an element 
        # check if both are qTerm -> unify them 
        # otherwise use qAtom_commute to check if they commute 
    end
end
            

function *(p1::qProd, p2::qProd)::qProd
    p = trivial_multiply(p1, p2)
    # now lets simplify by first sorting the qAbstract elements as much as possible and then simplifying subsequent qTerm terms.
end

function *(t1::qTerm, t2::qTerm, statespace::StateSpace)
    return multiply_qatom(t1, t2, statespace)
end
function *(Q1::qExpr, Q2::qExpr)::qExpr
    if Q1.statespace != Q2.statespace
        error("Cannot multiply qExpr’s from different statespaces.")
    end
    new_terms = qTerm[]
    for t1 in Q1.terms
        for t2 in Q2.terms
            prod_terms = multiply_qatom(t1, t2, Q1.statespace)
            append!(new_terms, prod_terms)
        end
    end
    return qExpr(new_terms, Q1.statespace)
end
function *(Q1::qExpr, num::Number)::qExpr
    return qExpr([qTerm(num * t.coeff, t.var_exponents, t.op_indices) for t in Q1.terms], Q1.statespace)
end
function *(s::Number, Q::qExpr)::qExpr
    new_terms = [qTerm(s * t.coeff, t.var_exponents, t.op_indices) for t in Q.terms]
    return qExpr(new_terms, Q.statespace)
end
function *(S::qSum, num::Number)::qSum
    return qSum(S.expr * num, S.indexes, S.subsystem_index, S.element_indexes, S.neq)
end
function *(num::Number, S::qSum)::qSum
    return qSum(S.expr * num, S.indexes, S.subsystem_index, S.element_indexes, S.neq)
end
function *(S::qSum, Q::qExpr)::qSum
    if S.expr.statespace != Q.statespace
        error("Cannot multiply qSum and qExpr from different statespaces.")
    end
    new_terms = S.expr * Q
    return qExpr([qSum(new_terms, S.indexes, S.subsystem_index, S.element_indexes)], Q.statespace)
end
function *(Q::qExpr, S::qSum)::qSum
    if S.expr.statespace != Q.statespace
        error("Cannot multiply qSum and qExpr from different statespaces.")
    end
    new_terms = Q * S.expr
    return qSum(new_terms, S.indexes, S.subsystem_index, S.element_indexes)
end
function *(Q::qSum, S::qSum)::qSum
    if Q.expr.statespace != S.expr.statespace
        error("Cannot multiply qSum and qExpr from different statespaces.")
    end
    new_terms = Q.expr * S.expr
    inner_sum = qSum(new_terms, S.indexes, S.subsystem_index, S.element_indexes, S.neq)
    outer_sum = qSum(qExpr([inner_sum], Q.expr.statespace), Q.index, Q.subsystem_index, Q.element_indexes, Q.neq)
    return qExpr([outer_sum], Q.expr.statespace)
end




##### Exponentiation ###################################################
function ^(Q::qExpr, n::Integer)
    if n < 0
        error("Negative exponent not defined for qExpr.")
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
# --- Define the commutator for qExpr ---
"""
    Commutator(Q1::qExpr, Q2::qExpr) -> qExpr
    Commutator(Q::qExpr, t::qSum) -> qExpr
    Commutator(t::qSum, Q::qExpr) -> qExpr

Computes the commutator [Q1, Q2] = Q1 * Q2 - Q2 * Q1.
Both qExpr's must share the same statespace.
"""
function Commutator(Q1::qExpr, Q2::qExpr)::qExpr
    # qExpr multiplication is already defined.
    return Q1 * Q2 - Q2 * Q1
end
function Commutator(Q::qExpr, t::qSum)::qSum
    return Q * t - t * Q
end
function Commutator(t::qSum, Q::qExpr)::qSum
    return t * Q - Q * t
end
function Commutator(Q::qSum, t::qSum)::qSum
    return Q * t - t * Q
end

function Commutator(t1::qTerm, t2::qTerm, statespace::StateSpace)::qExpr
    Q1 = qExpr([t1], statespace)
    Q2 = qExpr([t2], statespace)
    return Commutator(Q1, Q2)
end

function Identity(qspace::StateSpace)
    # Set variable exponents to zero.
    var_exponents = zeros(Int, length(qspace.vars))
    # Build a vector of neutral operator indexes.
    neutral_ops = [s.op_set.neutral_element for s in qspace.subspaces for key in s.keys]
    # Create a single-term qExpr.
    return qExpr(qspace, qProd(qspace, qspace.fone, qspace.neutral_op))
end

function +(Q::qExpr, exprs::Vector{qExpr})::qExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be added -> commutator.")
    end
    return Q + Commutator(exprs[1], exprs[2])
end
function -(Q::qExpr, exprs::Vector{qExpr})::qExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be added -> commutator.")
    end
    return Q - Commutator(exprs[1], exprs[2])
end
function *(Q::qExpr, exprs::Vector{qExpr})::qExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be added -> commutator.")
    end
    return Q * Commutator(exprs[1], exprs[2])
end
function +(d::diff_qEQ, exprs::Vector{qExpr})::diff_qEQ
    if length(exprs) != 2
        error("Only vectors of length 2 can be added -> commutator.")
    end
    return d + Commutator(exprs[1], exprs[2])
end
function -(d::diff_qEQ, exprs::Vector{qExpr})::diff_qEQ
    if length(exprs) != 2
        error("Only vectors of length 2 can be added -> commutator.")
    end
    return d - Commutator(exprs[1], exprs[2])
end
function *(d::diff_qEQ, exprs::Vector{qExpr})::diff_qEQ
    if length(exprs) != 2
        error("Only vectors of length 2 can be added -> commutator.")
    end
    return d * Commutator(exprs[1], exprs[2])
end


"""
    Dag(t::qTerm, qspace::StateSpace) -> qTerm
    Dag(t::qExpr) -> qExpr
    Dag(t::qComposite) -> qComposite
    
Returns the Hermitian conjugate (dagger) of a qTerm, qExpr or qSum.
Overloads the adjoint function, which can be called via `t′`.
"""
function Dag(t::qTerm, qspace::StateSpace)::Tuple{Vector{qTerm}, Vector{ComplexRational}}
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
    terms::Vector{qTerm} = qTerm[]
    coeffs::Vector{ComplexRational} = ComplexRational[]
    for combo in product(new_op_inds...)
        curr_inds::Vector{Is} = []
        curr_coeff::ComplexRational = new_coeff
        for i in 1:length(combo)
            curr_coeff *= combo[i][1]
            push!(curr_inds, combo[i][2])
        end
        push!(terms, qTerm(curr_inds))
        push!(coeffs, curr_coeff)
    end
    return terms, coeffs
end
function Dag(p::qProd)::qExpr 
    # continue here 
end

function Dag(Q::qExpr)::qExpr
    return sum([Dag(t, Q.statespace) for t in Q.terms])
end
function Dag(t::qSum)::qSum
    return qSum(Dag(t.expr), t.indexes, t.subsystem_index, t.element_indexes, t.neq)
end
function Dag(t::qSum, qspace::StateSpace)::qSum
    return qSum(Dag(t.expr), t.indexes, t.subsystem_index, t.element_indexes, t.neq)
end

adjoint(Q::qExpr) = Dag(Q)
adjoint(Q::qSum) = Dag(Q)

