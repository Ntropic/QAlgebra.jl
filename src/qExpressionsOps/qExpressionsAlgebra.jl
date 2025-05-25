export Dag, Commutator, is_numeric

""" 
    is_numeric(t::qTerm, qspace::StateSpace) -> Bool
    is_numeric(t::qAbstract, qspace::StateSpace) -> Bool
    is_numeric(p::qProd) -> Bool
    is_numeric(s::qSum) -> Bool
    is_numeric(expr::qEQ) -> Bool

Returns true if only the coefficient of the term(s) is non-zero.
"""
function is_numeric(e::qAtom, statespace::StateSpace)
    error("is_numeric (with given statespace) not implemented for qAbstract subtype $(typeof(e))")
end
function is_numeric(e::qComposite)
    error("is_numeric (without given statespace) not implemented for qComposite subtype $(typeof(e))")
end

function is_numeric(t::qTerm, qspace::StateSpace)::Bool
    return qspace.neutral_op == t.op_indices
end
function is_numeric(t::qAbstract, qspace::StateSpace)::Bool
    return false 
end

function is_numeric(p::qProd)::Bool
    return all(is_numeric(t, p.statespace) for t in p.terms)
end
function is_numeric(s::qSum)::Bool
    return is_numeric(s.expr)
end
function is_numeric(expr::qEQ)::Bool
    expr_s = simplify(expr)
    terms = expr_s.terms

    if isempty(terms)
        return true  # No terms = numeric 0
    elseif length(terms) == 1
        t = terms[1]
        if t isa qTerm
            return is_numeric(t, expr_s.statespace)
        else
            return false2
        end
    else
        return false
    end
end

function isapprox_num(x, y; atol=1e-12)
    return isapprox(x, y, atol=atol)
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

function ==(a::qEQ, b::qEQ)
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
function ==(expr::qEQ, n::Number)
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
function ==(n::Number, expr::qEQ)
    return expr == n  # Symmetric
end


### Basic Operations: 
function -(t::qTerm)::qTerm
    return qTerm(-t.coeff, t.var_exponents, copy(t.op_indices))
end
function -(t::qEQ)::qEQ
    return qEQ(.-t.terms, t.statespace)
end
function -(t::qSum)::qSum
    return qSum(-t.expr, t.indexes, t.subsystem_index, t.element_indexes, t.neq)
end

function +(Q1::qEQ, Q2::qEQ)::qEQ
    if Q1.statespace != Q2.statespace
        error("Cannot add qEQ’s from different statespaces.")
    end
    new_terms = vcat(Q1.terms, Q2.terms)
    # Optionally: group like terms here.
    return qEQ(new_terms, Q1.statespace)
end
function +(Q1::qEQ, term::qTerm)::qEQ
    new_terms = vcat(Q1.terms, [term])
    return qEQ(new_terms, Q1.statespace)
end
function +(term::qTerm, Q2::qEQ)::qEQ
    new_terms = vcat([term], Q2.terms)
    return qEQ(new_terms, Q2.statespace)
end
function +(Q::qEQ, S::qSum)::qEQ
    if Q.statespace != S.expr.statespace
        error("Cannot add qEQ’s from different statespaces.")
    end
    new_terms = vcat(Q.terms, S)
    return qEQ(new_terms, Q.statespace)
end
function +(S::qSum, Q::qEQ)::qEQ
    if Q.statespace != S.expr.statespace
        error("Cannot add qEQ’s from different statespaces.")
    end
    new_terms = vcat(S, Q.terms)
    return qEQ(new_terms, Q.statespace)
end
function +(S1::qSum, S2::qSum)::qEQ
    if S1.expr.statespace != S2.expr.statespace
        error("Cannot add qSum’s from different statespaces.")
    end
    # check if the following qSum parameters are the same: index, subsystem_index, element_index
    if S1.indexes == S2.indexes || S1.subsystem_index == S2.subsystem_index || S1.element_indexes == S2.element_indexes || S1.neq == S2.neq
        S_sum_expr = S1.expr + S2.expr
        S_sum = qSum(S_sum_expr, S1.indexes, S1.subsystem_index, S1.element_indexes, S1.neq)
        return qEQ([S_sum], S1.expr.statespace)
    end
    new_terms = vcat(S1, S2)
    return qEQ(new_terms, S1.expr.statespace)
end
# --- Define qEQ subtraction as addition of the negative ---
function -(Q1::qEQ, Q2::qEQ)
    Q2_minus = -Q2
    return Q1 + Q2_minus
end
function -(Q::qEQ, S::qSum)::qEQ
    S_minus = -S
    return Q + S_minus
end
function -(S::qSum, Q::qEQ)::qEQ
    Q_minus = -Q
    return S + Q_minus
end
function -(S1::qSum, S2::qSum)::qEQ
    S2_minus = -S2
    return S1 + S2_minus
end


"""
    Dag(t::qTerm, qspace::StateSpace) -> qTerm
    Dag(t::qEQ) -> qEQ
    Dag(t::qSum) -> qSum
    
Returns the Hermitian conjugate (dagger) of a qTerm, qEQ or qSum.
Overloads the adjoint function, which can be called via `t′`.
"""
function Dag(t::qTerm, qspace::StateSpace)::qEQ
    new_coeff = conj(t.coeff)
    new_exponents = copy(t.var_exponents)
    # Build a vector of the op_set for each operator factor, in the same order as t.op_indices.
    new_op_inds::Vector{Vector{Tuple{Number,Is}}} = []
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
    for combo in product(new_op_inds...)
        curr_inds::Vector{Is} = []
        curr_coeff::Number = new_coeff
        for i in 1:length(combo)
            curr_coeff *= combo[i][1]
            push!(curr_inds, combo[i][2])
        end
        push!(terms, qTerm(curr_coeff, copy(new_exponents), curr_inds))
    end
    return qEQ(terms, qspace)
end
function Dag(Q::qEQ)::qEQ
    return sum([Dag(t, Q.statespace) for t in Q.terms])
end
function Dag(t::qSum)::qSum
    return qSum(Dag(t.expr), t.indexes, t.subsystem_index, t.element_indexes, t.neq)
end
function Dag(t::qSum, qspace::StateSpace)::qSum
    return qSum(Dag(t.expr), t.indexes, t.subsystem_index, t.element_indexes, t.neq)
end

adjoint(Q::qEQ) = Dag(Q)
adjoint(Q::qSum) = Dag(Q)

"""
    multiply_qterm(t1::qTerm, t2::qTerm, statespace::StateSpace) -> Vector{qTerm}

Multiplies two qTerm’s from the same statespace:
  - Coefficients are multiplied.
  - var_exponents are added elementwise.
  - For each subspace, the corresponding operator indexes are combined using the subspace’s op_set.op_product.
    This returns a vector of tuples (factor, index). We then take the Cartesian product over subspaces;
    each combination yields a new qTerm with its coefficient multiplied by the product of the factors,
    and with op_indices given by the corresponding indexes.
"""
# no docstring for this function
function multiply_qterm(t1::qTerm, t2::qTerm, statespace::StateSpace)::Vector{qTerm}
    # Multiply coefficients and add variable exponents.
    new_coeff_base = t1.coeff * t2.coeff
    new_exponents = [a + b for (a, b) in zip(t1.var_exponents, t2.var_exponents)]

    op_indices1 = t1.op_indices
    op_indices2 = t2.op_indices

    # For the operator indexes, iterate over the subspaces.
    results_per_subspace::Vector{Vector{Tuple{Number,Is}}} = Vector{Tuple{Number,Is}}[]
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
    # Iterate over the Cartesian product.
    for combo in product(results_per_subspace...)
        # combo is a tuple, one element per subspace.
        new_factor = new_coeff_base
        new_term::Vector{Is} = []
        for c in combo
            new_factor *= c[1]
            push!(new_term, c[2])
        end
        push!(new_terms, qTerm(new_factor, copy(new_exponents), new_term))
    end

    return new_terms
end
function *(t1::qTerm, t2::qTerm, statespace::StateSpace)
    return multiply_qterm(t1, t2, statespace)
end
function *(Q1::qEQ, Q2::qEQ)::qEQ
    if Q1.statespace != Q2.statespace
        error("Cannot multiply qEQ’s from different statespaces.")
    end
    new_terms = qTerm[]
    for t1 in Q1.terms
        for t2 in Q2.terms
            prod_terms = multiply_qterm(t1, t2, Q1.statespace)
            append!(new_terms, prod_terms)
        end
    end
    return qEQ(new_terms, Q1.statespace)
end
function *(Q1::qEQ, num::Number)::qEQ
    return qEQ([qTerm(num * t.coeff, t.var_exponents, t.op_indices) for t in Q1.terms], Q1.statespace)
end
function *(s::Number, Q::qEQ)::qEQ
    new_terms = [qTerm(s * t.coeff, t.var_exponents, t.op_indices) for t in Q.terms]
    return qEQ(new_terms, Q.statespace)
end
function *(S::qSum, num::Number)::qSum
    return qSum(S.expr * num, S.indexes, S.subsystem_index, S.element_indexes, S.neq)
end
function *(num::Number, S::qSum)::qSum
    return qSum(S.expr * num, S.indexes, S.subsystem_index, S.element_indexes, S.neq)
end
function *(S::qSum, Q::qEQ)::qSum
    if S.expr.statespace != Q.statespace
        error("Cannot multiply qSum and qEQ from different statespaces.")
    end
    new_terms = S.expr * Q
    return qEQ([qSum(new_terms, S.indexes, S.subsystem_index, S.element_indexes)], Q.statespace)
end
function *(Q::qEQ, S::qSum)::qSum
    if S.expr.statespace != Q.statespace
        error("Cannot multiply qSum and qEQ from different statespaces.")
    end
    new_terms = Q * S.expr
    return qSum(new_terms, S.indexes, S.subsystem_index, S.element_indexes)
end
function *(Q::qSum, S::qSum)::qSum
    if Q.expr.statespace != S.expr.statespace
        error("Cannot multiply qSum and qEQ from different statespaces.")
    end
    new_terms = Q.expr * S.expr
    inner_sum = qSum(new_terms, S.indexes, S.subsystem_index, S.element_indexes, S.neq)
    outer_sum = qSum(qEQ([inner_sum], Q.expr.statespace), Q.index, Q.subsystem_index, Q.element_indexes, Q.neq)
    return qEQ([outer_sum], Q.expr.statespace)
end

# --- Define the commutator for qEQ ---
"""
    Commutator(Q1::qEQ, Q2::qEQ) -> qEQ
    Commutator(Q::qEQ, t::qSum) -> qEQ
    Commutator(t::qSum, Q::qEQ) -> qEQ

Computes the commutator [Q1, Q2] = Q1 * Q2 - Q2 * Q1.
Both qEQ's must share the same statespace.
"""
function Commutator(Q1::qEQ, Q2::qEQ)::qEQ
    # qEQ multiplication is already defined.
    return Q1 * Q2 - Q2 * Q1
end
function Commutator(Q::qEQ, t::qSum)::qSum
    return Q * t - t * Q
end
function Commutator(t::qSum, Q::qEQ)::qSum
    return t * Q - Q * t
end
function Commutator(Q::qSum, t::qSum)::qSum
    return Q * t - t * Q
end

function Commutator(t1::qTerm, t2::qTerm, statespace::StateSpace)::qEQ
    Q1 = qEQ([t1], statespace)
    Q2 = qEQ([t2], statespace)
    return Commutator(Q1, Q2)
end

function identity_qEQ(qspace::StateSpace)
    # Set variable exponents to zero.
    var_exponents = zeros(Int, length(qspace.vars))
    # Build a vector of neutral operator indexes.
    neutral_ops = [s.op_set.neutral_element for s in qspace.subspaces for key in s.keys]
    # Create a single-term qEQ.
    return qEQ([qTerm(1, var_exponents, neutral_ops)], qspace)
end

function ^(Q::qEQ, n::Integer)
    if n < 0
        error("Negative exponent not defined for qEQ.")
    elseif n == 0
        return identity_qEQ(Q.statespace)
    end

    result = identity_qEQ(Q.statespace)
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



# Addition of a diff_qEQ with a qEQ.
function +(d::diff_qEQ, Q::Union{qEQ,qExpr})::diff_qEQ
    new_rhs = simplify(d.right_hand_side + Q)
    return diff_qEQ(d.left_hand_side, new_rhs, d.statespace; braket=d.braket, do_sigma=d.do_sigma)
end
function -(d::diff_qEQ, Q::Union{qEQ,qExpr})::diff_qEQ
    new_rhs = simplify(d.right_hand_side - Q)
    return diff_qEQ(d.left_hand_side, new_rhs, d.statespace; braket=d.braket, do_sigma=d.do_sigma)
end
function *(d::diff_qEQ, Q::Union{qEQ,qExpr})::diff_qEQ
    new_rhs = simplify(d.right_hand_side * Q)
    return diff_qEQ(d.left_hand_side, new_rhs, d.statespace; braket=d.braket, do_sigma=d.do_sigma)
end

##### Commutators ###################################################
function +(Q::qEQ, exprs::Vector{qEQ})::qEQ
    if length(exprs) != 2
        error("Only vectors of length 2 can be added -> commutator.")
    end
    return Q + Commutator(exprs[1], exprs[2])
end
function -(Q::qEQ, exprs::Vector{qEQ})::qEQ
    if length(exprs) != 2
        error("Only vectors of length 2 can be added -> commutator.")
    end
    return Q - Commutator(exprs[1], exprs[2])
end
function *(Q::qEQ, exprs::Vector{qEQ})::qEQ
    if length(exprs) != 2
        error("Only vectors of length 2 can be added -> commutator.")
    end
    return Q * Commutator(exprs[1], exprs[2])
end
function +(d::diff_qEQ, exprs::Vector{qEQ})::diff_qEQ
    if length(exprs) != 2
        error("Only vectors of length 2 can be added -> commutator.")
    end
    return d + Commutator(exprs[1], exprs[2])
end
function -(d::diff_qEQ, exprs::Vector{qEQ})::diff_qEQ
    if length(exprs) != 2
        error("Only vectors of length 2 can be added -> commutator.")
    end
    return d - Commutator(exprs[1], exprs[2])
end
function *(d::diff_qEQ, exprs::Vector{qEQ})::diff_qEQ
    if length(exprs) != 2
        error("Only vectors of length 2 can be added -> commutator.")
    end
    return d * Commutator(exprs[1], exprs[2])
end