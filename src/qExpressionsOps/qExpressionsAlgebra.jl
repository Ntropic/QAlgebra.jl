export Dag, Commutator, is_numeric, same_statespace

""" 
    is_numeric(t::qTerm, qspace::StateSpace) -> Bool
    is_numeric(t::qAbstract, qspace::StateSpace) -> Bool
    is_numeric(p::qAtomProduct) -> Bool
    is_numeric(s::qSum) -> Bool
    is_numeric(expr::qExpr) -> Bool

Returns true if only the coefficient of the term(s) is non-zero.
"""
function is_numeric(op_indices::Vector{Is}, statespace::StateSpace)::Bool
    return statespace.neutral_op == op_indices
end
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
function is_numeric(p::qAtomProduct)::Bool
    return all(is_numeric(t, p.statespace) for t in p.expr)
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
function ==(a::qAtomProduct, b::qAtomProduct)
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
function -(t::qAtomProduct)::qAtomProduct
    return qAtomProduct(t.statespace, -t.coeff_fun, copy(t.expr))
end
function -(t::qExpr)::qExpr
    return qExpr(t.statespace, .-t.terms)  
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
function trivial_multiply(Q1::qAtomProduct, Q2::qAtomProduct)::qAtomProduct
    return qAtomProduct(Q1.statespace, Q1.coeff_fun*Q2.coeff_fun, vcat(Q1.expr, Q2.expr))
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
function *(t1::qTerm, t2::qTerm, statespace::StateSpace)  # should never be called by the user!!
    return multiply_qterm(t1, t2, statespace)
end

#### Main Multiplication Functions ################################################
function *(p1::qAtomProduct, p2::qAtomProduct)::Vector{qAtomProduct}
    p = trivial_multiply(p1, p2)  # append the terms of p1 and p2.
    return simplify(p)    # simplify the product. 
end
function *(p1::qAtomProduct, num::Number)::Vector{qAtomProduct}
    return [qAtomProduct(p1.statespace, p1.coeff_fun*num, copy(p1.expr))]
end
function *(num::Number, p1::qAtomProduct)::Vector{qAtomProduct}
    return [qAtomProduct(p1.statespace, p1.coeff_fun*num, copy(p1.expr))]
end
function *(p1::T1, p2::T2)::Vector{qCompositeProduct} where {T1<:qComposite, T2<:qComposite}
    return [qCompositeProduct(p1.statespace, [p1, p2])]
end

function *(Q1::qExpr, Q2::qExpr)::qExpr
    if Q1.statespace != Q2.statespace
        error("Cannot multiply qExpr’s from different statespaces.")
    end
    new_terms::AbstractVector{qComposite} = []
    for t1 in Q1.terms
        for t2 in Q2.terms
            append!(new_terms, t1*t2)   
        end
    end
    return qExpr(Q1.statespace, new_terms)
    # return simplify(qExpr(Q1.statespace, new_terms))
end
function *(Q1::qExpr, Q2::qComposite)::qExpr
    Q_new = copy(Q1)
    terms = qComposite[]  # or Vector{qComposite}()
    for q in Q1.terms
        append!(terms, q * num)
    end
    Q_new.terms = terms
    return simplify(Q_new)
end
function *(Q1::qExpr, num::Number)::qExpr
    if num == 0
        return qExpr(Q1.statespace, qComposite[])  # or Vector{qComposite}() if preferred
    end
    terms = qComposite[]  # or Vector{qComposite}()
    for q in Q1.terms
        append!(terms, q * num)
    end
    return qExpr(Q1.statespace, terms)
end
function *(num::Number, Q1::qExpr)::qExpr
    return Q1*num
end

function *(Q1::qComposite, Q2::qExpr)::qComposite
    Q_new = copy(Q1)
    Q_new.expr = Q1.expr * Q2
    return Q_new
end
function *(Q1::qComposite, num::Number)::Vector{qComposite}
    Q_new = copy(Q1)
    Q_new.expr = Q1.expr * num
    return [Q_new]
end
function *(num::Number, Q2::qComposite)::Vector{qComposite}
    return [Q2 * num]
end

##### Exponentiation ###################################################
function ^(Q::Union{qExpr, qComposite}, n::Integer)
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
function Commutator(Q1::Union{qExpr, qComposite}, Q2::Union{qExpr, qComposite})
    # qExpr multiplication is already defined.
    return Q1 * Q2 - Q2 * Q1
end
function +(Q::qExpr, exprs::Vector{qExpr})::qExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return c + Q
end
function +(exprs::Vector{qExpr}, Q::qExpr)::qExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return c + Q 
end
function -(Q::qExpr, exprs::Vector{qExpr})::qExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return Q - c
end
function -(exprs::Vector{qExpr}, Q::qExpr)::qExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return c - Q
end
function *(Q::qExpr, exprs::Vector{qExpr})::qExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return Q * c
end
function *(exprs::Vector{qExpr}, Q::qExpr)::qExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return c * Q 
end



function Identity(qspace::StateSpace)
    # Set variable exponents to zero.
    var_exponents = zeros(Int, length(qspace.vars))
    # Build a vector of neutral operator indexes.
    neutral_ops = [s.op_set.neutral_element for s in qspace.subspaces for key in s.keys]
    # Create a single-term qExpr.
    return qExpr(qspace, qAtomProduct(qspace, qspace.fone, qAtom[qTerm(qspace.neutral_op)]))
end

"""
    Dag(qspace::StateSpace, t::qTerm) -> qTerm
    Dag(t::qExpr) -> qExpr
    Dag(t::qComposite) -> qComposite
    
Returns the Hermitian conjugate (dagger) of a qTerm, qExpr or qSum.
Overloads the adjoint function, which can be called via `t′`.
"""
function Dag(qspace::StateSpace, t::qTerm)::Vector{Tuple{qTerm, ComplexRational}}
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
        curr_coeff::ComplexRational = ComplexRational(1,0,1)
        for i in 1:length(combo)
            curr_coeff *= combo[i][1]
            push!(curr_inds, combo[i][2])
        end
        push!(terms, qTerm(curr_inds))
        push!(coeffs, curr_coeff)
    end
    return collect(zip(terms, coeffs))
end
function Dag(qspace::StateSpace, t::qAbstract)::Vector{Tuple{qAbstract, ComplexRational}} 
    # check if is hermitian
    if t.operator_type.hermitian
        return Tuple{qAbstract, ComplexRational}[(t, one(ComplexRational))]
    end 
    new_t = copy(t) 
    new_t.dag = !t.dag 
    return Tuple{qAbstract, ComplexRational}[(new_t, one(ComplexRational))]
end

function Dag(p::qAtomProduct)::Vector{qAtomProduct} 
    # reverse order elementwise Dag 
    terms::Vector{Vector{Tuple{qAtom, ComplexRational}}} = []
    for t in reverse(p.expr)
        push!(terms, Dag(p.statespace, t))
    end
    coeff_fun = p.coeff_fun
    new_atom_products::Vector{qAtomProduct} = []
    for combo in Iterators.product(terms...)
        curr_atoms = [a[1] for a in combo]
        curr_coeff = [a[2] for a in combo]
        coeff = reduce(*, curr_coeff)
        if !iszero(coeff)
            push!(new_atom_products, qAtomProduct(p.statespace, coeff_fun*coeff, curr_atoms))
        end
    end
    return new_atom_products 
end

function Dag(Q::qExpr)::qExpr
    q_comps::Vector{qComposite} = []
    for q in Q.terms
        append!(q_comps, Dag(q))
    end
    return qExpr(q_comps)
end
function Dag(t::qComposite)::Vector{qComposite}
    t_new = copy(t) 
    t_new.expr = Dag(t.expr)
    return [t_new]
end
function Dag(t::qMultiComposite)::Vector{qMultiComposite}
    t_new = copy(t) 
    dag_exprs::Vector{qExpr} = []
    for expr in reverse(t.expr)
        append!(dag_exprs, Dag(expr))
    end
    t_new.expr = dag_exprs 
    return [t_new]
end

adjoint(Q::qExpr) = Dag(Q)
adjoint(Q::qComposite) = Dag(Q)
