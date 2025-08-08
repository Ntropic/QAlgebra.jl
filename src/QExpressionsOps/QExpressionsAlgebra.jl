export Dag, Commutator, is_numeric, same_statespace

""" 
    is_numeric(t::QTerm, qspace::StateSpace) -> Bool
    is_numeric(t::QAbstract, qspace::StateSpace) -> Bool
    is_numeric(p::QAtomProduct) -> Bool
    is_numeric(s::QSum) -> Bool
    is_numeric(expr::QExpr) -> Bool

Returns true if only the coefficient of the term(s) is non-zero.
"""
function is_numeric(op_indices::Vector{Is}, statespace::StateSpace)::Bool
    return statespace.neutral_op == op_indices
end
function is_numeric(e::QAtom, statespace::StateSpace)
    error("is_numeric (with given statespace) not implemented for QAbstract subtype $(typeof(e))")
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

function is_numeric(t::QTerm, statespace::StateSpace)::Bool
    return statespace.neutral_op == t.op_indices
end
function is_numeric(t::QAbstract, statespace::StateSpace)::Bool
    return false 
end
function is_numeric(s::QSum)::Bool
    return is_numeric(s.expr)
end
function is_numeric(expr::QExpr)::Bool
    expr_s = simplify(expr)
    terms = expr_s.terms

    if isempty(terms)
        return true  # No terms = numeric 0
    else
        return all(is_numeric(t) for t in terms)
    end
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
    if is_numeric(q) 
        c = q.coeff_fun
        if isnumeric(c)
            if isa(c, CAtom)
                return isone(c.coeff)
            else 
                return isone(simplify(c))
            end
        end
    end
    return false
end
function isone(q::QExpr)::Bool
    if isaQAtomProduct(q)
        isone(q.terms[1])
    else
        return false
    end
end

function where_neutral(q::QTerm, statespace::StateSpace)::Vector{Bool}
    return [op==neut for (op, neut) in zip(q.op_indices, statespace.neutral_op)]
end
function where_neutral(q::QAbstract, statespace::StateSpace)::Vector{Bool}
    return copy(q.operator_type.subspaces)
end
function where_acting(q::QTerm, statespace::StateSpace)::Vector{Bool}
    return [op!=neut for (op, neut) in zip(q.op_indices, statespace.neutral_op)]
end
function where_acting(q::QAbstract, statespace::StateSpace)::Vector{Bool}
    return copy(q.operator_type.non_subspaces)
end

function qAtom_commute(q1::QAtom, q2::QAtom, statespace::StateSpace)::Bool
    # check if all elements of where neutral are NAND
    act_q1 = where_acting(q1, statespace)
    act_q2 = where_acting(q2, statespace)
    return all([nand(a,b) for (a,b) in zip(act_q1, act_q2)])
end

function same_statespace(a::S, b::T)::Bool where {S<:QComposite, T<:QComposite}
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
    return a.coeff_fun == b.coeff_fun &&  all([ai == bi for (ai, bi) in zip(a.expr, b.expr)])
end

function ==(a::QExpr, b::QExpr)
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
    t_new = copy(t)
    t_new.expr = -t_new.expr  # Negate the expression inside the composite.
    return t_new 
end

#### Binary + ####################################################################
function +(Q1::QExpr, Q2::QExpr)::QExpr
    new_terms = vcat(Q1.terms, Q2.terms)
    # Optionally: group like terms here.
    return QExpr(Q1.statespace, new_terms)
end
function +(Q1::QExpr, Q2::T)::QExpr where T<:QComposite
    new_terms = copy(Q1.terms)
    push!(new_terms, Q2)
    return QExpr(Q1.statespace, new_terms)
end
function +(Q1::QExpr, N::Number)::QExpr 
    new_terms = copy(Q1.terms)
    push!(new_terms, QAtomProduct(Q1.statespace, Q1.statespace.fone*N, QTerm(Q1.statespace.neutral_op)))
    return QExpr(Q1.statespace, new_terms)
end
function +(N::Number, Q1::QExpr)::QExpr 
    new_terms = copy(Q1.terms)
    push!(new_terms, QAtomProduct(Q1.statespace, Q1.statespace.fone*N, QTerm(Q1.statespace.neutral_op)))
    return QExpr(Q1.statespace, new_terms)
end

#### Binary - ####################################################################
function -(Q1::QExpr, Q2::QExpr)::QExpr
    new_terms = vcat(Q1.terms, .-Q2.terms)
    # Optionally: group like terms here.
    return QExpr(Q1.statespace, new_terms)
end
function -(Q1::QExpr, Q2::T)::QExpr where T<:QComposite
    new_terms = copy(Q1.terms)
    push!(new_terms, .-Q2)
    return QExpr(Q1.statespace, new_terms)
end
function -(Q1::QExpr, N::Number)::QExpr 
    new_terms = copy(Q1.terms)
    push!(new_terms, QAtomProduct(Q1.statespace, -Q1.statespace.fone*N, QTerm(Q1.statespace.neutral_op)))
    return QExpr(Q1.statespace, new_terms)
end
function -(N::Number, Q1::QExpr)::QExpr 
    new_terms = copy((-Q1).terms)
    push!(new_terms, QAtomProduct(Q1.statespace, Q1.statespace.fone*N, QTerm(Q1.statespace.neutral_op)))
    return QExpr(Q1.statespace, new_terms)
end

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
    return QAtomProduct(Q1.statespace, Q1.coeff_fun*Q2.coeff_fun, newexpr)
end
# Multiplies two QTerm’s from the same statespace. Returns a vector of QTerm’s that are the result of this multiplication and corresponding ComplexRational coefficients. 
function multiply_qterm(t1::QTerm, t2::QTerm, statespace::StateSpace)::Tuple{Vector{QTerm}, Vector{ComplexRational}}
    # Multiply coefficients and add variable exponents.
    coeff_base = ComplexRational(1,0,1)

    # For the operator indexes, iterate over the subspaces.
    results_per_subspace::Vector{Vector{Tuple{ComplexRational,Is}}} = Vector{Tuple{ComplexRational,Is}}[]
    # We assume the length of op_indices equals the number of subspaces.
    i = 0
    for subspace in statespace.subspaces
        op_set = subspace.op_set
        for _ in eachindex(subspace.statespace_inds)
            i += 1
            op1 = t1.op_indices[i]
            op2 = t2.op_indices[i]
            # op_product should return a vector of tuples (factor, index).
            res = op_set.op_product(op1, op2)
            push!(results_per_subspace, res)
        end
    end

    new_terms = QTerm[]
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
        push!(new_terms, QTerm(new_term))
        push!(new_coeffs, new_factor)
    end
    return new_terms, new_coeffs
end            
function *(t1::QTerm, t2::QTerm, statespace::StateSpace)  # should never be called by the user!!
    return multiply_qterm(t1, t2, statespace)
end

#### Main Multiplication Functions ################################################
function *(p1::QAtomProduct, p2::QAtomProduct)::Vector{QAtomProduct}
    p = trivial_multiply(p1, p2)  # append the terms of p1 and p2.
    return simplifyqAtomProduct(p)    # simplify the product. 
end
function *(p1::QAtomProduct, num::Number)::Vector{QAtomProduct}
    return [QAtomProduct(p1.statespace, p1.coeff_fun*num, copy(p1.expr))]
end
function *(num::Number, p1::QAtomProduct)::Vector{QAtomProduct}
    return [QAtomProduct(p1.statespace, p1.coeff_fun*num, copy(p1.expr))]
end
function *(p1::T1, p2::T2)::Vector{QCompositeProduct} where {T1<:QComposite, T2<:QComposite}
    return [QCompositeProduct(vcat(p1, p2))]
end
function *(p1::QSum, p2::T2)::Vector{QSum} where T2<:QComposite
    new_expr::Vector{QComposite} = []
    for t in p1.expr
       append!(new_expr, t*p2)
    end
    return [QSum(p1.statespace, QExpr(new_expr), p1.indexes, p1.subsystem_index, p1.element_indexes, p1.neq)]
end
function *(p1::QCompositeProduct, p2::T2)::Vector{QCompositeProduct} where T2<:QComposite
    curr_comp = p1.expr
    last_prod = p1.expr[end] * p2
    return_vec::Vector{QCompositeProduct} = [] 
    for last in last_prod
        if isa(last, QCompositeProduct)
            push!(return_vec, QCompositeProduct(vcat(curr_comp, QExpr([p2]))))
        else 
            push!(return_vec, QCompositeProduct(vcat(curr_comp[1:end-1], QExpr([last]))))
        end
    end
    return return_vec
end
function *(p2::T2, p1::QCompositeProduct)::Vector{QCompositeProduct} where T2<:QComposite
    curr_comp = p1.expr
    first_prod = p2* p1.expr[1]
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
            append!(new_terms, t1*t2)   
        end
    end
    return QExpr(Q1.statespace, new_terms)
    # return simplify(QExpr(Q1.statespace, new_terms))
end
function *(Q1::QExpr, Q2::T)::QExpr where T<:QComposite
    terms = QComposite[]  # or Vector{QComposite}()
    for q in Q1.terms
        append!(terms, q * Q2)
    end
    return QExpr(Q1.statespace, terms)
end
function *(Q1::T, Q2::QExpr)::QExpr where T<:QComposite
    terms = QComposite[]  # or Vector{QComposite}()
    for q in Q2.terms
        append!(terms, Q1 * q)
    end
    return QExpr(Q1.statespace, terms)
end

function *(Q1::QExpr, num::Number)::QExpr
    if num == 0
        return QExpr(Q1.statespace, QComposite[])  # or Vector{QComposite}() if preferred
    end
    terms = QComposite[]  # or Vector{QComposite}()
    for q in Q1.terms
        append!(terms, q * num)
    end
    return QExpr(Q1.statespace, terms)
end
function *(num::Number, Q1::QExpr)::QExpr
    return Q1*num
end


function *(Q1::T, num::Number)::Vector{QComposite} where T<:QComposite
    Q_new = copy(Q1)
    Q_new.expr = Q1.expr * num
    return [Q_new]
end
function *(num::Number, Q2::T)::Vector{QComposite} where T<:QComposite
    return [Q2 * num]
end

##### Exponentiation ###################################################
function ^(Q::Union{QExpr, T}, n::Integer) where T<:QComposite
    if n < 0
        error("Negative exponent not defined for QExpr.")
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
    Commutator(Q1::QExpr, Q2::QExpr) -> QExpr
    Commutator(Q::QExpr, t::QSum) -> QExpr
    Commutator(t::QSum, Q::QExpr) -> QExpr

Computes the commutator [Q1, Q2] = Q1 * Q2 - Q2 * Q1.
Both QExpr's must share the same statespace.
"""
function Commutator(Q1::Union{QExpr, T}, Q2::Union{QExpr, QComposite}) where T<:QComposite
    # QExpr multiplication is already defined.
    return Q1 * Q2 - Q2 * Q1
end
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
    # Set variable exponents to zero.
    var_exponents = zeros(Int, length(qspace.vars))
    # Build a vector of neutral operator indexes.
    neutral_ops = [s.op_set.neutral_element for s in qspace.subspaces for key in s.keys]
    # Create a single-term QExpr.
    return QExpr(qspace, QAtomProduct(qspace, qspace.fone, QAtom[QTerm(qspace.neutral_op)]))
end

"""
    Dag(qspace::StateSpace, t::QTerm) -> QTerm
    Dag(t::QExpr) -> QExpr
    Dag(t::T) -> T where T <: QComposite
    
Returns the Hermitian conjugate (dagger) of a QTerm, QExpr or QSum.
Overloads the adjoint function, which can be called via `t′`.
"""
function Dag(qspace::StateSpace, t::QTerm)::Vector{Tuple{QTerm, ComplexRational}}
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
        curr_inds::Vector{Is} = []
        curr_coeff::ComplexRational = ComplexRational(1,0,1)
        for i in 1:length(combo)
            curr_coeff *= combo[i][1]
            push!(curr_inds, combo[i][2])
        end
        push!(terms, QTerm(curr_inds))
        push!(coeffs, curr_coeff)
    end
    return collect(zip(terms, coeffs))
end
function Dag(qspace::StateSpace, t::QAbstract)::Vector{Tuple{QAbstract, ComplexRational}} 
    # check if is hermitian
    if t.operator_type.hermitian
        return Tuple{QAbstract, ComplexRational}[(t, one(ComplexRational))]
    end 
    new_t = copy(t) 
    new_t.dag = !t.dag 
    return Tuple{QAbstract, ComplexRational}[(new_t, one(ComplexRational))]
end

function Dag(p::QAtomProduct)::Vector{QAtomProduct} 
    # reverse order elementwise Dag 
    terms::Vector{Vector{Tuple{QAtom, ComplexRational}}} = []
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
            push!(new_atom_products, QAtomProduct(p.statespace, coeff_fun*coeff, curr_atoms))
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
    t_new = copy(t) 
    dag_exprs::Vector{QExpr} = []
    for expr in reverse(t.expr)
        append!(dag_exprs, Dag(expr))
    end
    t_new.expr = dag_exprs 
    return [t_new]
end

adjoint(Q::QExpr) = Dag(Q)
adjoint(Q::T) where {T<:QComposite} = Dag(Q)
