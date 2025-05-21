export Dag, Commutator, isnumeric

""" 
    isnumeric(t::qTerm, statespace::StateSpace) -> Bool
    isnumeric(expr::qEQ) -> Bool
    isnumeric(s::qSum) -> Bool

Returns true if only the coefficient of the term(s) is non-zero.
"""
function isnumeric(t::qTerm, statespace::StateSpace)::Bool
    if statespace.neutral_op == t.op_indices && all(t.var_exponents .== 0)
        return true
    end
    return false
end
function isnumeric(expr::qEQ)::Bool
    expr_s = simplify(expr)
    terms = expr_s.terms

    if isempty(terms)
        return true  # No terms = numeric 0
    elseif length(terms) == 1
        t = terms[1]
        if t isa qTerm
            return isnumeric(t, expr_s.statespace)
        else
            return false
        end
    else
        return false
    end
end
function isnumeric(s::qSum)::Bool
    return false
end

import Base: ==
function isapprox_num(x, y; atol=1e-12)
    return isapprox(x, y, atol=atol)
end
function ==(a::qTerm, b::qTerm)
    return isapprox_num(a.coeff, b.coeff) &&
           a.var_exponents == b.var_exponents &&
           a.op_indices == b.op_indices
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
    if isnumeric(simple_expr)
        if length(simple_expr.terms) == 0
            return isapprox_num(0, n)
        else
            return isapprox_num(simple_expr.terms[1].coeff, n)
        end
    end
    return false
end
function ==(n::Number, expr::qEQ)::Bool
    return expr == n  # Symmetric
end

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
