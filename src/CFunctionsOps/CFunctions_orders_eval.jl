export max_exponents, evaluate

import ..ConcreteIndexes
import ..CFunctions
using ..CFunctions: has_indexed_parameters, ParameterValues, set_param!, value, recompute_functions!

"""
    max_exponents(f::CFunction) -> Vector{Int}

Return the maximum absolute exponent for each parameter appearing in `f`. Composite
expressions propagate the per-variable maxima across their children; atoms yield the
absolute values of their sparse exponent vector.
"""
function max_exponents(f::CFunction)::Vector{Int}
    int_vec::Vector{Int} = zeros(Int, f.param_info.dims) 
    for leaf in leaf_iter(f) # iterate over the var_exponents::SparseVector{Int,Int}
        for i in leaf.var_exponents.nzind
            int_vec[i] = max(int_vec[i], abs(leaf.var_exponents[i]))
        end
    end
    return int_vec
end


# Rational exponent on a numeric base
@inline function _pow_r(y, q::Rational{Int})
    if denominator(q) == 1
        return y ^ Int(q)
    else
        return y ^ float(q)
    end
end

ctimes(c::ComplexRational, d::T) where T <: Number = (c.a + im*c.b) / c.c * d

function _evaluate_atom(a::CAtom, pv::ParameterValues, indexes::Union{Nothing,ConcreteIndexes})
    prod_val = one(Float64)
    first_term = true
    for idx in a.var_exponents.nzind
        val = value(pv, idx, indexes)
        exp = a.var_exponents[idx]
        term = val ^ exp
        if first_term
            prod_val = term
            first_term = false
        else
            prod_val *= term
        end
    end
    return ctimes(a.coeff, prod_val)
end

function _evaluate_atom(a::CAtomIndexed, pv::ParameterValues)
    prod_val = one(Float64)
    first_term = true
    for idx in a.var_exponents.nzind
        val = value(pv, idx, a.indexes)
        exp = a.var_exponents[idx]
        term = val ^ exp
        if first_term
            prod_val = term
            first_term = false
        else
            prod_val *= term
        end
    end
    return ctimes(a.coeff, prod_val)
end

function _evaluate(f::CFunction, pv::ParameterValues, indexes::Union{Nothing,ConcreteIndexes})
    if f isa CAtom
        if has_indexed_parameters(f) && indexes === nothing
            error("Cannot evaluate indexed atom without concrete indexes. Provide indexes via evaluate(...; indexes=...) or convert to CAtomIndexed.")
        end
        return _evaluate_atom(f::CAtom, pv, indexes)
    elseif f isa CAtomIndexed
        return _evaluate_atom(f::CAtomIndexed, pv)
    elseif f isa CAbstract
        error("CAbstract requires substitution before evaluation.")
    elseif f isa CIntegral
        error("CIntegral evaluation requires a dedicated backend.")
    elseif f isa CSum
        return sum(_evaluate(term, pv, indexes) for term in f.expr)
    elseif f isa CProd
        return ctimes(f.coeff, prod(_evaluate(term, pv, indexes) for term in f.expr))
    elseif f isa CRational
        return _evaluate(f.numer, pv, indexes) / _evaluate(f.denom, pv, indexes)
    elseif f isa CExp
        return ctimes(f.coeff, exp(_evaluate(f.expr, pv, indexes)))
    elseif f isa CLog
        return ctimes(f.coeff, log(_evaluate(f.expr, pv, indexes)))
    elseif f isa CPower
        return ctimes(f.coeff, _pow_r(_evaluate(f.expr, pv, indexes), f.exponent))
    elseif f isa CVector
        return [ctimes(f.coeff, _evaluate(term, pv, indexes)) for term in f.expr]
    elseif f isa CMatrix
        m, n = size(f.expr)
        return reshape([ctimes(f.coeff, _evaluate(term, pv, indexes)) for term in f.expr], m, n)
    elseif f isa CCustomType
        arg_values = [_evaluate(term, pv, indexes) for term in f.expr]
        val = evaluate(f.ctype_def.fun, arg_values, f.ctype_def.index_map)
        return ctimes(f.coeff, val)
    elseif f isa CCustomTypeIndexed
        arg_values = [_evaluate(term, pv, indexes) for term in f.expr]
        val = evaluate(f.ctype_def.fun, arg_values, f.ctype_def.index_map)
        return ctimes(f.coeff, val)
    else
        error("Evaluation not implemented for $(typeof(f)).")
    end
end

"""
    evaluate(f::CFunction, values::ParameterValues; indexes=nothing)

Evaluate a CFunction using the definitions in values. 
"""
function evaluate(f::CFunction, pv::ParameterValues; indexes::Union{Nothing,ConcreteIndexes}=nothing)
    recompute_functions!(pv)
    return _evaluate(f, pv, indexes)
end

### Needs to be rewritten