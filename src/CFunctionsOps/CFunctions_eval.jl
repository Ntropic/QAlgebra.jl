export evaluate

import ..ConcreteIndexes
import ..CFunctions
using ..CFunctions: ParameterValues #, value, recompute_functions!


# Rational exponent on a numeric base
@inline function _pow_r(y, q::Rational{Int})
    if denominator(q) == 1
        return y ^ Int(q)
    else
        return y ^ float(q)
    end
end

ctimes(c::ComplexRational, d::T) where T <: Number = (c.a + im*c.b) / c.c * d

@inline function _atom_product(pv::ParameterValues, exps, indexes::ConcreteIndexes)
    prod_val = 1.0
    first_term = true
    @inbounds for idx in exps.nzind
        val = value(pv, idx, indexes)
        exp = exps[idx]
        term = exp == 1 ? val : val ^ exp
        if first_term
            prod_val = term
            first_term = false
        else
            prod_val *= term
        end
    end
    return first_term ? 1.0 : prod_val
end

@inline function _evaluate_atom(a::CAtom, pv::ParameterValues, indexes::ConcreteIndexes)
    return ctimes(a.coeff, _atom_product(pv, a.var_exponents, indexes))
end

@inline function _evaluate_atom(a::CAtomIndexed, pv::ParameterValues)
    return ctimes(a.coeff, _atom_product(pv, a.var_exponents, a.indexes))
end

@inline function _evaluate(f::CAtom, pv::ParameterValues, indexes::ConcreteIndexes)
    return _evaluate_atom(f, pv, indexes)
end

@inline function _evaluate(f::CAtomIndexed, pv::ParameterValues, ::ConcreteIndexes)
    return _evaluate_atom(f, pv)
end

@inline function _evaluate(::CAbstract, ::ParameterValues, ::ConcreteIndexes)
    error("CAbstract requires substitution before evaluation.")
end

@inline function _evaluate(::CIntegral, ::ParameterValues, ::ConcreteIndexes)
    error("CIntegral evaluation requires a dedicated backend.")
end

@inline function _evaluate(f::CSum, pv::ParameterValues, indexes::ConcreteIndexes)
    return sum(_evaluate(term, pv, indexes) for term in f.expr)
end

@inline function _evaluate(f::CProd, pv::ParameterValues, indexes::ConcreteIndexes)
    return ctimes(f.coeff, prod(_evaluate(term, pv, indexes) for term in f.expr))
end

@inline function _evaluate(f::CRational, pv::ParameterValues, indexes::ConcreteIndexes)
    return _evaluate(f.numer, pv, indexes) / _evaluate(f.denom, pv, indexes)
end

@inline function _evaluate(f::CExp, pv::ParameterValues, indexes::ConcreteIndexes)
    return ctimes(f.coeff, exp(_evaluate(f.expr, pv, indexes)))
end

@inline function _evaluate(f::CLog, pv::ParameterValues, indexes::ConcreteIndexes)
    return ctimes(f.coeff, log(_evaluate(f.expr, pv, indexes)))
end

@inline function _evaluate(f::CPower, pv::ParameterValues, indexes::ConcreteIndexes)
    return ctimes(f.coeff, _pow_r(_evaluate(f.expr, pv, indexes), f.exponent))
end

@inline function _evaluate(f::CVector, pv::ParameterValues, indexes::ConcreteIndexes)
    entries = [ctimes(f.coeff, _evaluate(term, pv, indexes)) for term in f.expr]
    if f.row
        n = length(entries)
        return reshape(entries, 1, n)
    end
    return entries
end

@inline function _evaluate(f::CMatrix, pv::ParameterValues, indexes::ConcreteIndexes)
    m, n = size(f.expr)
    data = [ctimes(f.coeff, _evaluate(term, pv, indexes)) for term in f.expr]
    return reshape(data, m, n)
end

@inline function _evaluate(f::CCustomType, pv::ParameterValues, indexes::ConcreteIndexes)
    arg_values = [_evaluate(term, pv, indexes) for term in f.expr]
    val = evaluate(f.ctype_def.fun, arg_values, f.ctype_def.index_map)
    return ctimes(f.coeff, val)
end

@inline function _evaluate(f::CCustomTypeIndexed, pv::ParameterValues, indexes::ConcreteIndexes)
    arg_values = [_evaluate(term, pv, indexes) for term in f.expr]
    val = evaluate(f.ctype_def.fun, arg_values, f.ctype_def.index_map)
    return ctimes(f.coeff, val)
end

@inline function _evaluate(f::CFunction, ::ParameterValues, ::ConcreteIndexes)
    error("Evaluation not implemented for $(typeof(f)).")
end

"""
    evaluate(f::CFunction, pv::ParameterValues, indexes::ConcreteIndexes)

Evaluate `f` using parameter values stored in `pv` and the provided concrete
indexes. All atoms expect fully resolved indexes; callers are responsible for
ensuring index tuples point to valid samples.
"""
function evaluate(f::CFunction, pv::ParameterValues, indexes::ConcreteIndexes)
    recompute_functions!(pv)
    return _evaluate(f, pv, indexes)
end

function evaluate(f::CFunction, pv::ParameterValues; indexes::ConcreteIndexes)
    return evaluate(f, pv, indexes)
end
