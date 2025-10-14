export evaluate, abstract_evaluate

import ..ConcreteIndexes
import ..CFunctions
using ..CFunctions: ParameterValues, ParameterInfo #, value, recompute_functions!


# Rational exponent on a numeric base
@inline function _pow_r(y, q::Rational{Int})
    if denominator(q) == 1
        return y ^ Int(q)
    else
        return y ^ float(q)
    end
end

ctimes(c::ComplexRational, d::T) where T <: Number = (c.a + im*c.b) / c.c * d

@inline function _sample_indices(param, indexes::ConcreteIndexes)
    ens_indexes = param.ensemble_indexes
    count = length(ens_indexes)
    count == 0 && return ()
    return ntuple(pos -> begin
        ens_idx = ens_indexes[pos]
        return indexes.indexes[ens_idx.outer][ens_idx.inner]
    end, count)
end

@inline function _parameter_value(pinfo::ParameterInfo,
                                  pv::ParameterValues{SampleIndexMode},
                                  idx::Int,
                                  indexes::ConcreteIndexes)
    param = @inbounds pinfo.params[idx]
    group_idx = param.group_index
    time_idx = param.coords[1] - 1
    sample_indices = _sample_indices(param, indexes)
    return pv[group_idx, time_idx, sample_indices...]
end

@inline function _parameter_value(pinfo::ParameterInfo,
                                  pv::ParameterValues{AbstractIndexMode},
                                  idx::Int)
    param = @inbounds pinfo.params[idx]
    group_idx = param.group_index
    coords = param.coords
    time_idx = coords[1] - 1
    if length(coords) == 1
        return pv[group_idx, time_idx]
    else
        return pv[group_idx, time_idx, coords[2:end]...]
    end
end

@inline function _atom_product(pinfo::ParameterInfo,
                               pv::ParameterValues{SampleIndexMode},
                               exps,
                               indexes::ConcreteIndexes)
    prod_val = 1.0
    first_term = true
    @inbounds for idx in exps.nzind
        val = _parameter_value(pinfo, pv, idx, indexes)
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

@inline function _atom_product(pinfo::ParameterInfo,
                               pv::ParameterValues{AbstractIndexMode},
                               exps)
    prod_val = 1.0
    first_term = true
    @inbounds for idx in exps.nzind
        val = _parameter_value(pinfo, pv, idx)
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

# Sample-index evaluation -------------------------------------------------------

@inline function _evaluate_atom(a::CAtom,
                                pv::ParameterValues{SampleIndexMode},
                                indexes::ConcreteIndexes)
    return ctimes(a.coeff, _atom_product(a.param_info, pv, a.var_exponents, indexes))
end

@inline function _evaluate_atom(a::CAtomIndexed, pv::ParameterValues{SampleIndexMode})
    return ctimes(a.coeff, _atom_product(a.param_info, pv, a.var_exponents, a.indexes))
end

@inline function _atom_product(a::CAtomReferenced{Mode}) where {Mode<:ParameterValuesMode}
    prod_val = 1.0
    first_term = true
    @inbounds for (ptr, anchor, exp) in zip(a.ptrs, a.anchors, a.exponents)
        val = _factor_value(ptr, anchor)
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

@inline function _evaluate_atom(a::CAtomReferenced{Mode}) where {Mode<:ParameterValuesMode}
    return ctimes(a.coeff, _atom_product(a))
end

@inline function _evaluate(f::CAtom,
                           pv::ParameterValues{SampleIndexMode},
                           indexes::ConcreteIndexes)
    return _evaluate_atom(f, pv, indexes)
end

@inline function _evaluate(f::CAtomIndexed,
                           pv::ParameterValues{SampleIndexMode},
                           ::ConcreteIndexes)
    return _evaluate_atom(f, pv)
end

@inline function _evaluate(f::CEval,
                           ::ParameterValues{Mode},
                           ::ConcreteIndexes) where {Mode<:ParameterValuesMode}
    return f.value
end

@inline function _evaluate(f::CAtomReferenced{Mode},
                           pv::ParameterValues{Mode},
                           ::ConcreteIndexes) where {Mode<:ParameterValuesMode}
    pv === f.pv || error("CAtomReferenced is bound to a different ParameterValues instance.")
    return _evaluate_atom(f)
end

@inline function _evaluate(::CAbstract, ::ParameterValues, ::ConcreteIndexes)
    error("CAbstract requires substitution before evaluation.")
end

@inline function _evaluate(::CIntegral, ::ParameterValues, ::ConcreteIndexes)
    error("CIntegral evaluation requires a dedicated backend.")
end

@inline function _evaluate(f::CSum,
                           pv::ParameterValues{SampleIndexMode},
                           indexes::ConcreteIndexes)
    return sum(_evaluate(term, pv, indexes) for term in f.expr)
end

@inline function _evaluate(f::CProd,
                           pv::ParameterValues{SampleIndexMode},
                           indexes::ConcreteIndexes)
    return ctimes(f.coeff, prod(_evaluate(term, pv, indexes) for term in f.expr))
end

@inline function _evaluate(f::CRational,
                           pv::ParameterValues{SampleIndexMode},
                           indexes::ConcreteIndexes)
    return _evaluate(f.numer, pv, indexes) / _evaluate(f.denom, pv, indexes)
end

@inline function _evaluate(f::CExp,
                           pv::ParameterValues{SampleIndexMode},
                           indexes::ConcreteIndexes)
    return ctimes(f.coeff, exp(_evaluate(f.expr, pv, indexes)))
end

@inline function _evaluate(f::CLog,
                           pv::ParameterValues{SampleIndexMode},
                           indexes::ConcreteIndexes)
    return ctimes(f.coeff, log(_evaluate(f.expr, pv, indexes)))
end

@inline function _evaluate(f::CPower,
                           pv::ParameterValues{SampleIndexMode},
                           indexes::ConcreteIndexes)
    return ctimes(f.coeff, _pow_r(_evaluate(f.expr, pv, indexes), f.exponent))
end

@inline function _evaluate(f::CVector,
                           pv::ParameterValues{SampleIndexMode},
                           indexes::ConcreteIndexes)
    entries = [ctimes(f.coeff, _evaluate(term, pv, indexes)) for term in f.expr]
    if f.row
        n = length(entries)
        return reshape(entries, 1, n)
    end
    return entries
end

@inline function _evaluate(f::CMatrix,
                           pv::ParameterValues{SampleIndexMode},
                           indexes::ConcreteIndexes)
    m, n = size(f.expr)
    data = [ctimes(f.coeff, _evaluate(term, pv, indexes)) for term in f.expr]
    return reshape(data, m, n)
end

@inline function _evaluate(f::CCustomType,
                           pv::ParameterValues{SampleIndexMode},
                           indexes::ConcreteIndexes)
    arg_values = [_evaluate(term, pv, indexes) for term in f.expr]
    val = evaluate(f.ctype_def.fun, arg_values, f.ctype_def.index_map)
    return ctimes(f.coeff, val)
end

@inline function _evaluate(f::CCustomTypeIndexed,
                           pv::ParameterValues{SampleIndexMode},
                           indexes::ConcreteIndexes)
    arg_values = [_evaluate(term, pv, indexes) for term in f.expr]
    val = evaluate(f.ctype_def.fun, arg_values, f.ctype_def.index_map)
    return ctimes(f.coeff, val)
end

# Abstract-index evaluation -----------------------------------------------------

@inline function _evaluate_atom(a::CAtom, pv::ParameterValues{AbstractIndexMode})
    return ctimes(a.coeff, _atom_product(a.param_info, pv, a.var_exponents))
end

@inline function _abstract_evaluate(::CAbstract, ::ParameterValues{AbstractIndexMode})
    error("CAbstract requires substitution before evaluation.")
end

@inline function _abstract_evaluate(::CIntegral, ::ParameterValues{AbstractIndexMode})
    error("CIntegral evaluation requires a dedicated backend.")
end

@inline function _abstract_evaluate(f::CAtom, pv::ParameterValues{AbstractIndexMode})
    return _evaluate_atom(f, pv)
end

@inline function _abstract_evaluate(f::CEval, ::ParameterValues{AbstractIndexMode})
    return f.value
end

@inline function _abstract_evaluate(::CAtomIndexed, ::ParameterValues{AbstractIndexMode})
    error("CAtomIndexed cannot be evaluated with abstract index parameter values.")
end

@inline function _abstract_evaluate(f::CSum, pv::ParameterValues{AbstractIndexMode})
    return sum(_abstract_evaluate(term, pv) for term in f.expr)
end

@inline function _abstract_evaluate(f::CProd, pv::ParameterValues{AbstractIndexMode})
    return ctimes(f.coeff, prod(_abstract_evaluate(term, pv) for term in f.expr))
end

@inline function _abstract_evaluate(f::CRational, pv::ParameterValues{AbstractIndexMode})
    return _abstract_evaluate(f.numer, pv) / _abstract_evaluate(f.denom, pv)
end

@inline function _abstract_evaluate(f::CExp, pv::ParameterValues{AbstractIndexMode})
    return ctimes(f.coeff, exp(_abstract_evaluate(f.expr, pv)))
end

@inline function _abstract_evaluate(f::CLog, pv::ParameterValues{AbstractIndexMode})
    return ctimes(f.coeff, log(_abstract_evaluate(f.expr, pv)))
end

@inline function _abstract_evaluate(f::CPower, pv::ParameterValues{AbstractIndexMode})
    return ctimes(f.coeff, _pow_r(_abstract_evaluate(f.expr, pv), f.exponent))
end

@inline function _abstract_evaluate(f::CVector, pv::ParameterValues{AbstractIndexMode})
    entries = [ctimes(f.coeff, _abstract_evaluate(term, pv)) for term in f.expr]
    if f.row
        n = length(entries)
        return reshape(entries, 1, n)
    end
    return entries
end

@inline function _abstract_evaluate(f::CMatrix, pv::ParameterValues{AbstractIndexMode})
    m, n = size(f.expr)
    data = [ctimes(f.coeff, _abstract_evaluate(term, pv)) for term in f.expr]
    return reshape(data, m, n)
end

@inline function _abstract_evaluate(f::CCustomType, pv::ParameterValues{AbstractIndexMode})
    arg_values = [_abstract_evaluate(term, pv) for term in f.expr]
    val = evaluate(f.ctype_def.fun, arg_values, f.ctype_def.index_map)
    return ctimes(f.coeff, val)
end

@inline function _abstract_evaluate(f::CCustomTypeIndexed, pv::ParameterValues{AbstractIndexMode})
    arg_values = [_abstract_evaluate(term, pv) for term in f.expr]
    val = evaluate(f.ctype_def.fun, arg_values, f.ctype_def.index_map)
    return ctimes(f.coeff, val)
end

"""
    evaluate(f::CFunction, pv::ParameterValues{SampleIndexMode}, indexes::ConcreteIndexes)

Evaluate the concrete function `f` using sample-index parameter values stored in `pv`
and the concrete ensemble indexes provided via `indexes`. The caller is responsible
for supplying fully resolved indexes that match the ensemble layout.
"""
function evaluate(f::CFunction, pv::ParameterValues{SampleIndexMode}, indexes::ConcreteIndexes)
    return _evaluate(f, pv, indexes)
end

"""
    abstract_evaluate(f::CFunction, pv::ParameterValues{AbstractIndexMode})

Evaluate the concrete function `f` against abstract-index parameter values stored in `pv`.
Unlike `evaluate`, no concrete ensemble indexes are required because the parameter values
already encode the abstract index coordinates.
"""
function abstract_evaluate(f::CFunction, pv::ParameterValues{AbstractIndexMode})
    recompute_functions!(pv)
    return _abstract_evaluate(f, pv)
end
