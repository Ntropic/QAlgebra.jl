export max_exponents, evaluate

import ..ConcreteIndexes
import ..CFunctions
using ..CFunctions: has_indexed_parameters, ParameterValues, set_param!, value, recompute_functions!

function max_vec(a::Vector{Int}, b::Vector{Int})::Vector{Int}
    return max.(a, b)
end

max_exponents(a::CAtom) = Vector(abs.(a.var_exponents))
max_exponents(a::CAtomIndexed) = Vector(abs.(a.var_exponents))
max_exponents(a::CAbstract) = error("max_exponents cannot be computed for CAbstract.")
max_exponents(::CIntegral) = error("max_exponents is undefined for CIntegral.")

"""
    max_exponents(f::CFunction) -> Vector{Int}

Return the maximum absolute exponent for each parameter appearing in `f`. Composite
expressions propagate the per-variable maxima across their children; atoms yield the
absolute values of their sparse exponent vector.
"""
function max_exponents(f::CFunction)::Vector{Int}
    if f isa CAtom
        return max_exponents(f::CAtom)
    elseif f isa CAtomIndexed
        return max_exponents(f::CAtomIndexed)
    elseif f isa CSum
        dims = dims(f)
        m = zeros(Int, dims)
        for term in f.expr
            m = max_vec(m, max_exponents(term))
        end
        return m
    elseif f isa CProd
        dims = dims(f)
        m = zeros(Int, dims)
        for term in f.expr
            m = max_vec(m, max_exponents(term))
        end
        return m
    elseif f isa CPower
        return max_exponents(f.expr)
    elseif f isa CLog
        return max_exponents(f.expr)
    elseif f isa CExp
        return max_exponents(f.expr)
    elseif f isa CRational
        return max_vec(max_exponents(f.numer), max_exponents(f.denom))
    elseif f isa CVector
        dims = dims(f)
        m = zeros(Int, dims)
        for term in f.expr
            m = max_vec(m, max_exponents(term))
        end
        return m
    elseif f isa CMatrix
        dims = dims(f)
        m = zeros(Int, dims)
        for term in f.expr
            m = max_vec(m, max_exponents(term))
        end
        return m
    elseif f isa CCustomType
        dims = dims(f)
        m = zeros(Int, dims)
        for term in f.expr
            m = max_vec(m, max_exponents(term))
        end
        return m
    elseif f isa CCustomTypeIndexed
        dims = dims(f)
        m = zeros(Int, dims)
        for term in f.expr
            m = max_vec(m, max_exponents(term))
        end
        return m
    else
        error("max_exponents not implemented for $(typeof(f)).")
    end
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
    evaluate(f::CFunction, x::AbstractVector)

Numerically evaluate a coefficient expression. The first method consumes a
[`ParameterValues`](@ref) table (optionally restricting evaluation to concrete
index selections). The second method is a convenience wrapper that fills a fresh
`ParameterValues` instance from the numeric vector `x`.
"""
function evaluate(f::CFunction, pv::ParameterValues; indexes::Union{Nothing,ConcreteIndexes}=nothing)
    recompute_functions!(pv)
    return _evaluate(f, pv, indexes)
end

function evaluate(f::CFunction, x::AbstractVector{<:Number})
    pv = ParameterValues(f.param_info)
    n = min(length(x), length(pv.param_info.params_name))
    for idx in 1:n
        set_param!(pv, idx, x[idx])
    end
    return evaluate(f, pv)
end
