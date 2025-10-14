import ..ConcreteIndexes

"""
    Indexed(f::CFunction, indexes::ConcreteIndexes)

Return a copy of `f` where every coefficient leaf is converted into its indexed
counterpart using the supplied concrete ensemble indexes.
"""
function Indexed(f::CFunction, indexes::ConcreteIndexes)
    return _indexed(f, indexes)
end

function Indexed(::CFunction, indexes)
    throw(ArgumentError("Indexed expects indexes as a ConcreteIndexes object; construct one with ConcreteIndexes(param_info, ...)") )
end

struct CCustomTypeIndexed <: CFunction
    param_info::ParameterInfo
    coeff::ComplexRational
    expr::Vector{CFunction}
    ctype_def::CTypeDefinition
    indexes::ConcreteIndexes
    function CCustomTypeIndexed(param_info::ParameterInfo, coeff::ComplexRational,
                                expr::Vector{CFunction}, ctype_def::CTypeDefinition, indexes)
        concrete = indexes isa ConcreteIndexes ? indexes : ConcreteIndexes(param_info, indexes)
        return new(param_info, coeff, expr, ctype_def, concrete)
    end
end

CCustomTypeIndexed(c::CCustomType, indexes) =
    CCustomTypeIndexed(c.param_info, c.coeff, copy(c.expr), c.ctype_def, indexes)

var_exponents(c::CCustomTypeIndexed) = spzeros(Int, c.param_info.dims)
coeff(c::CCustomTypeIndexed) = [c.coeff]
length(::CCustomTypeIndexed) = 1
modify_expr(c::CCustomTypeIndexed, new_expr::Vector{CFunction}) =
    CCustomTypeIndexed(c.param_info, c.coeff, new_expr, c.ctype_def, c.indexes)
modify_coeff(c::CCustomTypeIndexed, coeff::ComplexRational) =
    CCustomTypeIndexed(c.param_info, coeff, c.expr, c.ctype_def, c.indexes)
function modify_indexes(c::CCustomTypeIndexed, indexes)
    concrete = indexes isa ConcreteIndexes ? indexes : ConcreteIndexes(c.param_info, indexes)
    return CCustomTypeIndexed(c.param_info, c.coeff, c.expr, c.ctype_def, concrete)
end

function _indexed(f::CFunction, ::ConcreteIndexes)
    return f
end

_indexed(a::CAtom, indexes::ConcreteIndexes) = CAtomIndexed(a, indexes)
_indexed(a::CAtomIndexed, indexes::ConcreteIndexes) = modify_indexes(a, indexes)

function _indexed(c::CCustomType, indexes::ConcreteIndexes)
    new_args = [_indexed(arg, indexes) for arg in c.expr]
    return CCustomTypeIndexed(c.param_info, c.coeff, new_args, c.ctype_def, indexes)
end

function _indexed(c::CCustomTypeIndexed, indexes::ConcreteIndexes)
    new_args = [_indexed(arg, indexes) for arg in c.expr]
    return CCustomTypeIndexed(c.param_info, c.coeff, new_args, c.ctype_def, indexes)
end

function _indexed(s::CSum, indexes::ConcreteIndexes)
    new_terms = [_indexed(term, indexes) for term in s.expr]
    return CSum(s.param_info, new_terms)
end

function _indexed(p::CProd, indexes::ConcreteIndexes)
    new_terms = [_indexed(term, indexes) for term in p.expr]
    return CProd(p.param_info, p.coeff, new_terms, Val(:nosimp))
end

function _indexed(r::CRational, indexes::ConcreteIndexes)
    new_num = _indexed(r.numer, indexes)
    new_den = _indexed(r.denom, indexes)
    return CRational(r.param_info, r.coeff, new_num, new_den)
end

function _indexed(exp::CExp, indexes::ConcreteIndexes)
    return CExp(exp.param_info, exp.coeff, _indexed(exp.expr, indexes))
end

function _indexed(logc::CLog, indexes::ConcreteIndexes)
    return CLog(logc.param_info, logc.coeff, _indexed(logc.expr, indexes))
end

function _indexed(cust::CCustomTypeIndexed, indexes::AbstractVector{<:AbstractVector{<:Integer}})
    return _indexed(cust, ConcreteIndexes(cust.param_info, indexes))
end

function _indexed(f::CFunction, indexes::AbstractVector{<:AbstractVector{<:Integer}})
    return _indexed(f, ConcreteIndexes(f.param_info, indexes))
end
