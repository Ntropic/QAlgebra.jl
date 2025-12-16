import ..ConcreteIndexes

"""
    Indexed(f::CFunction, indices::ConcreteIndexes)

Return a copy of `f` where every coefficient leaf is converted into its indexed
counterpart using the supplied concrete ensemble indices.
"""
function Indexed(f::CFunction, indices::ConcreteIndexes)
    return _indexed(f, indices)
end

function Indexed(::CFunction, indices)
    throw(ArgumentError("Indexed expects indices as a ConcreteIndexes object; construct one with ConcreteIndexes(param_info, ...)") )
end

struct CCustomTypeIndexed <: CFunction
    param_info::ParameterInfo
    coeff::ComplexRational
    expr::Vector{CFunction}
    ctype_def::CTypeDefinition
    indices::ConcreteIndexes
    function CCustomTypeIndexed(param_info::ParameterInfo, coeff::ComplexRational,
                                expr::Vector{CFunction}, ctype_def::CTypeDefinition, indices)
        concrete = indices isa ConcreteIndexes ? indices : ConcreteIndexes(param_info, indices)
        return new(param_info, coeff, expr, ctype_def, concrete)
    end
end

CCustomTypeIndexed(c::CCustomType, indices) =
    CCustomTypeIndexed(c.param_info, c.coeff, copy(c.expr), c.ctype_def, indices)

var_exponents(c::CCustomTypeIndexed) = spzeros(Int, c.param_info.dims)
coeff(c::CCustomTypeIndexed) = [c.coeff]
length(::CCustomTypeIndexed) = 1
modify_expr(c::CCustomTypeIndexed, new_expr::Vector{CFunction}) =
    CCustomTypeIndexed(c.param_info, c.coeff, new_expr, c.ctype_def, c.indices)
modify_coeff(c::CCustomTypeIndexed, coeff::ComplexRational) =
    CCustomTypeIndexed(c.param_info, coeff, c.expr, c.ctype_def, c.indices)
function modify_indices(c::CCustomTypeIndexed, indices)
    concrete = indices isa ConcreteIndexes ? indices : ConcreteIndexes(c.param_info, indices)
    return CCustomTypeIndexed(c.param_info, c.coeff, c.expr, c.ctype_def, concrete)
end

function _indexed(f::CFunction, ::ConcreteIndexes)
    return f
end

_indexed(a::CAtom, indices::ConcreteIndexes) = CAtomIndexed(a, indices)
_indexed(a::CAtomIndexed, indices::ConcreteIndexes) = modify_indices(a, indices)

function _indexed(c::CCustomType, indices::ConcreteIndexes)
    new_args = [_indexed(arg, indices) for arg in c.expr]
    return CCustomTypeIndexed(c.param_info, c.coeff, new_args, c.ctype_def, indices)
end

function _indexed(c::CCustomTypeIndexed, indices::ConcreteIndexes)
    new_args = [_indexed(arg, indices) for arg in c.expr]
    return CCustomTypeIndexed(c.param_info, c.coeff, new_args, c.ctype_def, indices)
end

function _indexed(s::CSum, indices::ConcreteIndexes)
    new_terms = [_indexed(term, indices) for term in s.expr]
    return CSum(s.param_info, new_terms)
end

function _indexed(p::CProd, indices::ConcreteIndexes)
    new_terms = [_indexed(term, indices) for term in p.expr]
    return CProd(p.param_info, p.coeff, new_terms, Val(:nosimp))
end

function _indexed(r::CRational, indices::ConcreteIndexes)
    new_num = _indexed(r.numer, indices)
    new_den = _indexed(r.denom, indices)
    return CRational(r.param_info, r.coeff, new_num, new_den)
end

function _indexed(exp::CExp, indices::ConcreteIndexes)
    return CExp(exp.param_info, exp.coeff, _indexed(exp.expr, indices))
end

function _indexed(logc::CLog, indices::ConcreteIndexes)
    return CLog(logc.param_info, logc.coeff, _indexed(logc.expr, indices))
end

function _indexed(cust::CCustomTypeIndexed, indices::AbstractVector{<:AbstractVector{<:Integer}})
    return _indexed(cust, ConcreteIndexes(cust.param_info, indices))
end

function _indexed(f::CFunction, indices::AbstractVector{<:AbstractVector{<:Integer}})
    return _indexed(f, ConcreteIndexes(f.param_info, indices))
end
