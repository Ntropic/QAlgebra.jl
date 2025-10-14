import ..ConcreteIndexes
import ..StringUtils: indexes2str
import ..Sampler: QDistribution, QEnsembleFunction

const _IndexLookupValue = Union{Number, AbstractVector, AbstractArray, Dict, Function}

function has_indexed_parameters(a::CAtom)::Bool
    params = a.param_info.params
    for ind in a.var_exponents.nzind
        params[ind].indexed_param && return true
    end
    return false
end

"""
    _validate_concrete_indexes(param_info, indexes)

Ensure the provided concrete index layout matches the ensemble sizes attached to
`param_info`. Accepts either a [`ConcreteIndexes`](@ref) instance or any
vector-of-vectors of integers.
"""
function _validate_concrete_indexes(param_info::ParameterInfo, indexes::ConcreteIndexes)
    expected = param_info.how_many_by_ensemble
    indexes.expected_lengths == expected ||
        error("Concrete indexes do not match ensemble sizes of the provided ParameterInfo.")
    return indexes
end

function _validate_concrete_indexes(param_info::ParameterInfo, indexes::AbstractVector{<:AbstractVector{<:Integer}})
    vectors = [Vector{Int}(idxs) for idxs in indexes]
    return _validate_concrete_indexes(param_info, ConcreteIndexes(param_info.how_many_by_ensemble, vectors))
end

ConcreteIndexes(param_info::ParameterInfo) =
    ConcreteIndexes(param_info.how_many_by_ensemble)
ConcreteIndexes(param_info::ParameterInfo, indexes::AbstractVector{<:AbstractVector{<:Integer}}) =
    _validate_concrete_indexes(param_info, indexes)
ConcreteIndexes(param_info::ParameterInfo, indexes::ConcreteIndexes) =
    _validate_concrete_indexes(param_info, indexes)

function _collect_index_values(indexes::ConcreteIndexes, ens_indexes::AbstractVector)
    values = Int[]
    for ens_idx in ens_indexes
        ensemble = getproperty(ens_idx, :outer)
        inner = getproperty(ens_idx, :inner)
        ensemble <= length(indexes.indexes) ||
            error("Concrete indexes missing ensemble $(ensemble).")
        ensemble_entries = indexes.indexes[ensemble]
        inner <= length(ensemble_entries) ||
            error("Concrete indexes missing entry $(inner) in ensemble $(ensemble).")
        idx = ensemble_entries[inner]
        idx > 0 || error("Concrete index for ensemble $(ensemble) inner $(inner) not set.")
        push!(values, idx)
    end
    return values
end

function _resolve_indexed_value(param_value, idxs::Vector{Int})
    isempty(idxs) && return param_value
    if param_value isa Function
        return param_value(idxs...)
    elseif param_value isa AbstractDict
        key = length(idxs) == 1 ? idxs[1] : Tuple(idxs)
        return param_value[key]
    elseif param_value isa AbstractArray
        return param_value[idxs...]
    elseif param_value isa AbstractVector
        length(idxs) == 1 ||
            error("Expected one index for vector parameter value, got $(length(idxs)).")
        return param_value[idxs[1]]
    elseif param_value isa Number
        return param_value
    elseif param_value === nothing
        error("Parameter value for indexed parameter is unspecified (nothing).")
    else
        error("Unsupported parameter storage type $(typeof(param_value)) for indexed evaluation.")
    end
end

"""
    CAtomIndexed(param_info, coeff, var_exponents, indexes)

Coefficient atom that keeps concrete ensemble indexes next to its sparse
exponent vector. `indexes` may be a [`ConcreteIndexes`](@ref) instance or any
vector of integer vectors aligned with the ensemble layout of `param_info`.
"""
struct CAtomIndexed <: CAtomic
    param_info::ParameterInfo
    coeff::ComplexRational
    var_exponents::SparseVector{Int,Int}
    indexes::ConcreteIndexes
    function CAtomIndexed(param_info::ParameterInfo, coeff::ComplexRational,
                          var_exponents::SparseVector{Int,Int}, indexes::ConcreteIndexes)
        return new(param_info, coeff, var_exponents, _validate_concrete_indexes(param_info, indexes))
    end
end

CAtomIndexed(param_info::ParameterInfo, var_exponents, indexes) =
    CAtomIndexed(CAtom(param_info, var_exponents), indexes)
CAtomIndexed(param_info::ParameterInfo, coeff::Number, var_exponents, indexes::AbstractVector{<:AbstractVector{<:Integer}}) =
    CAtomIndexed(CAtom(param_info, coeff, var_exponents), indexes)
CAtomIndexed(param_info::ParameterInfo, coeff::ComplexRational, var_exponents, indexes::AbstractVector{<:AbstractVector{<:Integer}}) =
    CAtomIndexed(CAtom(param_info, coeff, var_exponents), indexes)
CAtomIndexed(param_info::ParameterInfo, coeff::Complex, var_exponents, indexes::AbstractVector{<:AbstractVector{<:Integer}}) =
    CAtomIndexed(CAtom(param_info, coeff, var_exponents), indexes)
CAtomIndexed(param_info::ParameterInfo, coeff::ComplexRational, var_exponents::SparseVector{Int,Int}, indexes::AbstractVector{<:AbstractVector{<:Integer}}) =
    CAtomIndexed(param_info, coeff, var_exponents, _validate_concrete_indexes(param_info, indexes))

function CAtomIndexed(atom::CAtom, indexes)
    checked = _validate_concrete_indexes(atom.param_info, indexes)
    return CAtomIndexed(atom.param_info, atom.coeff, atom.var_exponents, checked)
end

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

var_exponents(a::CAtomIndexed) = a.var_exponents
coeff(a::CAtomIndexed)::Vector{ComplexRational} = [a.coeff]
length(::CAtomIndexed) = 1

function modify_exponents(a::CAtomIndexed, var_exponents)
    return CAtomIndexed(a.param_info, a.coeff, var_exponents, a.indexes)
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

modify_coeff(a::CAtomIndexed, coeff::ComplexRational) =
    CAtomIndexed(a.param_info, coeff, a.var_exponents, a.indexes)

function modify_indexes(a::CAtomIndexed, indexes)
    checked = _validate_concrete_indexes(a.param_info, indexes)
    return CAtomIndexed(a.param_info, a.coeff, a.var_exponents, checked)
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

function _indexed_parameter_name(atom::CAtomIndexed, idx::Int, do_latex::Bool)::String
    param = atom.param_info.params[idx]
    default = do_latex ? param.param_latex : param.param_str
    ens_indexes = param.ensemble_indexes
    isempty(ens_indexes) && return default
    indexes = atom.indexes
    values = Int[]
    for ens_idx in ens_indexes
        ensemble = ens_idx.outer
        inner = ens_idx.inner
        ensemble ≤ length(indexes.indexes) || error("Concrete indexes missing ensemble $(ensemble).")
        entries = indexes.indexes[ensemble]
        inner ≤ length(entries) || error("Concrete indexes missing entry $(inner) in ensemble $(ensemble).")
        push!(values, string(entries[inner]))
    end
    return default * indexes2str(values; do_latex=do_latex)
end

function _indexed_parameter_names(atom::CAtomIndexed, do_latex::Bool)::Vector{String}
    params = atom.param_info.params
    defaults = do_latex ? [p.param_latex for p in params] : [p.param_str for p in params]
    isempty(atom.indexes.indexes) && return defaults
    for idx in eachindex(defaults)
        isempty(params[idx].ensemble_indexes) && continue
        defaults[idx] = _indexed_parameter_label(atom.param_info, idx, atom.indexes, do_latex)
    end
    return defaults
end


function _indexed_parameter_label(param_info::ParameterInfo, param_index::Int, indexes::ConcreteIndexes, do_latex::Bool)
    param = param_info.params[param_index]
    base_str = do_latex ? param.param_latex : param.param_str
    ens_indexes = param.ensemble_indexes
    isempty(ens_indexes) && return base_str
    return _indexed_parameter_label_from_indexes(base_str, ens_indexes, indexes, do_latex)
end

function _indexed_parameter_label_from_indexes(base::String, ens_indexes::AbstractVector, indexes::ConcreteIndexes, do_latex::Bool)
    values = String[]
    for ens_idx in ens_indexes
        ensemble = getproperty(ens_idx, :outer)
        inner = getproperty(ens_idx, :inner)
        ensemble ≤ length(indexes.indexes) || error("Concrete indexes missing ensemble $(ensemble).")
        entries = indexes.indexes[ensemble]
        inner ≤ length(entries) || error("Concrete indexes missing entry $(inner) in ensemble $(ensemble).")
        push!(values, string(entries[inner]))
    end
    return base * indexes2str(values; do_latex=do_latex)
end

function _indexes_suffix(indexes::ConcreteIndexes, do_latex::Bool)
    flat = String[]
    for ensemble in indexes.indexes
        for idx in ensemble
            idx == 0 && continue
            push!(flat, string(idx))
        end
    end
    return indexes2str(flat; do_latex=do_latex)
end

function _indexed_parameter_value(param_info::ParameterInfo, param_index::Int, indexes::ConcreteIndexes, values::AbstractVector)
    length(values) == param_info.dims ||
        throw(DimensionMismatch("Expected a value vector of length $(param_info.dims), got $(length(values))."))
    param = param_info.params[param_index]
    ens_indexes = param.ensemble_indexes
    idxs = _collect_index_values(indexes, ens_indexes)
    raw_value = values[param_index]
    raw_value === nothing &&
        error("No numeric value supplied for parameter index $(param_index).")
    return _resolve_indexed_value(raw_value, idxs)
end
