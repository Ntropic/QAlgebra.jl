module CFunctions

using ..StringUtils
import ..SubSpaceIndex
import ..ConcreteIndexes
using ComplexRationals
using SparseArrays
using ..SparsePermutationTools: SparsePermutation
using ..QAlgebra: get_default, FLIP_IF_FIRST_TERM_NEGATIVE, DO_BRACED
using ..Sampler: QInterpolator, ContinuousSamples, integrate_node_funs
using ..ParameterGroups: ParameterGroup, ParameterGroupLike

export CFunction, CAbstractDefinition, CTypeDefinition, CIntegralDefinition, ParameterInfo
export define_cabstract, define_ctype, define_cintegral
export CAbstract, CIntegral, CCustomType, CCustomTypeIndexed, CAtom, CAtomIndexed, CAtomReferenced, CEval, CSum, CRational, CProd, CExp, CLog, CPower, CVector, CMatrix
export CMatrix, CVector, CPower
export coeff, var_exponents, unique_first_terms
export contains_non_simple_CFunction, Indexed
export list_cabstracts, list_ctypes, list_cintegrals
export where_acting, where_acting!, which_params_acting, which_params_acting!, param_index_tuples
export which_ensemble_acting, which_ensemble_acting!, substitute, separate_by_cond
export ParameterValues, update_t!, resolve_param!, resolve_ensemble_values!, refresh_ensemble_values!
export compute_integral_weights!

import Base: copy, exp, log, length, getindex, iterate, size
import ComplexRationals: isonelike
import ..QAlgebra: vecvec_or, vecvec_or!, sort_unique!, variants_C

const CR_ZERO = ComplexRational(0,0,1)
const CR_ONE  = ComplexRational(1,0,1)


"""
    CFunction

Abstract supertype for symbolic functions representing atoms (`CAtom`), sums (`CSum`), and rationals (`CRational`).
"""
abstract type CFunction end
abstract type CAtomic <: CFunction end
abstract type CComposite <: CFunction end 
abstract type CMultiComposite <: CFunction end 
abstract type CDef end

# =======================> Abstract CFun Definitions <===================================================================
abstract type AbstractCAbstract <: CAtomic end   # define here as a resesrvation, to concretely define later, for circular dependencies.
abstract type AbstractParameterInfo end
abstract type AbstractParameter end
"""
    CAbstractDef

Defines an abstract symbol, specifying its string and latex string, its index and referencing it to a CDefinitionsDB.
Contains:
- `symbol`     : unique tag
- `index`      : identifies which abstract definition in DB in is
- `plain`      : string to print for it 
- `latex`      : latexstring to print for i
This type is not to be used for manipulations of equations. Instead it is used to define new CTypes from existing ones. 
"""
struct CAbstractDefinition <: CDef
    symbol::Symbol
    name::String   # e.g. "A₁"
    latex::String  # e.g. "A_{1}"
    index::Int     # index in param_info
    sortkey::Int
    param_info::AbstractParameterInfo
end

"""
    CTypeDefinition

Row in the *functions* cluster (parametric defs like cos, sinh).
- `name`     : unique tag
- `plain`    : "cos"
- `latex`    : raw"\\cos"
- `sortkey`  : Int used by your sorter
- `fun` : (param_info, x)::CFunction → definitional expansion (outer coeff is not included)
"""
struct CTypeDefinition <: CDef
    name::Symbol
    type_symbols::NTuple{5,Symbol}
    plain::String
    latex::String
    index::Int    # index in param_info
    sortkey::Int
    fun::CFunction
    has_abstract::Bool
    abstract_parameters::Vector{AbstractCAbstract}
    index_map::Vector{Int} # maps the Cabstract.index to our abstractvector
    param_info::AbstractParameterInfo
end

struct ParameterIndexes # Helps find the indexes (ensemble and time indexes) associated with the parameters 
    labels::Vector{String}

    t_labels::Vector{String}
    t_labels_latex::Vector{String}
    label_parameter_indexes::Vector{Vector{Int}}
    label_parameter_t_indexes::Vector{Vector{Int}}

    all_indexes::Vector{Int}   # for has_indexes -> all parameters that would lead to an index being present 

    function ParameterIndexes(labels::Vector{String}, t_labels::Vector{String}, t_labels_latex::Vector{String}, label_parameter_indexes::Vector{Vector{Int}}, label_parameter_t_indexes::Vector{Vector{Int}})
        all_indexes = sort_unique!(vcat(vcat(label_parameter_indexes...),vcat(label_parameter_t_indexes...)))
        new(labels, t_labels, t_labels_latex, label_parameter_indexes, label_parameter_t_indexes, all_indexes)
    end
end
"""
    CIntegralDefinition

Container describing a coefficient integral definition stored within a
[`ParameterInfo`](@ref). It records the defining integrand `expr` and the
subsystem indexes integrated over.
"""
struct CIntegralDefinition{N} <: CDef
    index::Int
    sortkey::Int
    expr::CFunction
    indexes::Vector{Vector{SubSpaceIndex}}
    param_info::AbstractParameterInfo
    interpolator::QInterpolator
    axis_lengths::Vector{Int}
    pdfs::Vector{Function}
    assignments::Vector{Vector{Int}}  # cached parameter index lists for assigning interpolated values
    values::Array{ComplexF64,N}
    function CIntegralDefinition(index::Int,
                                 sortkey::Int,
                                 expr::CFunction,
                                 indexes::Vector{Vector{SubSpaceIndex}},
                                 qspace)
        helpers = getfield(parentmodule(@__MODULE__), :QExpressions)
        param_info = qspace.param_info
        interp, axis_lengths, pdfs, dim_info = helpers._build_integral_interpolator(qspace, indexes)
        assignments = helpers._build_integral_assignments(param_info, dim_info)
        N = length(axis_lengths)
        N > 0 || error("CIntegralDefinition requires at least one integration dimension.")
        vals = Array{ComplexF64}(undef, axis_lengths...)
        return new{N}(index, sortkey, expr, indexes, param_info, interp, axis_lengths, pdfs, assignments, vals)
    end
end
const AnyCIntegralDefinition = CIntegralDefinition{N} where N

"""
    ParameterInfo

Holds both clusters and the dimension of the polynomial variable space.
"""
struct ParameterInfo <: AbstractParameterInfo
    dims::Int
    outer_labels_symbols::Vector{Symbol}
    inner_labels_symbols_flat::Vector{Symbol}

    outer_labels::Vector{String}
    outer_labels_str::Vector{String}
    outer_labels_latex::Vector{String}

    subspace_index_maps::Vector{Array{SparsePermutation,2}}
    t_index_transform::Array{SparsePermutation,2}
    indexes_by_t_index::Vector{Vector{Int}}
    indexes_of_t::Vector{Int}

    how_many_by_ensemble::Vector{Int}

    subspace_info::Any
    param_indexes::ParameterIndexes
    param_groups::Vector{ParameterGroupLike}
    params::Vector{AbstractParameter}
    abstract_definitions::Vector{CAbstractDefinition}
    custom_ctype::Vector{CTypeDefinition}
    integral_definitions::Vector{AnyCIntegralDefinition}

    function ParameterInfo(
        outer_labels_symbols::Vector{Symbol}, inner_labels_symbols_flat::Vector{Symbol}, outer_labels::Vector{String},
        outer_labels_str::Vector{String}, outer_labels_latex::Vector{String},
        subspace_index_maps::Vector{Array{SparsePermutation,2}}, t_index_transform::Array{SparsePermutation,2}, indexes_by_t_index::Vector{Vector{Int}},
        indexes_of_t::Vector{Int}, how_many_by_ensemble::Vector{Int},
        subspace_info::Any, param_indexes::ParameterIndexes, param_groups::Vector{ParameterGroupLike}, params::AbstractVector{<:AbstractParameter})
        dims = length(inner_labels_symbols_flat)
        stored_params = AbstractParameter[params...]
        new(dims, outer_labels_symbols, inner_labels_symbols_flat, outer_labels,
            outer_labels_str, outer_labels_latex,
            subspace_index_maps, t_index_transform,
            indexes_by_t_index, indexes_of_t, how_many_by_ensemble,
            subspace_info, param_indexes, param_groups, stored_params, CAbstractDefinition[], CTypeDefinition[], AnyCIntegralDefinition[])
    end
end
function Base.show(io::IO, info::ParameterInfo)
    params = info.params
    labels = String[]
    for group_idx in 1:length(info.param_groups)
        idx = findfirst(p -> p.group_index == group_idx, params)
        idx === nothing && continue
        param = params[idx]
        push!(labels, param.param_str)
    end
    print(io, "ParameterInfo([", join(labels, ","), "])")
end

include("CFunctionsOps/ParameterDicts.jl")


"""
    param_index_tuples(param_info::ParameterInfo, param_index::Int)

Return the cached [`EnsembleIndex`] entries describing where parameter
`param_index` acts. Non-indexed parameters yield an empty vector.
"""
function param_index_tuples(param_info::ParameterInfo, param_index::Int)
    1 ≤ param_index ≤ length(param_info.params) ||
        error("Parameter index $(param_index) out of bounds.")
    return param_info.params[param_index].ensemble_indexes
end

"""
    coeff(f::CFunction) -> Vector{ComplexRational}

Return the scalar coefficients present in `f`. For atomic objects this is the
single leading coefficient; for structured expressions the result collects the
scalars contributed by each branch.
"""
function coeff end

"""
    var_exponents(f::CFunction) -> Vector{Int}

Return the polynomial exponents associated with each variable in `f`. Composite
objects delegate to their children, while purely numeric constructs return a
zero vector.
"""
function var_exponents end

include("CFunctionsDefinitions/CIntegral.jl")
include("CFunctionsDefinitions/CAbstract.jl")
include("CFunctionsDefinitions/CCustom.jl")

# ===================> MAIN TYPES <==========================================================================================
function modify_expr(f::CFunction, new_expr::Vector{CFunction})
    error("modify_expr not implemented for type $(typeof(f)).")
end

@inline function _sparse_exponents(param_info::ParameterInfo, exps)::SparseVector{Int}
    exps isa AbstractVector || return _sparse_exponents(param_info, collect(exps))
    length(exps) == param_info.dims || throw(DimensionMismatch("expected $(param_info.dims) exponents, got $(length(exps))"))
    return SparseVector{Int}(exps)
end
"""
    CAtom(param_info::ParameterInfo, var_exponents)
    CAtom(param_info::ParameterInfo, coeff::Number, var_exponents)

Single polynomial atom with a complex-rational coefficient and integer exponents per
variable. Sparse storage keeps zero exponents implicit.

Notes:
- Pure exponent constructor assumes coefficient 1.
- Numeric coefficients are promoted to `ComplexRational`.
- `var_exponents[j]` stores the power of variable _j_ defined in `param_info`.
"""
struct CAtom <: CAtomic
    param_info::ParameterInfo
    coeff::ComplexRational
    var_exponents::SparseVector{Int,Int}
    function CAtom(param_info::ParameterInfo, var_exponents)
        c = ComplexRational(1, 0, 1)
        return new(param_info, c, _sparse_exponents(param_info, var_exponents))
    end
    function CAtom(param_info::ParameterInfo, coeff::Int, var_exponents)
        c = ComplexRational(coeff, 0, 1)
        return new(param_info, c, _sparse_exponents(param_info, var_exponents))
    end
    function CAtom(param_info::ParameterInfo, coeff::Rational, var_exponents)
        c = ComplexRational(numerator(coeff), 0, denominator(coeff))
        return new(param_info, c, _sparse_exponents(param_info, var_exponents))
    end
    function CAtom(param_info::ParameterInfo, coeff::Complex, var_exponents)
        c = crationalize(coeff)
        return new(param_info, c, _sparse_exponents(param_info, var_exponents))
    end
    function CAtom(param_info::ParameterInfo, coeff::ComplexRational, var_exponents)
        return new(param_info, coeff, _sparse_exponents(param_info, var_exponents))
    end
    function CAtom(param_info::ParameterInfo, coeff::Number, var_exponents)
        c = crationalize(coeff + 0im)
        return new(param_info, c, _sparse_exponents(param_info, var_exponents))
    end
end
@inline function _sparse_exponents(param_info::ParameterInfo, exps::Tuple)
    return _sparse_exponents(param_info, collect(exps))
end
@inline function zero_catom(param_info::ParameterInfo)
    return CAtom(param_info, CR_ZERO, spzeros(Int, param_info.dims))
end

coeff(a::CAtom)::Vector{ComplexRational} = [a.coeff]
modify_coeff(a::CAtom, coeff::ComplexRational)::CAtom = CAtom(a.param_info, coeff, a.var_exponents)
var_exponents(a::CAtom) = a.var_exponents
length(a::CAtom) = 1
function repartition(f::CAtom, var_tuples::Vector{Tuple{Int, Int}})::CAtom 
    curr_var_exponents = copy(f.var_exponents)
    @inbounds for (i, tar) in var_tuples
        curr_var_exponents[tar] += curr_var_exponents[i]
        curr_var_exponents[i] = 0 
    end 
    CAtom(f.param_info, f.coeff, curr_var_exponents)
end

"""
    CSum(expr::AbstractVector{<:CFunction})

Constructs a sum of `CFunction` expr.
- Flattens any nested `CSum` automatically.
- Variadic form `CSum(a, b, c)` is provided for convenience.
"""
struct CSum <: CMultiComposite
    param_info::ParameterInfo
    expr::Vector{CFunction}
end
function _CSum(param_info::ParameterInfo, ts::AbstractVector{<:CFunction}) 
    if length(ts) == 1
        return ts[1]
    end
    for t in ts
        if t isa CSum
            error("Shouldn't have a CSum in a CSum!")  # remove this loop later on 
        end
    end
    return simplify_CSum(param_info, collect(ts))
end
function _CSum(param_info::ParameterInfo, ts::AbstractVector{<:CFunction}, ::Val{:nosimp})
    return CSum(param_info, ts)
end
function modify_expr(f::CSum, new_expr::Vector{CFunction})
    return CSum(f.param_info, new_expr)
end
function modify_coeff(f::CSum, coeff::ComplexRational)
    iszero(coeff) && return zero_catom(f.param_info)
    coeff == CR_ONE && return f
    return f * coeff
end
coeff(x::CSum) = [ComplexRational(1,0,1)] #error("Sums don't have a coeff, you likely have a sum in a sum, this shouldn't happen. Please inform the developers. ")
length(q::CSum) = length(q.expr)
repartition(f::CSum, var_tuples::Vector{Tuple{Int, Int}}) = _CSum(f.param_info, repartition.(f.expr, Ref(var_tuples)) )
var_exponents(a::CSum) = min.(var_exponents.(a.expr)...)

"""
    CProd(expr::AbstractVector{<:CFunction})
    CProd(coeff::ComplexRational, expr::AbstractVector{<:CFunction})

Product of coefficient expressions. The constructor performs light
simplification (flattening nested products and combining scalars) unless the
`Val(:nosimp)` variant is used.
"""
struct CProd <: CMultiComposite
    param_info::ParameterInfo
    coeff::ComplexRational
    expr::Vector{CFunction}
    function CProd(param_info::ParameterInfo, coeff::ComplexRational, expr::AbstractVector{<:CFunction})
        if length(expr) == 1
            return expr[1] * coeff
        end
        return simplify_CProd(param_info, coeff, collect(expr))
    end
    function CProd(param_info::ParameterInfo, coeff::ComplexRational, expr::AbstractVector{<:CFunction}, ::Val{:nosimp})
        return new(param_info, coeff, copy(expr))
    end
end
function CProd(param_info::ParameterInfo, expr::AbstractVector{<:CFunction})
    CProd(param_info, ComplexRational(1, 0, 1), collect(expr))
end
function modify_expr(f::CProd, new_expr::Vector{CFunction})
    return CProd(f.param_info, f.coeff, new_expr, Val(:nosimp))
end
function modify_coeff(f::CProd, coeff::ComplexRational)::CFunction
    iszero(coeff) && return zero_catom(f.param_info)
    return CProd(f.param_info, coeff, f.expr, Val(:nosimp))
end
coeff(x::CProd) = [x.coeff]
length(q::CProd) = max(length.(q.expr)...)
repartition(f::CProd, var_tuples::Vector{Tuple{Int, Int}})= CProd(f.param_info, f.coeff, repartition.(f.expr, Ref(var_tuples)) )
function var_exponents(a::CProd) 
    if length(a.expr) > 0 
        return var_exponents(a.expr[1])
    else 
        return spzeros(Int, a.param_info.dims)
    end
end

"""
    CRational(numer::CSum, denom::CSum)

Represents a rational function with numerator `numer` and denominator `denom`, both sums of `CFunction` expr.
"""
struct CRational <: CFunction   # special case
    param_info::ParameterInfo
    numer::CFunction
    denom::CFunction
    function CRational(param_info::ParameterInfo, numer::T, denom::S) where {T <: CFunction, S <: CFunction}  
        return simplify_CRational(param_info, numer, denom)
    end
    function CRational(param_info::ParameterInfo, numer::T, denom::S,  ::Val{:nosimp}) where {T <: CFunction, S <: CFunction}  
        return new(param_info, numer, denom)
    end
end
coeff(x::CRational) = coeff(x.numer) #/coeff(x.denom)
length(q::CRational) = max(length(q.numer), length(q.denom))
repartition(q::CRational, var_tuples::Vector{Tuple{Int, Int}}) = CRational(q.param_info, repartition(q.numer, var_tuples), repartition(q.denom, var_tuples))
function modify_coeff(r::CRational, coeff::ComplexRational)::CFunction
    iszero(coeff) && return zero_catom(r.param_info)
    numer = modify_coeff(r.numer, coeff)
    return CRational(r.param_info, numer, r.denom, Val(:nosimp))
end
var_exponents(a::CRational) = var_exponents(a.numer)


"""
    CExp(expr::CFunction)
    CExp(coeff::ComplexRational, expr::CFunction)

Symbolic exponential `coeff * exp(expr)` used for closed-form coefficient
expressions. Construction simplifies simple logarithmic inverses automatically.
"""
struct CExp <: CComposite
    param_info::ParameterInfo
    coeff::ComplexRational
    expr::CFunction
    function CExp(param_info::ParameterInfo, coeff::ComplexRational, expr::T,  ::Val{:nosimp}) where T <: CFunction
        new(param_info, coeff, expr)
    end
    function CExp(param_info::ParameterInfo, coeff::ComplexRational, expr::T) where T <: CFunction
        return simplify_CExp(param_info, coeff, expr)
    end
end
function CExp(param_info::ParameterInfo, expr::CFunction)
    CExp(param_info, ComplexRational(1,0,1), expr)
end
function exp(param_info::ParameterInfo, expr::CFunction)
    return CExp(param_info, expr) 
end
function exp(expr::CFunction)
    return CExp(expr.param_info, expr) 
end
function modify_expr(f::CExp, new_expr::Vector{CFunction})
    @assert length(new_expr) == 1
    return CExp(f.param_info, f.coeff, new_expr[1], Val(:nosimp))
end
function modify_coeff(f::CExp, coeff::ComplexRational)::CFunction
    iszero(coeff) && return zero_catom(f.param_info)
    return CExp(f.param_info, coeff, f.expr, Val(:nosimp))
end
coeff(x::CExp) = [x.coeff]
length(q::CExp) = 1
repartition(q::CExp, var_tuples::Vector{Tuple{Int, Int}}) = CExp(q.param_info, q.coeff, repartition(q.expr, var_tuples))
var_exponents(a::CExp) = spzeros(Int, a.param_info.dims)


"""
    CLog(expr::CFunction)
    CLog(coeff::ComplexRational, expr::CFunction)

Symbolic logarithm `coeff * log(expr)` with light simplification (e.g.
`log(exp(x)) → x`).
"""
struct CLog <: CComposite
    param_info::ParameterInfo
    coeff::ComplexRational
    expr::CFunction
    function CLog(param_info::ParameterInfo, coeff::ComplexRational, expr::T,  ::Val{:nosimp}) where T <: CFunction
        new(param_info, coeff, expr)
    end
    function CLog(param_info::ParameterInfo, coeff::ComplexRational, expr::T)  where T <: CFunction
        return simplify_CLog(param_info, coeff, expr)
    end
end
function CLog(param_info::ParameterInfo, expr::CFunction)
    CLog(param_info, ComplexRational(1,0,1), expr)
end
function log(param_info::ParameterInfo, expr::CFunction)
    return CLog(param_info, expr) 
end
function log(expr::CFunction)
    return CLog(expr.param_info, expr) 
end
function modify_expr(f::CLog, new_expr::Vector{CFunction})
    @assert length(new_expr) == 1
    return CLog(f.param_info, f.coeff, new_expr[1], Val(:nosimp))
end
function modify_coeff(f::CLog, coeff::ComplexRational)::CFunction
    iszero(coeff) && return zero_catom(f.param_info)
    return CLog(f.param_info, coeff, f.expr, Val(:nosimp))
end
coeff(x::CLog) = [x.coeff] 
length(q::CLog) = 1
repartition(q::CLog, var_tuples::Vector{Tuple{Int, Int}}) = CLog(q.param_info, q.coeff, repartition(q.expr, var_tuples))
var_exponents(a::CLog) = spzeros(Int, a.param_info.dims)

"""
    CPower(coeff::ComplexRational, x::CFunction, exponent::Rational{Int})
    CPower(x::CFunction, exponent::Rational{Int})
    CPower(x::CFunction, exponent::Int)

Symbolic power with a **rational** exponent: `coeff * x^(p//q)`.
Use `x ^ (p//q)` or `sqrt(x)` (which maps to `x^(1//2)`).
"""
struct CPower <: CComposite
    param_info::ParameterInfo
    coeff::ComplexRational
    expr::CFunction
    exponent::Rational{Int}
    # inner :nosimp constructor that *only* wraps
    function CPower(param_info::ParameterInfo, coeff::ComplexRational, expr::T, exponent::Rational{Int}, ::Val{:nosimp}) where {T<:CFunction}
        new(param_info, coeff, expr, exponent)
    end    
end

# thin outer constructors that delegate to simplify
function CPower(param_info::ParameterInfo, coeff::ComplexRational, expr::CFunction, exponent::Rational{Int})
    simplify_CPower(param_info, coeff, expr, exponent)
end
CPower(param_info::ParameterInfo, expr::CFunction, n::Int)       = CPower(param_info, ComplexRational(1,0,1), expr, n//1)
CPower(param_info::ParameterInfo, expr::CFunction, q::Rational{Int}) = CPower(param_info, ComplexRational(1,0,1), expr, q)
function modify_expr(f::CPower, new_expr::Vector{CFunction})
    @assert length(new_expr) == 1
    return CPower(f.param_info, f.coeff, new_expr[1], f.exponent, Val(:nosimp))
end
function modify_coeff(f::CPower, coeff::ComplexRational)::CFunction
    iszero(coeff) && return zero_catom(f.param_info)
    return CPower(f.param_info, coeff, f.expr, f.exponent, Val(:nosimp))
end
coeff(p::CPower) = [p.coeff]
length(::CPower) = 1
repartition(p::CPower, var_tuples::Vector{Tuple{Int,Int}}) = CPower(p.param_info, p.coeff, repartition(p.expr, var_tuples), p.exponent)
var_exponents(a::CPower) = spzeros(Int, a.param_info.dims)


"""
    CVector(entries::AbstractVector{<:CFunction}; row::Bool=false)
    CVector(coeff::ComplexRational, entries::AbstractVector{<:CFunction}; row::Bool=false)

An oriented vector of `CFunction`s.
- `row=false` ⇒ n×1 (column, Julia's default)
- `row=true`  ⇒ 1×n (row)
"""
struct CVector <: CFunction
    param_info::ParameterInfo
    coeff::ComplexRational
    expr::Vector{CFunction}
    row::Bool  # false => n×1 (column), true => 1×n (row)
end
CVector(param_info::ParameterInfo, expr::AbstractVector{<:CFunction}; row::Bool=false) = CVector(param_info, ComplexRational(1,0,1), collect(expr), row)
CVector(param_info::ParameterInfo, coeff::ComplexRational, expr::AbstractVector{<:CFunction}; row::Bool=false) = CVector(param_info, coeff, collect(expr), row)
modify_exprs(f::CVector, new_expr::Vector{CFunction}) = CVector(f.param_info, f.coeff, new_expr; row=f.row)
function modify_coeff(v::CVector, coeff::ComplexRational)::CVector
    return CVector(v.param_info, coeff, v.expr; row=v.row)
end
coeff(v::CVector) = isempty(v.expr) ? ComplexRational[] : vcat(coeff.(v.expr)...)
length(v::CVector) = length(v.expr)
size(v::CVector) = v.row ? (1, length(v.expr)) : (length(v.expr), 1)
getindex(v::CVector, i::Int) = v.expr[i]
iterate(v::CVector, st::Int=1) = st > length(v.expr) ? nothing : (v.expr[st], st+1)
repartition(v::CVector, var_tuples::Vector{Tuple{Int,Int}}) = CVector(v.param_info, v.coeff, repartition.(v.expr, Ref(var_tuples)); row=v.row)
var_exponents(a::CVector) = spzeros(Int, a.param_info.dims)


"""
    CMatrix(entries::AbstractMatrix{<:CFunction})
    CMatrix(coeff::ComplexRational, entries::AbstractMatrix{<:CFunction})

A matrix of `CFunction`s. 
"""
struct CMatrix <: CFunction
    param_info::ParameterInfo
    coeff::ComplexRational
    expr::Matrix{CFunction}
end
CMatrix(param_info::ParameterInfo, expr::AbstractMatrix{<:CFunction}) = CMatrix(param_info, ComplexRational(1,0,1), Matrix{CFunction}(expr))
CMatrix(param_info::ParameterInfo, coeff::ComplexRational, expr::AbstractMatrix{<:CFunction}) = CMatrix(param_info, coeff, Matrix{CFunction}(expr))
modify_exprs(f::CMatrix, new_expr::Matrix{CFunction}) = CMatrix(f.param_info, f.coeff, new_expr)
function modify_coeff(M::CMatrix, coeff::ComplexRational)::CMatrix
    return CMatrix(M.param_info, coeff, M.expr)
end
coeff(M::CMatrix) = [M.coeff]
length(M::CMatrix) = length(M.expr)         # number of elements (m*n)
size(M::CMatrix) = size(M.expr)
getindex(M::CMatrix, i::Int, j::Int) = M.expr[i, j]
repartition(M::CMatrix, var_tuples::Vector{Tuple{Int,Int}}) = CMatrix(M.param_info, M.coeff, reshape(repartition.(M.expr[:], Ref(var_tuples)), size(M.expr)))
var_exponents(a::CMatrix) = spzeros(Int, a.param_info.dims)

#### Some basic functions ##############################################################################################



import Base: length, getindex, iterate, deleteat!, reverse
length(p::CFunction)::Int = 1

getindex(p::CSum, i::Int) = p.expr[i]
iterate(p::CSum, state=1) = state > length(p.expr) ? nothing : (p.expr[state], state + 1)
deleteat!(p::CSum, i::Int) = _CSum(deleteat!(p.expr, i))
reverse(q::CSum) = CSum(reverse(q.expr))

"""
    contains_non_simple_CFunction(c::CFunction) -> Bool 

Return `true` if `c` contains non-trivial constructs such as `CExp`, `CLog`, or
`CProd`, and `false` for plain atoms.
"""
function contains_non_simple_CFunction end
contains_non_simple_CFunction(c::T) where {T<: CFunction} = true

"""
    stringer(f::CFunction; kwargs...) -> Tuple{Bool,String}

Internal formatter returning the sign flag and body string for a coefficient
expression. Used by [`to_stringer`](@ref) and [`to_string`](@ref).
"""
function stringer end

"""
    to_stringer(f::CFunction; kwargs...) -> Tuple{Bool,String}

Wrapper around [`stringer`](@ref) that applies default formatting options.
"""
function to_stringer end

"""
    to_string(f::CFunction; kwargs...) -> String

Render `f` as a plain-text string using the coefficient formatting preferences.
"""
function to_string end
contains_non_simple_CFunction(c::CAtom)::Bool = false 
include("CFunctionsOps/ParameterValues.jl")
include("CFunctionsOps/CAtoms_variants.jl")
include("CFunctionsOps/CFunctions_Indexed.jl")
include("CFunctionsOps/CFunctions_algebra.jl")
include("CFunctionsOps/CFunctions_sort.jl")
include("CFunctionsOps/CFunctions_substitute.jl")
include("CFunctionsOps/CFunctions_simplify.jl")
include("CFunctionsOps/CFunctions_eval.jl")
include("CFunctionsOps/CFunctions_expand.jl")
include("CFunctionsOps/CFunctions_helper.jl")
include("CFunctionsOps/CFunctions_separate.jl")
include("CFunctionsOps/CFunctions_print.jl")

contains_non_simple_CFunction(c::CAtomIndexed)::Bool = false
contains_non_simple_CFunction(c::CAtomReferenced)::Bool = false
contains_non_simple_CFunction(c::CEval)::Bool = false
contains_non_simple_CFunction(c::CSum)::Bool = any(contains_non_simple_CFunction, c.expr)
# Not sure if CRational should be counted here?! -> Design choices 

end # module CFunctions
