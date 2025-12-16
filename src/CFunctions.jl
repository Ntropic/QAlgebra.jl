module CFunctions

using ..StringUtils
using ..QIndexes: QIndex, AbstractIndex, TimeIndex
import ..SubSpaceIndex
import ..ConcreteIndexes
using ComplexRationals
using SparseArrays
using ..SparsePermutationTools: SparsePermutation
using ..QAlgebra: get_default, FLIP_IF_FIRST_TERM_NEGATIVE, DO_BRACED
using ..Sampler: QInterpolator, QIntegrator
using ..ParameterGroups: ParameterGroup, ParameterGroupDistribution, ParameterGroupStorageUnion, ParameterGroupLike, AbstractSubSpaceInfo

export CFunction, CAbstractDefinition, CTypeDefinition, CIntegralDefinition, ParameterInfo
export define_cabstract, define_ctype, define_cintegral
export CAbstract, CIntegral, CCustomType, CCustomTypeIndexed, CParticle, CAtom, CAtomIndexed, CAtomReferenced, CEval, CSum, CRational, CProd, CExp, CLog, CPower, CVector, CMatrix
export CMatrix, CVector, CPower
export coeff, var_exponents, unique_first_terms
export contains_non_simple_CFunction, Indexed
export list_cabstracts, list_ctypes, list_cintegrals
export where_acting, where_acting!, which_params_acting, which_params_acting!
export which_ensemble_acting, which_ensemble_acting!, substitute, separate_by_cond
export ParameterValues, update_t!, resolve_param!, resolve_ensemble_values!, refresh_ensemble_values!
export compute_integral_weights!

import Base: copy, exp, log, length, getindex, iterate, size
import ComplexRationals: isonelike
import ..QAlgebra: vecvec_or, vecvec_or!, sort_unique!, variants_C

const CR_ZERO = ComplexRational(0,0,1)
const CR_ONE  = ComplexRational(1,0,1)

const SORTKEY_BASE_CABSTRACT = 1_000_000
const SORTKEY_BASE_CTYPE     = 2_000_000
const SORTKEY_BASE_CINTEGRAL = 3_000_000


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
"""
    CAbstractDefinition

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
    param_info::Any
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
    param_info::Any
end

struct ParameterIndexes end
"""
    CIntegralDefinition

Container describing a coefficient integral definition stored within a
[`ParameterInfo`](@ref). It records the defining integrand `expr`, the
subsystem indices integrated over, and the ensemble parameter groups relevant
for the integration.
"""
struct CIntegralDefinition <: CDef
    param_info::Any
    index::Int
    sortkey::Int
    expr::CFunction
    indices::Vector{Vector{SubSpaceIndex}}
    parameter_group_indices::Vector{Vector{Int}}
    interpolator::Union{Nothing,QInterpolator}
    integrator::Union{Nothing,QIntegrator}
end

"""
    ParameterInfo

Holds both clusters and the dimension of the polynomial variable space.
"""
struct ParameterInfo
    params_symbols::Vector{Symbol}
    params_raw::Vector{String}
    params_str::Vector{String}
    params_latex::Vector{String}

    subspace_info::AbstractSubSpaceInfo
    param_groups::Vector{ParameterGroupLike}

    abstract_definitions::Vector{CAbstractDefinition}
    custom_ctype::Vector{CTypeDefinition}
    integral_definitions::Vector{CIntegralDefinition}

    function ParameterInfo(
        params_symbols::Vector{Symbol},
        params_raw::Vector{String},
        params_str::Vector{String},
        params_latex::Vector{String},
        subspace_info::AbstractSubSpaceInfo,
        param_groups::Vector{ParameterGroupLike})
        new(params_symbols, params_raw, params_str, params_latex,
            subspace_info, param_groups,
            CAbstractDefinition[], CTypeDefinition[], CIntegralDefinition[])
    end
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

# ===================> MAIN TYPES <=========================================================================================

function modify_expr(f::CFunction, new_expr::Vector{CFunction})
    error("modify_expr not implemented for type $(typeof(f)).")
end


"""
    CParticle{T}(group_index, exponent, indices, time_index)

Defines a single parameter, with its exponent and indices.
"""
struct CParticle{T<:QIndex}
    group_index::Int
    exponent::Int
    abstract_indices::Vector{T}
    time_index::TimeIndex
end

include("CFunctionsOps/CParticles_helper.jl")


"""
    CAtom(param_info::ParameterInfo, particles)
    CAtom(param_info::ParameterInfo, coeff::Number, particles)

Single polynomial atom with a complex-rational coefficient and explicit particle
content.
"""
struct CAtom <: CAtomic
    param_info::ParameterInfo
    coeff::ComplexRational
    particles::Vector{CParticle{AbstractIndex}}
    function CAtom(param_info::ParameterInfo, particles::Vector{CParticle{AbstractIndex}})
        return new(param_info, CR_ONE, particles)
    end
    function CAtom(param_info::ParameterInfo, coeff::Int, particles::Vector{CParticle{AbstractIndex}})
        rational_coeff = ComplexRational(coeff, 0, 1)
        return new(param_info, rational_coeff, particles)
    end
    function CAtom(param_info::ParameterInfo, coeff::Rational, particles::Vector{CParticle{AbstractIndex}})
        rational_coeff = ComplexRational(numerator(coeff), 0, denominator(coeff))
        return new(param_info, rational_coeff, particles)
    end
    function CAtom(param_info::ParameterInfo, coeff::Complex, particles::Vector{CParticle{AbstractIndex}})
        rational_coeff = crationalize(coeff)
        return new(param_info, rational_coeff, particles)
    end
    function CAtom(param_info::ParameterInfo, coeff::ComplexRational, particles::Vector{CParticle{AbstractIndex}})
        return new(param_info, coeff, particles)
    end
    function CAtom(param_info::ParameterInfo, coeff::Number, particles::Vector{CParticle{AbstractIndex}})
        rational_coeff = crationalize(coeff + 0im)
        return new(param_info, rational_coeff, particles)
    end
end
@inline function zero_catom(param_info::ParameterInfo)
    return CAtom(param_info, CR_ZERO, CParticle{AbstractIndex}[])
end

coeff(a::CAtom)::Vector{ComplexRational} = [a.coeff]
modify_coeff(a::CAtom, coeff::ComplexRational)::CAtom = CAtom(a.param_info, coeff, a.particles)
var_exponents(a::CAtom)::Vector{Int} = Int[p.exponent for p in a.particles]
length(a::CAtom) = 1
function repartition(f::CAtom, var_tuples::Vector{Tuple{Int, Int}})::CAtom 
    new_particles = copy(f.particles)
    exps = Int[p.exponent for p in new_particles]
    @inbounds for (source, target) in var_tuples
        exps[target] += exps[source]
        exps[source] = 0
    end
    for (idx, exp) in enumerate(exps)
        part = new_particles[idx]
        new_particles[idx] = CParticle(part.group_index, exp, part.abstract_indices, part.time_index)
    end
    filter!(p -> p.exponent != 0, new_particles)
    return CAtom(f.param_info, f.coeff, new_particles)
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
