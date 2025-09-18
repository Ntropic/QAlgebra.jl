module CFunctions

using ..StringUtils
using ComplexRationals
using ..QAlgebra: get_default, FLIP_IF_FIRST_TERM_NEGATIVE, DO_BRACED

export CFunction, CAbstractDefinition, CTypeDefinition, ParameterInfo, add_cabstract!, add_ctype!, CAbstract, CCustomType, CAtom, CSum, CRational, CProd, CExp, CLog, CPower, CVector, CMatrix
export CMatrix, CVector, CPower
export isnumeric, coeff, var_exponents
export contains_non_simple_CFunction
export define_cabstract, define_ctype, list_cabstracts, list_ctypes
export which_ensemble_acting

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
    ParameterInfo

Holds both clusters and the dimension of the polynomial variable space.
"""
struct ParameterInfo <: AbstractParameterInfo
    dims::Int
    outer_labels_symbols::Vector{Symbol}
    inner_labels_symbols_flat::Vector{Symbol}
    
    outer_labels::Vector{String}
    params_name::Vector{String}
    params_str::Vector{String}
    params_latex::Vector{String}

    param_of_indexes::BitVector   
    outer_group_by_index::Vector{Int}   
    t_index_by_index::Vector{Int}       # -1 for parameters that aren't of t. 
    ss_ensemble_indexes_by_group::Vector{Vector{Int}}    # which ss ensembles are used for indexing in each group. 
    ss_ensemble_present_by_group::Vector{BitVector}   # which ss ensembles are present in each group.

    indexed_parameter_indexes::Vector{Int}                  # which parameters have indexes?
    where_acting_by_parameter::Vector{Vector{BitVector}}  # for each variable, where are they acting. 

    # Maps indexes for index transformation, once for switching subsystem indexes and once for time indexes
    subspace_index_maps::Vector{Array{Vector{Int},2}}
    t_index_transform::Array{Vector{Int},2}
    indexes_by_t_index::Vector{Vector{Int}}   # for each t_index which indexes have it? 
    indexes_of_t::Vector{Int}

    how_many_by_ensemble::Vector{Int}

    param_of_t::BitVector
    param_is_t::BitVector
    param_values::Vector   # specifies for example values or functions or vectors for the parameters (vectors of the index values), functions of time ...
    
    param_indexes::ParameterIndexes
    abstract_definitions::Vector{CAbstractDefinition}
    custom_ctype::Vector{CTypeDefinition}

    function ParameterInfo(
        outer_labels_symbols::Vector{Symbol}, inner_labels_symbols_flat::Vector{Symbol}, outer_labels::Vector{String}, params_name::Vector{String},
        params_str::Vector{String}, params_latex::Vector{String}, param_of_indexes::BitVector, outer_group_by_index::Vector{Int},
        t_index_by_index::Vector{Int}, ss_ensemble_indexes_by_group::Vector{Vector{Int}}, ss_ensemble_present_by_group::Vector{BitVector}, indexed_parameter_indexes::Vector{Int},
        where_acting_by_parameter::Vector{Vector{BitVector}}, subspace_index_maps::Vector{Array{Vector{Int},2}}, t_index_transform::Array{Vector{Int},2}, indexes_by_t_index::Vector{Vector{Int}},
        indexes_of_t::Vector{Int}, how_many_by_ensemble::Vector{Int}, param_of_t::BitVector, param_is_t::BitVector, param_values::Vector, param_indexes::ParameterIndexes)
        dims = length(inner_labels_symbols_flat)
        new(dims, outer_labels_symbols, inner_labels_symbols_flat, outer_labels,
            params_name, params_str, params_latex, param_of_indexes,
            outer_group_by_index, t_index_by_index, ss_ensemble_indexes_by_group, ss_ensemble_present_by_group,
            indexed_parameter_indexes, where_acting_by_parameter, subspace_index_maps, t_index_transform,
            indexes_by_t_index, indexes_of_t, how_many_by_ensemble, param_of_t, param_is_t, param_values, param_indexes, [],[])
    end
end

######################################################################################################################################################
function list_cabstracts(param_info::ParameterInfo)
    return param_info.abstract_definitions
end
function define_cabstract(param_info::ParameterInfo,  name::Union{Symbol, String})::CAbstractDefinition
    name_str, name_latex = symbol2formatted(String(name))
    for (i, abstract_def) in enumerate(param_info.abstract_definitions)
        if abstract_def.name == String(name)
            error("Abstract with name $(String(name)) already defined.")
        end
    end
    index = length(param_info.abstract_definitions) + 1
    sortkey = index + 15
    c_abstract = CAbstractDefinition(Symbol(name), name_str, name_latex, index, sortkey, param_info)
    push!(param_info.abstract_definitions, c_abstract)
    return c_abstract
end

function c_abstract_exists(param_info::ParameterInfo, name::Union{Symbol, String})::Bool 
    sym_name = Symbol(name)
    for (i, abstract_def) in enumerate(param_info.abstract_definitions)
        if abstract_def.symbol == sym_name 
            return true 
        end
    end
    return false 
end

function list_ctypes(param_info::ParameterInfo)
    return param_info.custom_ctype
end
function define_ctype(param_info::ParameterInfo, name::Union{Symbol,String}, fun::CFunction)::CTypeDefinition
    CName, Name, base = variants_C(name)
    name_sym = Symbol(base)
    # check if name_sym is already present in custom_ctype 
    if any(x -> x.name == name_sym, param_info.custom_ctype) 
        error("Cannot define $name_sym, because it already exists in ParameterInfo.")
    end
    plain, latex = symbol2formatted(String(base))

    index   = length(param_info.custom_ctype) + 1
    sortkey = index + 10^6

    abstract_parameters = abstract_from_abstractdef.(contains_which_abstracts(fun))                # defined below
    abstract_indexes    = [c.index for c in abstract_parameters]
    index_map = isempty(abstract_indexes) ? Int[] : begin
        m = maximum(abstract_indexes)
        im = zeros(Int, m)
        for (j, ind) in enumerate(abstract_indexes)
            im[ind] = j
        end
        im
    end
    has_abstract  = !isempty(abstract_parameters)
    if has_abstract && has_indexes(fun)
        error("CCustomType functions either require no arguments (i.e. are deifned free of CAbstracts) or have no indexes or time dependences in their definition.")
    end
    type_symbols  = (Symbol(CName), Symbol(Name), Symbol(base), :Any, :any)
    c_type_def = CTypeDefinition(name_sym, type_symbols, plain, latex, index, sortkey, fun, has_abstract, abstract_parameters, index_map, param_info)
    push!(param_info.custom_ctype, c_type_def)
    return c_type_def
end



# =====================================================> CFunction Types <=====================================================================================================
"""
    CAbstract

Abstract symbol instance (optionally daggered and/or with an integer/rational power)
with a complex-rational coefficient:

    coeff * A_index^(exponent)  (daggered if dag=true)

Fields
- `param_info` : ParameterInfo
- `coeff`      : ComplexRational
- `index`      : Int (1-based index into `param_info.abstract_definitions`)
- `exponent`   : Rational{Int} (use `n//1` for integer n)
- `dag`        : Bool
- `abstract_def` : CAbstractDefinition (back-reference convenience)
"""
struct CAbstract <: AbstractCAbstract
    param_info::ParameterInfo
    coeff::ComplexRational
    index::Int
    exponent::Rational{Int}
    dag::Bool
    abstract_def::CAbstractDefinition

    # Core inner constructors
    function CAbstract(param_info::ParameterInfo, coeff::ComplexRational, index::Int, exponent::Rational{Int}=1//1, dag::Bool=false)
        return new(param_info, coeff, index, exponent, dag, param_info.abstract_definitions[index])
    end
end
coeff(a::CAbstract) = [a.coeff]
exponent(a::CAbstract) = a.exponent
var_exponents(a::CAbstract) = zeros(Int, a.param_info.dims)
isdag(a::CAbstract) = a.dag
modify_coeff(a::CAbstract, c::ComplexRational) = CAbstract(a.param_info, c, a.index, a.exponent, a.dag)
modify_exponent(a::CAbstract, q::Rational{Int}) = CAbstract(a.param_info, a.coeff, a.index, q, a.dag)
modify_exponent(a::CAbstract, n::Integer) = modify_exponent(a, n//1)
modify_dag(a::CAbstract, d::Bool=true) = CAbstract(a.param_info, a.coeff, a.index, a.exponent, d)
toggle_dag(a::CAbstract) = modify_dag(a, !a.dag)
repartition(::CAbstract, ::Vector{Tuple{Int,Int}}) = error("You should not repartition abstract parameters! Remove them before repartitioning.")

"""
    CCustomType(param_info, def_id, coeff, x)

Instance of a parametric custom function: `coeff * name(x)`.
"""
struct CCustomType <: CFunction
    param_info::ParameterInfo
    coeff::ComplexRational
    expr::Vector{CFunction}
    ctype_def::CTypeDefinition
end

function repartition(f::CCustomType, var_tuples::Vector{Tuple{Int, Int}})::CCustomType
    new_parameters = repartition.(f.expr, Ref(var_tuples))
    return CCustomType(f.param_info, f.coeff, new_parameters, f.ctype_def)
end
function modify_expr(f::CCustomType, new_expr::Vector{CFunction})
    return CCustomType(f.param_info, f.coeff, new_expr, f.ctype_def)
end
var_exponents(a::CCustomType) = zeros(Int, a.param_info.dims)
coeff(f::CCustomType) = [f.coeff]
length(f:: CCustomType) = 1

# ===================> MAIN TYPES <==========================================================================================
function modify_expr(f::CFunction, new_expr::Vector{CFunction})
    error("modify_expr not implemented for type $(typeof(f)).")
end

"""
    CAtom(coeff::Int, var_exponents::Vector{Int})
    CAtom(coeff::Rational, var_exponents::Vector{Int})
    CAtom(coeff::ComplexRational, var_exponents::Vector{Int})

A single term with a complex‐rational coefficient and integer exponents for each variable.
- The `Int` and `Rational` constructors wrap the coefficient into a `ComplexRational`.
- `var_exponents[j]` is the exponent of variable _j_.
"""
struct CAtom <: CAtomic
    param_info::ParameterInfo
    coeff::ComplexRational
    var_exponents::Vector{Int}
    function CAtom(param_info::ParameterInfo, var_exponents::Vector{Int})
        c = ComplexRational(1, 0, 1)
        return new(param_info, c, copy(var_exponents))
    end
    function CAtom(param_info::ParameterInfo, coeff::Int, var_exponents::Vector{Int})
        c = ComplexRational(coeff, 0, 1)
        return new(param_info, c, copy(var_exponents))
    end
    function CAtom(param_info::ParameterInfo, coeff::Rational, var_exponents::Vector{Int})
        c = ComplexRational(numerator(coeff), 0, denominator(coeff))
        return new(param_info, c, copy(var_exponents))
    end
    function CAtom(param_info::ParameterInfo, coeff::Complex, var_exponents::Vector{Int})
        c = crationalize(coeff)
        return new(param_info, c, copy(var_exponents))
    end
    function CAtom(param_info::ParameterInfo, coeff::ComplexRational, var_exponents::Vector{Int})
        return new(param_info, coeff, copy(var_exponents))
    end
    function CAtom(param_info::ParameterInfo, coeff::Number, var_exponents::Vector{Int})
        c = crationalize(coeff+0im)
        return new(param_info, c, copy(var_exponents))
    end
end
coeff(a::CAtom)::Vector{ComplexRational} = [a.coeff]
modify_exponents(a::CAtom, var_exponents::Vector{Vector{Int}})::CAtom = CAtom(a.param_info, a.coeff, var_exponents)
modify_coeff(a::CAtom, coeff::ComplexRational)::CAtom = CAtom(a.param_info, coeff, a.var_exponents)
modify_coeff_exponents(a::CAtom, coeff::ComplexRational, var_exponents::Vector{Vector{Int}}) = CAtom(a.param_info, coeff, var_exponents)
var_exponents(a::CAtom) = a.var_exponents
length(a::CAtom) = 1
function repartition(f::CAtom, var_tuples::Vector{Tuple{Int, Int}})::CAtom 
    curr_var_exponents = f.var_exponents
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
coeff(x::CSum) = [ComplexRational(1,0,1)] #error("Sums don't have a coeff, you likely have a sum in a sum, this shouldn't happen. Please inform the developers. ")
length(q::CSum) = length(q.expr)
repartition(f::CSum, var_tuples::Vector{Tuple{Int, Int}}) = _CSum(f.param_info, repartition.(f.expr, Ref(var_tuples)) )
var_exponents(a::CSum) = min.(var_exponents.(a.expr)...)

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
coeff(x::CProd) = [x.coeff]
length(q::CProd) = max(length.(q.expr)...)
repartition(f::CProd, var_tuples::Vector{Tuple{Int, Int}})= CProd(f.param_info, f.coeff, repartition.(f.expr, Ref(var_tuples)) )
function var_exponents(a::CProd) 
    if length(a.expr) > 0 
        return var_exponents(a.expr[1])
    else 
        return zeros(Int, a.param_info.dims)
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
var_exponents(a::CRational) = var_exponents(a.numer)


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
coeff(x::CExp) = [x.coeff]
length(q::CExp) = 1
repartition(q::CExp, var_tuples::Vector{Tuple{Int, Int}}) = CExp(q.param_info, q.coeff, repartition(q.expr, var_tuples))
var_exponents(a::CExp) = zeros(Int, a.param_info.dims)


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
coeff(x::CLog) = [x.coeff] 
length(q::CLog) = 1
repartition(q::CLog, var_tuples::Vector{Tuple{Int, Int}}) = CLog(q.param_info, q.coeff, repartition(q.expr, var_tuples))
var_exponents(a::CExp) = zeros(Int, a.param_info.dims)

"""
    CPower(coeff::ComplexRational, x::CFunction, exponent::Rational{Int})
    CPower(x::CFunction, exponent::Rational{Int})
    CPower(x::CFunction, exponent::Integer)

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
CPower(param_info::ParameterInfo, expr::CFunction, n::Integer)       = CPower(param_info, ComplexRational(1,0,1), expr, n//1)
CPower(param_info::ParameterInfo, expr::CFunction, q::Rational{Int}) = CPower(param_info, ComplexRational(1,0,1), expr, q)
function modify_expr(f::CPower, new_expr::Vector{CFunction})
    @assert length(new_expr) == 1
    return CPower(f.param_info, f.coeff, new_expr[1], f.exponent, Val(:nosimp))
end
coeff(p::CPower) = [p.coeff]
length(::CPower) = 1
repartition(p::CPower, var_tuples::Vector{Tuple{Int,Int}}) = CPower(p.param_info, p.coeff, repartition(p.expr, var_tuples), p.exponent)
var_exponents(a::CExp) = zeros(Int, a.param_info.dims)


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
coeff(v::CVector) = isempty(v.expr) ? ComplexRational[] : vcat(coeff.(v.expr)...)
length(v::CVector) = length(v.expr)
size(v::CVector) = v.row ? (1, length(v.expr)) : (length(v.expr), 1)
getindex(v::CVector, i::Int) = v.expr[i]
iterate(v::CVector, st::Int=1) = st > length(v.expr) ? nothing : (v.expr[st], st+1)
repartition(v::CVector, var_tuples::Vector{Tuple{Int,Int}}) = CVector(v.param_info, v.coeff, repartition.(v.expr, Ref(var_tuples)); row=v.row)
var_exponents(a::CExp) = zeros(Int, a.param_info.dims)


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
coeff(M::CMatrix) = [M.coeff]
length(M::CMatrix) = length(M.expr)         # number of elements (m*n)
size(M::CMatrix) = size(M.expr)
getindex(M::CMatrix, i::Int, j::Int) = M.expr[i, j]
repartition(M::CMatrix, var_tuples::Vector{Tuple{Int,Int}}) = CMatrix(M.param_info, M.coeff, reshape(repartition.(M.expr[:], Ref(var_tuples)), size(M.expr)))
var_exponents(a::CExp) = zeros(Int, a.param_info.dims)


#### Some basic functions ##############################################################################################



import Base: length, getindex, iterate, deleteat!, reverse
length(p::CFunction)::Int = 1

getindex(p::CSum, i::Int) = p.expr[i]
iterate(p::CSum, state=1) = state > length(p.expr) ? nothing : (p.expr[state], state + 1)
deleteat!(p::CSum, i::Int) = _CSum(deleteat!(p.expr, i))
reverse(q::CSum) = CSum(reverse(q.expr))

"""
    contains_non_simple_CFunction(c::CFunction) -> Bool 

Does the expression contain non simple classical functions, such as CExp, CLog, CProd? 
"""
contains_non_simple_CFunction(c::T) where {T<: CFunction} = true
contains_non_simple_CFunction(c::CAtom)::Bool = false 
contains_non_simple_CFunction(c::CSum)::Bool = any(contains_non_simple_CFunction, c.expr)
# Not sure if CRational should be counted here?! -> Design choices 


include("CFunctionsOps/CFunctions_algebra.jl")
include("CFunctionsOps/CFunctions_sort.jl")
include("CFunctionsOps/CFunctions_simplify.jl")
include("CFunctionsOps/CFunctions_orders_eval.jl")
include("CFunctionsOps/CFunctions_expand.jl")
include("CFunctionsOps/CFunctions_helper.jl")
include("CFunctionsOps/CFunctions_print.jl")

end # module CFunctions