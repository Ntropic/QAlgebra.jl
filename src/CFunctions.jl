module CFunctions

using ..StringUtils
import ..SubSpaceIndex
import ..ConcreteIndexes
using ComplexRationals
using SparseArrays
using ..SparsePermutationTools: SparsePermutation
using ..QAlgebra: get_default, FLIP_IF_FIRST_TERM_NEGATIVE, DO_BRACED
using ..QDistributions: QDistribution, QEnsembleFunction

export CFunction, CAbstractDefinition, CTypeDefinition, CIntegralDefinition, ParameterInfo
export define_cabstract, define_ctype, define_cintegral
export CAbstract, CIntegral, CCustomType, CCustomTypeIndexed, CAtom, CAtomIndexed, CSum, CRational, CProd, CExp, CLog, CPower, CVector, CMatrix
export CMatrix, CVector, CPower
export coeff, var_exponents, unique_first_terms
export contains_non_simple_CFunction, Indexed, has_indexed_parameters
export list_cabstracts, list_ctypes, list_cintegrals
export where_acting, where_acting!, which_params_acting, which_params_acting!, param_index_tuples
export which_ensemble_acting, which_ensemble_acting!, substitute, separate_by_cond
export ParameterValues, set_param!, set_time!, get_parameter_index, param_value, ensure_functions!

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
    CIntegralDefinition

Container describing a coefficient integral definition stored within a
[`ParameterInfo`](@ref). It records the defining integrand `expr` and the
subsystem indexes integrated over.
"""
struct CIntegralDefinition <: CDef
    index::Int
    sortkey::Int
    expr::CFunction
    indexes::Vector{Vector{SubSpaceIndex}}
    param_info::AbstractParameterInfo
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
    outer_labels_str::Vector{String}
    outer_labels_latex::Vector{String}
    params_name::Vector{String}
    params_str::Vector{String}
    params_latex::Vector{String}

    param_of_indexes::BitVector
    param_group_by_index::Vector{Int}
    t_index_by_index::Vector{Int}
    ss_ensemble_indexes_by_group::Vector{Vector{Int}}
    ss_ensemble_present_by_group::Vector{BitVector}

    indexed_parameter_indexes::Vector{Int}
    where_acting_by_parameter::Vector{Vector{BitVector}}
    params_acting_by_index::Vector{Vector{BitVector}}
    param_index_tuples::Vector{Vector{Tuple{Int,Int}}}

    subspace_index_maps::Vector{Array{SparsePermutation,2}}
    t_index_transform::Array{SparsePermutation,2}
    indexes_by_t_index::Vector{Vector{Int}}
    indexes_of_t::Vector{Int}

    how_many_by_ensemble::Vector{Int}
    param_of_t::BitVector
    param_is_t::BitVector

    group_distributions::Vector{Union{Nothing,QDistribution}}
    group_functions::Vector{Union{Nothing,QEnsembleFunction}}
    function_param_refs::Vector{Union{Nothing,Vector{Int}}}
    group_time_counts::Vector{Int}
    group_index_sizes::Vector{Vector{Int}}
    param_coords::Vector{Vector{Int}}
    params_by_group::Vector{Vector{Int}}

    subspace_info::Any
    param_indexes::ParameterIndexes
    abstract_definitions::Vector{CAbstractDefinition}
    custom_ctype::Vector{CTypeDefinition}
    integral_definitions::Vector{CIntegralDefinition}

    function ParameterInfo(
        outer_labels_symbols::Vector{Symbol}, inner_labels_symbols_flat::Vector{Symbol}, outer_labels::Vector{String},
        outer_labels_str::Vector{String}, outer_labels_latex::Vector{String}, params_name::Vector{String},
        params_str::Vector{String}, params_latex::Vector{String}, param_of_indexes::BitVector, param_group_by_index::Vector{Int},
        t_index_by_index::Vector{Int}, ss_ensemble_indexes_by_group::Vector{Vector{Int}}, ss_ensemble_present_by_group::Vector{BitVector}, indexed_parameter_indexes::Vector{Int},
        where_acting_by_parameter::Vector{Vector{BitVector}}, params_acting_by_index::Vector{Vector{BitVector}}, param_index_tuples::Vector{Vector{Tuple{Int,Int}}},
        subspace_index_maps::Vector{Array{SparsePermutation,2}}, t_index_transform::Array{SparsePermutation,2}, indexes_by_t_index::Vector{Vector{Int}},
        indexes_of_t::Vector{Int}, how_many_by_ensemble::Vector{Int}, param_of_t::BitVector, param_is_t::BitVector,
        group_distributions::Vector{Union{Nothing,QDistribution}}, group_functions::Vector{Union{Nothing,QEnsembleFunction}}, function_param_refs::Vector{Union{Nothing,Vector{Int}}},
        group_time_counts::Vector{Int}, group_index_sizes::Vector{Vector{Int}}, param_coords::Vector{Vector{Int}}, params_by_group::Vector{Vector{Int}},
        subspace_info::Any, param_indexes::ParameterIndexes)
        dims = length(inner_labels_symbols_flat)
        new(dims, outer_labels_symbols, inner_labels_symbols_flat, outer_labels,
            outer_labels_str, outer_labels_latex,
            params_name, params_str, params_latex, param_of_indexes,
            param_group_by_index, t_index_by_index, ss_ensemble_indexes_by_group, ss_ensemble_present_by_group,
            indexed_parameter_indexes, where_acting_by_parameter, params_acting_by_index, param_index_tuples,
            subspace_index_maps, t_index_transform,
            indexes_by_t_index, indexes_of_t, how_many_by_ensemble, param_of_t, param_is_t,
            group_distributions, group_functions, function_param_refs,
            group_time_counts, group_index_sizes, param_coords, params_by_group,
            subspace_info, param_indexes, CAbstractDefinition[], CTypeDefinition[], CIntegralDefinition[])
    end
end

"""
    param_index_tuples(param_info::ParameterInfo, param_index::Int)

Return the cached `(ensemble, inner)` tuples describing where parameter
`param_index` acts. Non-indexed parameters yield an empty vector.
"""
function param_index_tuples(param_info::ParameterInfo, param_index::Int)
    1 ≤ param_index ≤ length(param_info.param_index_tuples) ||
        error("Parameter index $(param_index) out of bounds.")
    return param_info.param_index_tuples[param_index]
end

######################################################################################################################################################
"""
    list_cintegrals(param_info::ParameterInfo) -> Vector{CIntegralDefinition}

Return all integral definitions registered in the provided [`ParameterInfo`](@ref).
"""
function list_cintegrals(param_info::ParameterInfo)
    return param_info.integral_definitions
end

"""
    define_cintegral(param_info::ParameterInfo, expr, indexes)

Register a coefficient integral with integrand `expr` and integration indexes
`indexes` (outer vector per ensemble, inner vector per integrated subsystem).

`QExpressions` offers the following forwarding helpers built on top of this
method:

```
define_cintegral(qspace::QSpace, expr::QExpr, indexes::Vector{Vector{SubSpaceIndex}})
define_cintegral(qspace::QSpace, expr::QExpr)
define_cintegral(expr::QExpr, indexes::Vector{Vector{SubSpaceIndex}})
define_cintegral(expr::QExpr)
```

The overloads accepting `QExpr` perform the necessary neutrality checks and, in
the variant without explicit indexes, derive the integration blocks from the
ensembles touched by `expr` before delegating back here.
"""
function define_cintegral(param_info::ParameterInfo, expr::CFunction, indexes::Vector)::CIntegralDefinition
    expr.param_info === param_info ||
        error("Integral integrand belongs to a different ParameterInfo.")

    subspace_info = param_info.subspace_info
    subspace_info === nothing &&
        error("ParameterInfo is missing subspace information required for integrals.")

    acting = which_ensemble_acting(expr)
    n_ensembles = length(acting)
    n_expected = length(subspace_info.where_ensembles)
    n_ensembles == n_expected ||
        error("Mismatch between acting ensembles and stored subspace information.")
    length(indexes) == n_expected ||
        error("Expected $(n_expected) index groups, got $(length(indexes)).")

    coerced = Vector{Vector{SubSpaceIndex}}(undef, n_expected)
    for ensemble_idx in 1:n_expected
        group = indexes[ensemble_idx]
        coerced_group = Vector{SubSpaceIndex}(undef, length(group))
        expected_outer = subspace_info.where_ensembles[ensemble_idx]
        bits = acting[ensemble_idx]
        bitlen = length(bits)

        for (j, idx) in enumerate(group)
            idx isa SubSpaceIndex ||
                error("Index #$(j) in ensemble #$(ensemble_idx) is not a SubSpaceIndex.")
            idx.outer == expected_outer ||
                error("Index ensemble mismatch: expected outer $(expected_outer), got $(idx.outer).")
            1 ≤ idx.inner ≤ bitlen ||
                error("Index inner position $(idx.inner) out of range for ensemble #$(ensemble_idx).")
            bits[idx.inner] ||
                error("Integral index $(idx) does not appear in the integrand.")
            coerced_group[j] = idx
        end

        for inner in findall(bits)
            any(idx.inner == inner for idx in coerced_group) ||
                error("Missing integration index for ensemble #$(ensemble_idx), position $(inner).")
        end

        coerced[ensemble_idx] = coerced_group
    end

    index = length(param_info.integral_definitions) + 1
    sortkey = index + 2_000_000

    c_def = CIntegralDefinition(index, sortkey, expr, coerced, param_info)
    push!(param_info.integral_definitions, c_def)
    return c_def
end

function define_cintegral(param_info::ParameterInfo, expr::CFunction)::CIntegralDefinition
    acting = which_ensemble_acting(expr)
    subspace_info = param_info.subspace_info
    subspace_info === nothing &&
        error("ParameterInfo is missing subspace information required for integrals.")

    indexes = Vector{Vector{SubSpaceIndex}}(undef, length(acting))
    for (ensemble_idx, bits) in enumerate(acting)
        outer = subspace_info.where_ensembles[ensemble_idx]
        idxs = SubSpaceIndex[]
        for inner in findall(bits)
            push!(idxs, SubSpaceIndex(outer, inner, subspace_info))
        end
        indexes[ensemble_idx] = idxs
    end
    return define_cintegral(param_info, expr, indexes)
end


"""
    list_cabstracts(param_info::ParameterInfo) -> Vector{CAbstractDefinition}

Return all abstract symbols registered in the provided [`ParameterInfo`](@ref).
Useful for inspection and documentation purposes.
"""
function list_cabstracts(param_info::ParameterInfo)
    return param_info.abstract_definitions
end
"""
    define_cabstract(param_info::ParameterInfo, name) -> CAbstractDefinition

Register a new abstract coefficient symbol identified by `name`. The symbol is
stored inside `param_info` and can later be referenced when constructing
`CAbstract` terms.

Convenience wrappers exposed via `QExpressions` forward to this implementation
and accept a `QSpace` directly:

```
define_cabstract(qspace::QSpace, name)
```
"""
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

"""
    list_ctypes(param_info::ParameterInfo) -> Vector{CTypeDefinition}

Return every custom coefficient type defined for the given parameter info.
Each entry describes the presentation and implementation of a registered
function such as `cos` or user-defined variants.
"""
function list_ctypes(param_info::ParameterInfo)
    return param_info.custom_ctype
end
"""
    define_ctype(param_info::ParameterInfo, name, fun) -> CTypeDefinition

Register a custom coefficient function `name` whose body is given by `fun`
(a `CFunction`). The new type is available for constructing `CCustomType`
instances and is tracked inside `param_info`.

In `QExpressions` the following helper methods are provided for convenience:

```
define_ctype(qspace::QSpace, name, expr::QExpr)
define_ctype(name, expr::QExpr)
```

The QExpr overloads require `expr` to be a single-term, operator-neutral
expression; the wrapper validates these constraints and converts to the
primitive `CFunction` before delegating here.
"""
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
"""
    CIntegral

Coefficient atom referencing a registered integral definition stored inside the
owning [`ParameterInfo`](@ref).
"""
struct CIntegral <: CAtomic
    param_info::ParameterInfo
    coeff::ComplexRational
    index::Int
    definition::CIntegralDefinition

    function CIntegral(param_info::ParameterInfo, coeff::ComplexRational, index::Int)
        1 ≤ index ≤ length(param_info.integral_definitions) ||
            error("Integral index $index out of bounds for supplied ParameterInfo.")
        return new(param_info, coeff, index, param_info.integral_definitions[index])
    end
end
CIntegral(param_info::ParameterInfo, index::Int) = CIntegral(param_info, ComplexRational(1,0,1), index)
CIntegral(def::CIntegralDefinition) = CIntegral(def.param_info, ComplexRational(1,0,1), def.index)

@inline integral_definition(int::CIntegral) = int.definition
@inline integral_indexes(int::CIntegral) = int.definition.indexes
@inline integral_expr(int::CIntegral) = int.definition.expr
"""
    coeff(f::CFunction) -> Vector{ComplexRational}

Return the scalar coefficients present in `f`. For atomic objects this is the
single leading coefficient; for structured expressions the result collects the
scalars contributed by each branch.
"""
function coeff end
coeff(a::CAbstract) = [a.coeff]
coeff(i::CIntegral) = [i.coeff]
exponent(a::CAbstract) = a.exponent
"""
    var_exponents(f::CFunction) -> Vector{Int}

Return the polynomial exponents associated with each variable in `f`. Composite
objects delegate to their children, while purely numeric constructs return a
zero vector.
"""
function var_exponents end
var_exponents(a::CAbstract) = spzeros(Int, a.param_info.dims)
var_exponents(i::CIntegral) = spzeros(Int, i.param_info.dims)
isdag(a::CAbstract) = a.dag
modify_coeff(a::CAbstract, c::ComplexRational) = CAbstract(a.param_info, c, a.index, a.exponent, a.dag)
modify_exponent(a::CAbstract, q::Rational{Int}) = CAbstract(a.param_info, a.coeff, a.index, q, a.dag)
modify_exponent(a::CAbstract, n::Int) = modify_exponent(a, n//1)
modify_dag(a::CAbstract, d::Bool=true) = CAbstract(a.param_info, a.coeff, a.index, a.exponent, d)
toggle_dag(a::CAbstract) = modify_dag(a, !a.dag)
repartition(::CAbstract, ::Vector{Tuple{Int,Int}}) = error("You should not repartition abstract parameters! Remove them before repartitioning.")
length(::CIntegral) = 1
modify_coeff(i::CIntegral, c::ComplexRational) = CIntegral(i.param_info, c, i.index)
repartition(::CIntegral, ::Vector{Tuple{Int,Int}}) = error("Cannot repartition integral definitions. Register a new integral if needed.")

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
function modify_coeff(f::CCustomType, coeff::ComplexRational)::CFunction
    iszero(coeff) && return zero_catom(f.param_info)
    return CCustomType(f.param_info, coeff, f.expr, f.ctype_def)
end
var_exponents(a::CCustomType) = spzeros(Int, a.param_info.dims)
coeff(f::CCustomType) = [f.coeff]
length(f:: CCustomType) = 1

# ===================> MAIN TYPES <==========================================================================================
function modify_expr(f::CFunction, new_expr::Vector{CFunction})
    error("modify_expr not implemented for type $(typeof(f)).")
end

"""
    CAtom(param_info::ParameterInfo, var_exponents::AbstractVector{<:Int})
    CAtom(param_info::ParameterInfo, coeff::Int, var_exponents::AbstractVector{<:Int})
    CAtom(param_info::ParameterInfo, coeff::Rational, var_exponents::AbstractVector{<:Int})

A single term with a complex‐rational coefficient and integer exponents for each variable.
- The `Int` and `Rational` constructors wrap the coefficient into a `ComplexRational`.
- Exponents are stored as a sparse vector to avoid keeping zero entries.
- `var_exponents[j]` is the exponent of variable _j_.
"""
@inline function _sparse_exponents(param_info::ParameterInfo, exps)::SparseVector{Int}
    exps isa AbstractVector || return _sparse_exponents(param_info, collect(exps))
    length(exps) == param_info.dims || throw(DimensionMismatch("expected $(param_info.dims) exponents, got $(length(exps))"))
    return SparseVector{Int}(exps)
end

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
modify_exponents(a::CAtom, var_exponents) = CAtom(a.param_info, a.coeff, var_exponents)
modify_coeff(a::CAtom, coeff::ComplexRational)::CAtom = CAtom(a.param_info, coeff, a.var_exponents)
modify_coeff_exponents(a::CAtom, coeff::ComplexRational, var_exponents) = CAtom(a.param_info, coeff, var_exponents)
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

include("CFunctionsOps/CFunctions_Indexed.jl")


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
contains_non_simple_CFunction(c::CAtomIndexed)::Bool = false
contains_non_simple_CFunction(c::CSum)::Bool = any(contains_non_simple_CFunction, c.expr)
# Not sure if CRational should be counted here?! -> Design choices 


include("CFunctionsOps/CFunctions_algebra.jl")
include("CFunctionsOps/CFunctions_sort.jl")
include("CFunctionsOps/CFunctions_substitute.jl")
include("CFunctionsOps/CFunctions_simplify.jl")
include("CFunctionsOps/ParameterValues.jl")
include("CFunctionsOps/CFunctions_orders_eval.jl")
include("CFunctionsOps/CFunctions_expand.jl")
include("CFunctionsOps/CFunctions_helper.jl")
include("CFunctionsOps/CFunctions_separate.jl")
include("CFunctionsOps/CFunctions_print.jl")

end # module CFunctions
