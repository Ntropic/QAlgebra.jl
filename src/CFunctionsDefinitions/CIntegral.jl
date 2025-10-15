# NOTE: `CIntegralDefinition` lives in `src/CFunctions.jl`. This file hosts the
# helper utilities and API built around that core type.

###################################################################################################
"""
    list_cintegrals(param_info::ParameterInfo) -> Vector{AnyCIntegralDefinition}

Return all integral definitions registered in the provided [`ParameterInfo`](@ref).
"""
function list_cintegrals(param_info::ParameterInfo)
    return param_info.integral_definitions
end

"""
    define_cintegral(qspace::QSpace, expr, indexes)

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
function define_cintegral(qspace, expr::CFunction, indexes::Vector)
    param_info = qspace.param_info
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

    c_def = CIntegralDefinition(index, sortkey, expr, coerced, qspace)
    push!(param_info.integral_definitions, c_def)
    return c_def
end

function define_cintegral(qspace, expr::CFunction)
    param_info = qspace.param_info
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
    return define_cintegral(qspace, expr, indexes)
end

function compute_integral_weights!(def::CIntegralDefinition, qspace; pv=nothing)
    def.param_info === qspace.param_info ||
        error("Supplied QSpace does not match the CIntegralDefinition ParameterInfo.")
    helpers = getfield(parentmodule(@__MODULE__), :QExpressions)
    pv_local = pv === nothing ? deepcopy(qspace.sample_index_param_values) : pv
    integrand = helpers._make_integrand(def.expr, pv_local, def.assignments)
    weights = integrate_node_funs(def.interpolator, def.pdfs; f=integrand)
    def.values .= ComplexF64.(weights)
    return def
end

"""
    CIntegral

Coefficient atom referencing a registered integral definition stored inside the
owning [`ParameterInfo`](@ref).
"""
struct CIntegral{N} <: CAtomic
    param_info::ParameterInfo
    coeff::ComplexRational
    index::Int
    definition::CIntegralDefinition{N}

    function CIntegral(param_info::ParameterInfo, coeff::ComplexRational, index::Int)
        1 ≤ index ≤ length(param_info.integral_definitions) ||
            error("Integral index $index out of bounds for supplied ParameterInfo.")
        def = param_info.integral_definitions[index]
        return new{length(def.axis_lengths)}(param_info, coeff, index, def)
    end
end
CIntegral(param_info::ParameterInfo, index::Int) = CIntegral(param_info, ComplexRational(1,0,1), index)
CIntegral(def::CIntegralDefinition{N}) where {N} = CIntegral(def.param_info, ComplexRational(1,0,1), def.index)

@inline integral_definition(int::CIntegral) = int.definition
@inline integral_indexes(int::CIntegral) = int.definition.indexes
@inline integral_expr(int::CIntegral) = int.definition.expr

coeff(i::CIntegral) = [i.coeff]
var_exponents(i::CIntegral) = spzeros(Int, i.param_info.dims)
length(::CIntegral) = 1
modify_coeff(i::CIntegral, c::ComplexRational) = CIntegral(i.param_info, c, i.index)
repartition(::CIntegral, ::Vector{Tuple{Int,Int}}) = error("Cannot repartition integral definitions. Register a new integral if needed.")
