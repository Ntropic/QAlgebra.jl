# NOTE: `CIntegralDefinition` lives in `src/CFunctions.jl`. This file hosts the
# helper utilities and API built around that core type.


"""
    list_cintegrals(param_info::ParameterInfo) -> Vector{CIntegralDefinition}

Return all integral definitions registered for the supplied `ParameterInfo`.
Useful for introspection, documentation, or tooling that needs to inspect the
available coefficient integrals.
"""
function list_cintegrals(param_info::ParameterInfo)
    return param_info.integral_definitions
end

"""
    define_cintegral(param_info::ParameterInfo, expr::CFunction, indices::Vector{Vector{SubSpaceIndex}}, parameter_group_indices::Vector{Vector{Int}}) -> CIntegralDefinition
    define_cintegral(qspace::QSpace, expr::CFunction, blocks::Vector{ConstrainedIndexBlock}) -> CIntegralDefinition

Register `expr` as a coefficient integral and return its definition. Provide
the bound ensemble indices via `indices`, with matching sampling groups in
`parameter_group_indices`.
"""
function define_cintegral(param_info::ParameterInfo, expr::CFunction, indices::Vector{Vector{SubSpaceIndex}}, parameter_group_indices::Vector{Vector{Int}})
    index = length(param_info.integral_definitions) + 1
    sortkey = index + SORTKEY_BASE_CINTEGRAL
    definition = CIntegralDefinition(param_info, index, sortkey, expr, indices, parameter_group_indices, nothing, nothing)
    push!(param_info.integral_definitions, definition)
    return definition
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
        def = param_info.integral_definitions[index]
        return new(param_info, coeff, index, def)
    end
end
CIntegral(param_info::ParameterInfo, index::Int) = CIntegral(param_info, ComplexRational(1,0,1), index)
CIntegral(def::CIntegralDefinition) = CIntegral(def.param_info::ParameterInfo, ComplexRational(1,0,1), def.index)

@inline integral_definition(int::CIntegral) = int.definition
@inline integral_indices(int::CIntegral) = int.definition.indices
@inline integral_expr(int::CIntegral) = int.definition.expr

coeff(i::CIntegral) = [i.coeff]
var_exponents(i::CIntegral) = spzeros(Int, i.param_info.dims)
length(::CIntegral) = 1
modify_coeff(i::CIntegral, c::ComplexRational) = CIntegral(i.param_info, c, i.index)
repartition(::CIntegral, ::Vector{Tuple{Int,Int}}) = error("Cannot repartition integral definitions. Register a new integral if needed.")
