###################################################################################################
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
function define_cabstract(param_info::ParameterInfo, name::Union{Symbol, String})::CAbstractDefinition
    name_str, name_latex = symbol2formatted(String(name))
    for (i, abstract_def) in enumerate(param_info.abstract_definitions)
        if abstract_def.name == String(name)
            error("Abstract with name $(String(name)) already defined.")
        end
    end
    index = length(param_info.abstract_definitions) + 1
    sortkey = index + SORTKEY_BASE_CABSTRACT
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
var_exponents(a::CAbstract) = spzeros(Int, a.param_info.dims)
modify_coeff(a::CAbstract, c::ComplexRational) = CAbstract(a.param_info, c, a.index, a.exponent, a.dag)
modify_exponent(a::CAbstract, q::Rational{Int}) = CAbstract(a.param_info, a.coeff, a.index, q, a.dag)
modify_exponent(a::CAbstract, n::Int) = modify_exponent(a, n//1)
modify_dag(a::CAbstract, d::Bool=true) = CAbstract(a.param_info, a.coeff, a.index, a.exponent, d)
repartition(::CAbstract, ::Vector{Tuple{Int,Int}}) = error("You should not repartition abstract parameters! Remove them before repartitioning.")
