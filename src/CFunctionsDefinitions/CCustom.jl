###################################################################################################
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
    sortkey = index + SORTKEY_BASE_CTYPE

    abstract_parameters = abstract_from_abstractdef.(contains_which_abstracts(fun))                # defined below
    abstract_indices    = [c.index for c in abstract_parameters]
    index_map = isempty(abstract_indices) ? Int[] : begin
        m = maximum(abstract_indices)
        im = zeros(Int, m)
        for (j, ind) in enumerate(abstract_indices)
            im[ind] = j
        end
        im
    end
    has_abstract  = !isempty(abstract_parameters)
    if has_abstract && has_indices(fun)
        error("CCustomType functions either require no arguments (i.e. are deifned free of CAbstracts) or have no indices or time dependences in their definition.")
    end
    type_symbols  = (Symbol(CName), Symbol(Name), Symbol(base), :Any, :any)
    c_type_def = CTypeDefinition(name_sym, type_symbols, plain, latex, index, sortkey, fun, has_abstract, abstract_parameters, index_map, param_info)
    push!(param_info.custom_ctype, c_type_def)
    return c_type_def
end

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
length(f::CCustomType) = 1
