module QAlgebra

using Preferences
# === Default Coefficient Preferences ===
const DEFAULT_COEFF_PREFS = Dict(
    :FLIP_IF_FIRST_TERM_NEGATIVE  => true,
    :DO_BRACED => true
    )

""" 
    get_default(name::Symbol)

Returns the current default value for the coefficient preference with the given name. If no such preference has been set 
it returns the default value from `DEFAULT_COEFF_PREFS`.
"""
function get_default(name::Symbol)
    return @load_preference(String(name), DEFAULT_COEFF_PREFS[name])
end
function set_default(name::Symbol, value)
    @set_preferences!(String(name) => value)
end
function default_if_nothing(x, key)
    return x === nothing ? get_default(key) : x
end

# --- Module-level global variables (not const) ---
# These are initialized immediately at module load time
FLIP_IF_FIRST_TERM_NEGATIVE  = get_default(:FLIP_IF_FIRST_TERM_NEGATIVE )
DO_BRACED = get_default(:DO_BRACED)

"""
    set_flip_if_first_term_negative(mode::Bool)
Sets a new default value for the first mode and saves it persistently.
First mode specifies whether braced terms with a leading negative are flipped or only if all terms are negative.
"""
function set_flip_if_first_term_negative(mode::Bool)
    set_default(:FLIP_IF_FIRST_TERM_NEGATIVE , mode)
    @eval $(Symbol(:FLIP_IF_FIRST_TERM_NEGATIVE)) = $mode
end

"""
    set_do_braced(b::Bool)
Sets a new default value for :DO_BRACED. Toggles whether terms are grouped when printing them, into groups with common coefficients. 
"""
function set_do_braced(mode::Bool)
    set_default(:DO_BRACED, mode)
    @eval $(Symbol(:DO_BRACED)) = $mode
end



export get_default, set_flip_if_first_term_negative, set_do_braced, FLIP_IF_FIRST_TERM_NEGATIVE , DO_BRACED

include("Helper.jl")

include("StringUtils.jl")
using .StringUtils
export symbol2formatted, str2sub, str2sup, brace, braket, indexes2str
export int_exponent2str, exponentdag2str

include("CFunctions.jl")
using .CFunctions
export CFunction, CAbstractDefinition, CTypeDefinition, ParameterInfo, add_cabstract!, add_ctype!, CAbstract, CCustomType, CAtom, CSum, CRational, CProd, CExp, CLog
export CMatrix, CVector, CPower
export repartition, max_exponents, build_xpows, evaluate, stringer, to_stringer, to_string, sort_key
export coeff, var_exponents, expand, substitute
export contains_non_simple_CFunction
export define_cabstract, define_ctype, list_cabstracts, list_ctypes
export which_ensemble_acting

include("QSpace.jl")
using .QSpace
export OperatorSet
export SubSpace, SubSpaceDefinitions 
export OperatorType, OperatorTypeInfo, OperatorDefinitions
export Parameter, ParameterInfo, ParameterDefinitions
export StateSpace
export QubitPauli, QubitPM, Ladder

include("QExpressions.jl")
using .QExpressions
export QObj, QAtom, QComposite, QCompositeN, QMultiComposite, QAbstract, QTerm, QExpr, QAtomProduct, QSum, Sum, ∑, QCompositeProduct, diff_QEq, d_dt
export QCommutator, QExp, QLog, QPower, power, QRoot, root #, simplify
export Dag, Commutator, same_statespace
export is_t_var, is_local, contains_non_simple_QObj, contains_non_simple, contains_abstract, contains_time, max_moment_of_terms, contains_which_t_indexes, where_acting
export is_unitary, is_hermitian, substitution_properties_fulfilled
export base_operators, QExprLookup 
export string, latex_string
export term
export @define, @define_basics, QExpr2CFunction

export contains_abstract, are_indexes_defined, which_summations_acting, which_summations_to_root
export Substitution, --> 
export repartition!,repartition, neq, flatsums, complexsums
export QExpr2string   # remove later 
end # module QAlgebra

