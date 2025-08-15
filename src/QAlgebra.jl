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

include("StringUtils.jl")
using .StringUtils
export symbol2formatted, str2sub, str2sup, brace

include("CFunctions.jl")
using .CFunctions
export CFunction, CAtom, CSum, CRational, CProd, CExp, CLog
export reorder, max_exponents, build_xpows, evaluate, stringer, to_stringer, to_string, sort_key
export isnumeric, coeff, var_exponents, expand

include("QSpace.jl")
using .QSpace
export OperatorSet, SubSpace, Parameter, OperatorType, StateSpace, string2operator_type, GLOBAL_STATE_SPACE
export QubitPauli, QubitPM, Ladder

include("QExpressions.jl")
using .QExpressions
export QEq, QObj, QAtom, QAbstract, QComposite, QCompositeN, QCompositeProduct, QMultiComposite, QTerm, QAtomProduct, QExpr, QSum, Sum, ∑, diff_QEq, base_operators, neq, d_dt
export QCommutator, QExp, QLog, QPower, power, QRoot, root, simplify, simplify_QAtomProduct
export Dag, Commutator, is_numeric, same_statespace
export reorder
export string, latex_string
export term

export contains_abstract, which_continuum_acting, are_indexes_defined
export substitute
export reorder!
export QExpr2string   # remove later 
end # module QAlgebra