module QAlgebra

using Preferences
# === Default Coefficient Preferences ===
const DEFAULT_COEFF_PREFS = Dict(
    :FLIP_IF_FIRST_TERM_NEGATIVE  => true,
    :DO_BRACED => true,
    :EXPAND_CUMULANTS => false
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
EXPAND_CUMULANTS = get_default(:EXPAND_CUMULANTS)

@inline function _update_pref!(name::Symbol, value)
    set_default(name, value)
    @eval $(Symbol(name)) = $value
end

"""
    set_flip_if_first_term_negative(mode::Bool)

Sets a new default value for the first mode and saves it persistently.
First mode specifies whether braced terms with a leading negative are flipped or only if all terms are negative.
"""
function set_flip_if_first_term_negative(mode::Bool)
    _update_pref!(:FLIP_IF_FIRST_TERM_NEGATIVE, mode)
end

"""
    set_do_braced(b::Bool)
Sets a new default value for :DO_BRACED. Toggles whether terms are grouped when printing them, into groups with common coefficients. 
"""
function set_do_braced(mode::Bool)
    _update_pref!(:DO_BRACED, mode)
end

"""
    set_expand_cumulants(mode::Bool)

Set a new default for whether cumulants print in expanded form (`true`) or compact form (`false`).
The preference persists via `Preferences.jl`.
"""
function set_expand_cumulants(mode::Bool)
    _update_pref!(:EXPAND_CUMULANTS, mode)
end



export get_default, set_flip_if_first_term_negative, set_do_braced, set_expand_cumulants,
       FLIP_IF_FIRST_TERM_NEGATIVE , DO_BRACED, EXPAND_CUMULANTS

include("Helper.jl")
include("OffsetArrays.jl")
export OffsetArray

include("StringUtils.jl")
using .StringUtils
export symbol2formatted, str2sub, str2sup, brace, braket, indexes2str
export int_exponent2str, exponentdag2str

include("CFunctions.jl")
using .CFunctions
export CFunction, CAbstractDefinition, CTypeDefinition, ParameterInfo, add_cabstract!, add_ctype!, CAbstract, CCustomType, CAtom, CSum, CRational, CProd, CExp, CLog
export CMatrix, CVector, CPower
export reorder, max_exponents, build_xpows, evaluate, stringer, to_stringer, to_string, sort_key
export coeff, var_exponents, expand, substitute
export contains_non_simple_CFunction
export define_cabstract, define_ctype, list_cabstracts, list_ctypes
export where_acting, where_acting!, which_ensemble_acting, which_ensemble_acting!
export which_params_acting, which_params_acting!, separate_by_cond

include("Cumulants.jl")
using .Cumulants

include("QSpace.jl")
using .QSpaces
export OperatorSet, max_operator_magnitude
export Ensemble, SubSpace, SubSpaceDefinitions 
export OperatorType, OperatorTypeInfo, OperatorDefinitions
export Parameter, ParameterInfo, ParameterDefinitions
export QSpace
export QubitPauli, QubitPM, Ladder

include("QExpressions.jl")
using .QExpressions
export QObj, QAtom, QComposite, QCompositeN, QMultiComposite, QAbstract, QTerm, QExpr, QCumulant, QAtomProduct, permutation, AbstractQSum, QSum, QIntegral, QInt, ∑, ∫, integral, QCompositeProduct, diffQEq, d_dt
export QCommutator, QExp, QLog, QPower, power, QRoot, root, Expectation
export Dag, Commutator, same_qspace
export is_t_var, is_local, contains_non_simple_QObj, contains_non_simple, contains_abstract, contains_time, max_order_of_terms, contains_which_t_indexes, iter_QAtomProducts, iter_QInts
export is_unitary, is_hermitian, substitution_properties_fulfilled
export base_operators, QExprLookup 
export string, latex_string
export term
export @define, @define_basics, QExpr2CFunction, Cumulant, cumulant_string

export contains_abstract, are_indexes_defined, which_summations_acting, which_summations_to_root, are_all_neq
export Substitution, --> 
export reorder, reorder_full, reorder_time, neq, flatsums, complexsums , Sum2Int
# Preindexing
export QAtomOrdered, QNeutral, OrderbyOperator

#include("Indexing.jl")
#include("Combinatorics.jl")

#include("QEqSets.jl")
#using .QEqSets
export diffQEqSet, diffQEqSetOrdered, diffQEqSet, OrderedDiffQEqSet
end # module QAlgebra
