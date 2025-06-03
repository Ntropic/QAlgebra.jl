module qAlgebra

using Preferences
# === Default Coefficient Preferences ===
const DEFAULT_COEFF_PREFS = Dict(
    :FLOAT_DIGITS => 2,
    :EXP_DIGITS => 2,
    :FIRST_MODE => true,
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

""" 
    set_float_digits(d::Int)

Sets a new default value for the digit precision of floating point coefficients and saves it persistently.
"""
function set_float_digits(d::Int)
    set_default(:FLOAT_DIGITS, d)
end
"""
    set_exp_digits(d::Int)
Sets a new default value for the digit precision of exponential coefficients and saves it persistently.
"""
function set_exp_digits(d::Int)
    set_default(:EXP_DIGITS, d)
end
"""
    set_first_mode(mode::Bool)
Sets a new default value for the first mode and saves it persistently.
First mode specifies whether braced terms with a leading negative are flipped or only if all terms are negative.
"""
function set_first_mode(mode::Bool)
    set_default(:FIRST_MODE, mode)
end

export get_default, set_float_digits, set_exp_digits, set_first_mode

include("StringUtils.jl")

include("FFunctions.jl")
using .FFunctions
export FFunction, FAtom, FSum, FRational, simplify, isnumeric, iszero, max_exponents, build_xpows, evaluate, stringer, to_stringer, to_string

include("qSpace.jl")
using .qSpace
export OperatorSet, SubSpace, Parameter, OperatorType, StateSpace, string2operator_type, GLOBAL_STATE_SPACE
export QubitPauli, QubitPM, Ladder

include("qExpressions.jl")
using .qExpressions
export qObj, qAtom, qComposite, qTerm, qProd, qExpr, qSum, Sum, ∑, diff_qEQ, term, base_operators, flatten, neq, d_dt
export Dag, Commutator, isnumeric
export string, latex_string

end # module qAlgebra