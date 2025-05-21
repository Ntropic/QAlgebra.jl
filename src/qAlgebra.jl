module qAlgebra

using Preferences
# === Default Coefficient Preferences ===
const DEFAULT_COEFF_PREFS = Dict(
    :FLOAT_DIGITS => 2,
    :EXP_DIGITS => 2,
    :FIRST_MODE => true,
    )

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
First mode specifies whether braced terms with a leading negative are flipped.
"""
function set_first_mode(mode::Bool)
    set_default(:FIRST_MODE, mode)
end

export set_float_digits, set_exp_digits, set_first_mode

include("StringUtils.jl")

include("FFunctions.jl")
using .FFunctions
export FFunction, FAtom, FSum, FRational, simplify, isnumeric, iszero, max_exponents, build_xpows, evaluate, to_string

include("qSpace.jl")
using .qSpace
export OperatorSet, SubSpace, Parameter, OperatorType, StateSpace, string2operator_type, GLOBAL_STATE_SPACE
export QubitPauli, QubitPM, Ladder

include("qExpressions.jl")
using .qExpressions
export qExpr, qAtom, qComposite, qTerm, qEQ, qSum, Sum, ∑, diff_qEQ, term, base_operators,simplify, flatten, neq, d_dt
export Dag, Commutator, isnumeric
export string, latex_string

end # module qAlgebra