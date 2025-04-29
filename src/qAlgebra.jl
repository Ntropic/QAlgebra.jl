module qAlgebra

using Preferences
# === Default Coefficient Preferences ===
const DEFAULT_COEFF_PREFS = Dict(
    :DEFAULT_DIGITS => 2,
    :FLOAT_DIGITS => 2,
    :EXP_DIGITS => 2)

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

export set_float_digits, set_exp_digits

include("StringUtils.jl")

include("qSpace.jl")
using .qSpace
export OperatorSet, SubSpace, Parameter, StateSpace
export QubitPauli, QubitPM, Ladder

include("qExpressions.jl")
using .qExpressions
export qExpr, qTerm, qEQ, qSum, diff_qEQ, term, simplify, base_operators, Sum, ∑, flatten, neq, d_dt, is_numeric
export Dag, Commutator
export string, latex_string

end # module qAlgebra