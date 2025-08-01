export QubitPauli

#### Cayley tables for Pauli Operators
const PAULI_TRANSFORM_TABLE = [
    4 3 2 1;
    3 4 1 2;
    2 1 4 3;
    1 2 3 4]
const PAULI_COEFF_TABLE = crationalize.([
    1 im -im 1;
    -im 1 im 1;
    im -im 1 1;
    1 1 1 1])

@doc raw"""
    QubitPauli() -> OperatorSet

Creates the OperatorSet for a qubit using Pauli operators (``\sigma_x``, ``\sigma_y``, ``\sigma_z``, ``\sigma_I``).
"""
function QubitPauli(symbol::String="")::OperatorSet
    ops = ["x", "y", "z", "I"]
    base_pauli = [1, 2, 3]
    non_base_ops::Dict{String, Vector{Tuple{ComplexRational, Is}}} = Dict("p"=> [(ComplexRational(1,0,1), 1),(ComplexRational(0,1,1), 2)], "m" => [(ComplexRational(1,0,1), 1), (ComplexRational(0,-1,1), 2)])
    
    symbol_str, symbol_latex = symbol2formatted(symbol, do_hat=true)
    do_symbol::Bool = length(symbol) > 0
    # Define the transformation function using the PAULI_* tables.
    function pauli_product(op1::Int, op2::Int)::Vector{Tuple{ComplexRational,Int}}
        # Look up the coefficient and new operator index.
        c = PAULI_COEFF_TABLE[op1, op2]
        new_index = PAULI_TRANSFORM_TABLE[op1, op2]
        # Return as a one-term sum.
        return [(c, new_index)]
    end
    function pauli_dag(op::Int)::Vector{Tuple{ComplexRational,Int}}
        return [(ComplexRational(1,0,1), op)]
    end
    function paulistr2ind(str::String)::Vector{Tuple{ComplexRational,Int}}
        if length(str) == 0
            return [(ComplexRational(1,0,1), 4)]
        end
        res = findfirst(==(str), ops)
        if isnothing(res)
            error("Invalid Pauli string: $str, must be one of $ops, or empty")
        end
        return [(ComplexRational(1,0,1), res)]
    end
    function pauli2str(ind::Int, sym::String=""; formatted::Bool=true)::String
        # create underscored string representation of sym using subscript_indexes
        if do_symbol
            if formatted
                return symbol_str * str2sup(ops[ind]) * str2sub(sym) 
            else
                return symbol_str * "_" * sym 
            end
        else
            if formatted
                return ops[ind] * str2sub(sym) 
            else
                return ops[ind] * "_" * sym
            end
        end
    end
    function pauli2latex(ind::Int, sym::String)::String
        # create underscored string representation of sym using subscript_indexes
        curr_str::String = raw""
        if do_symbol
            curr_str *= raw"{"*symbol_latex*raw"}^{" * ops[ind] * "}"
            if length(sym) > 0
                curr_str *= raw"_{" * sym * raw"}"
            end
        else
            curr_str *= raw"\hat{" * ops[ind] * "}"
            if length(sym) > 0 
                curr_str *= raw"_{" * sym * "}"
            end
        end
        return curr_str
    end
    return OperatorSet("Pauli Qubit", true, 1, 4, base_pauli, non_base_ops, ops, pauli_product, pauli_dag, pauli2str, pauli2latex)
end
## Test 
#q = QubitPauli()
#q.strs2ind("x")
#display(latexstring(q.op2latex(1, "i")))
#q.op_product(1, 2)
#q