export QubitPM

#### Cayley tables for PM Operators
const PM_TRANSFORM_TABLE1 = Int[4 4 1 1; # p 
    4 4 2 2; # m
    1 2 4 3; # z
    1 2 3 4] # I
#                               p,m,z,I
const PM_TRANSFORM_TABLE2 = Int[0 3 0 0; # p 
    3 0 0 0; # m
    0 0 0 0; # z
    0 0 0 0] # I

const PM_COEFF_TABLE1 = crationalize.([0.0 0.5 -1.0 1.0;
    0.5 0.0 1.0 1.0;
    1.0 -1.0 1.0 1.0;
    1.0 1.0 1.0 1.0])
const PM_COEFF_TABLE2 = crationalize.([0 0.5 0 0;
    -0.5 0 0 0;
    0 0 0 0;
    0 0 0 0])
const PM_HOW_MANY = Int[1 2 1 1; # p 
    2 1 1 1; # m
    1 1 1 1; # z
    1 1 1 1] # I

@doc raw""" 
    QubitPM() -> OperatorSet

Creates the OperatorSet for a qubit using Raising and Lowering operators (``\sigma_+``, ``\sigma_-``, ``\sigma_z``, ``\sigma_I``).
"""
function QubitPM(symbol::String="")::OperatorSet
    ops = ["p", "m", "z", "I"]
    ops_str = ["+", "-", "z", "I"]
    non_base_ops::Dict{String, Vector{Tuple{ComplexRational, Vector{Int}}}} = Dict("x"=> [(ComplexRational(1,0,1), [1]),(ComplexRational(1,0,1), [2])], "y" => [(ComplexRational(0,-1,1), [1]), (ComplexRational(0,1,1), [2])])
    base_pm = [[1], [2], [3]]
    pm_dag_inds = [[2], [1], [3], [4]]

    symbol_str, symbol_latex = symbol2formatted(symbol, do_hat=true)
    do_symbol::Bool = length(symbol) > 0
    # Define the transformation function using the PAULI_* tables.
    function pm_product(op1s::Vector{Int}, op2s::Vector{Int})::Vector{Tuple{ComplexRational,Vector{Int}}}
        # Look up the coefficient and new operator index.
        op1 = op1s[1]
        op2 = op2s[1]
        c1 = PM_COEFF_TABLE1[op1, op2]
        new_index1 = [PM_TRANSFORM_TABLE1[op1, op2]]
        if PM_HOW_MANY[op1, op2] == 2
            c2 = PM_COEFF_TABLE2[op1, op2]
            new_index2 = [PM_TRANSFORM_TABLE2[op1, op2]]
            return [(c1, new_index1), (c2, new_index2)]
        else
            # Return as a one-term sum.
            return [(c1, new_index1)]
        end
    end
    function pm_dag(op::Vector{Int})::Vector{Tuple{ComplexRational,Vector{Int}}}
        return [(ComplexRational(1,0,1), pm_dag_inds[op[1]])]
    end
    function pm2str(inds::Vector{Int}, sym::String=""; formatted::Bool=true)::String
        # create underscored string representation of sym using subscript_indexes
        ind = inds[1]
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
    function pm2latex(inds::Vector{Int}, sym::String)::String
        # create underscored string representation of sym using subscript_indexes
        ind = inds[1]
        curr_str::String = raw""
        if do_symbol
            curr_str *= symbol_latex*raw"^{" * ops_str[ind] * "}"
            if length(sym) > 0
                curr_str *= raw"_{" * sym * "}"
            end
        else
            curr_str *= raw"\hat{" * ops[ind] * "}"
            if length(sym) > 0 
                curr_str *= raw"_{" * sym * "}"
            end
        end
        return curr_str
    end
    function pmcommutes(op1::Vector{Int}, op2::Vector{Int})::Bool
        # everything commutes with 4, otherwise must be the same
        return (op1[1] == 4 || op2[1] == 4 || op1[1] == op2[1])
    end
    return OperatorSet("PM Qubit", "Fermion", 1, [4], base_pm, non_base_ops, ops, pm_product, pm_dag, pm2str, pm2latex, pmcommutes)
end
# Test 
#q = QubitPM()
#q.strs2ind("p")
#display(latexstring(q.op2latex(1, "i")))
#q.op_product(1, 2)