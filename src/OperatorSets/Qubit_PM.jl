export QubitPM

#### Cayley tables for PM Operators
#                               I, p, m, z
const PM_TRANSFORM_TABLE1 = Int[1 2 3 4; # I
                                2 1 1 2; # p
                                3 1 1 3; # m 
                                4 2 3 1] # z
const PM_TRANSFORM_TABLE2 = Int[0 0 0 0; # I 
                                0 0 4 0; # p
                                0 4 0 0; # m
                                0 0 0 0] # z

const PM_COEFF_TABLE1 = crationalize.([ 1.0  1.0  1.0  1.0;
                                        1.0  0.0  0.5 -1.0;
                                        1.0  0.5  0.0  1.0;
                                        1.0  1.0 -1.0  1.0])
const PM_COEFF_TABLE2 = crationalize.([ 0.0  0.0  0.0  0.0
                                        0.0  0.0  0.5  0.0
                                        0.0 -0.5  0.0  0.0
                                        0.0  0.0  0.0  0.0])
const PM_HOW_MANY = Int[1  1  1  1;
                        1  1  2  1;
                        1  2  1  1;
                        1  1  1  1] # I

@doc raw""" 
    QubitPM() -> OperatorSet

Creates the OperatorSet for a qubit using Raising and Lowering operators (``\sigma_+``, ``\sigma_-``, ``\sigma_z``, ``\sigma_I``).
"""
function QubitPM(symbol::String="")::OperatorSet
    ops = ["p", "m", "z"]
    ops_str = ["+", "-", "z"]
    base_pm = [[2], [3], [4]]
    pm_dag_inds = [[1], [3], [2], [4]]

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
                return symbol_str * str2sup(ops[ind-1]) * str2sub(sym) 
            else
                return symbol_str * "_" * sym 
            end
        else
            if formatted
                return ops[ind-1] * str2sub(sym) 
            else
                return ops[ind-1] * "_" * sym
            end
        end
    end
    function pm2latex(inds::Vector{Int}, sym::String)::String
        # create underscored string representation of sym using subscript_indexes
        ind = inds[1]
        curr_str::String = raw""
        if do_symbol
            curr_str *= symbol_latex*raw"^{" * ops_str[ind-1] * "}"
            if length(sym) > 0
                curr_str *= raw"_{" * sym * "}"
            end
        else
            curr_str *= raw"\hat{" * ops[ind-1] * "}"
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
    return OperatorSet("PM Qubit", "Fermion", 1, [1], base_pm, ops, pm_product, pm_dag, pm2str, pm2latex;
                      commutes=pmcommutes, min_ints=Int[1], max_ints=Int[4])
end
# Test 
#q = QubitPM()
#q.strs2ind("p")
#display(latexstring(q.op2latex(1, "i")))
#q.op_product(1, 2)


####################################### Helpers to transform Cayley Tables for new ordering: ##################################################################################################################
# apply the shift to a whole table of operator indices
#function shift_table_indices(tbl::Matrix{Int})
#    return [shift_index(x) for x in tbl]
#end
#
# permute matrix entries first by column [4,1,2,3], then rows the same way
#function permute_table(tbl::Matrix)
#    order = [4, 1, 2, 3]
#    return tbl[order, order]
#end
#function shiftandpermute(A::Matrix{Int})
#    return permute_table(shift_table_indices(A))
#end
