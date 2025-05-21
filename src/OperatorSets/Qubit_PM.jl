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

raw""" 
    QubitPM() -> OperatorSet

Creates the OperatorSet for a qubit using Raising and Lowering operators ($\sigma_+$, $\sigma_-$, $\sigma_z$, $\sigma_I$).
"""
function QubitPM()
    ops = ["p", "m", "z", "I"]
    ops_str = ["+", "-", "z", "I"]
    non_base_ops::Dict{String, Vector{Tuple{ComplexRational, Is}}} = Dict("x"=> [(ComplexRational(1,0,1), 1),(ComplexRational(1,0,1), 2)], "y" => [(ComplexRational(0,-1,1), 1), (ComplexRational(0,1,1), 2)])
    base_pm = [1, 2, 3]
    pm_dag_inds = [2, 1, 3, 4]
    # Define the transformation function using the PAULI_* tables.
    function pm_product(op1::Int, op2::Int)::Vector{Tuple{ComplexRational,Int}}
        # Look up the coefficient and new operator index.
        c1 = PM_COEFF_TABLE1[op1, op2]
        new_index1 = PM_TRANSFORM_TABLE1[op1, op2]
        if PM_HOW_MANY[op1, op2] == 2
            c2 = PM_COEFF_TABLE2[op1, op2]
            new_index2 = PM_TRANSFORM_TABLE2[op1, op2]
            return [(c1, new_index1), (c2, new_index2)]
        else
            # Return as a one-term sum.
            return [(c1, new_index1)]
        end
    end
    function pm_dag(op::Int)::Vector{Tuple{ComplexRational,Int}}
        return [(ComplexRational(1,0,1), pm_dag_inds[op])]
    end
    function pmstr2ind(str::String)::Vector{Tuple{ComplexRational,Int}}
        if length(str) == 0
            return [(ComplexRational(1,0,1), 4)]
        end
        res = findfirst(==(str), ops)
        if isnothing(res)
            error("Invalid Pauli string: $str, must be one of $ops, or empty")
        end
        return [(ComplexRational(1,0,1), res)]
    end
    function pmstr2ind(str::Vector{String})::Vector{Tuple{ComplexRational,Int}}
        # iteratively take paulistr2ind for each element and unify them iteratively using pauli_product 
        if length(str) == 0
            return [(ComplexRational(1,0,1), 4)]
        end
        inds = [pmstr2ind(s)[1][2] for s in str]
        coeff = ComplexRational(1,0,1)
        while length(inds) > 1
            new_coeff, new_ind = pm_product(inds[end-1], inds[end])[1]
            coeff *= new_coeff
            inds[end-1] = new_ind
            deleteat!(inds, length(inds))
        end
        return [(coeff, inds[1])]
    end
    function pm2str(ind::Int, sym::String="")::String
        # create underscored string representation of sym using subscript_indexes
        return ops[ind] * str2sub(sym)
    end
    function pm2latex(ind::Int, sym::String, do_sigma::Bool=false)::String
        # create underscored string representation of sym using subscript_indexes
        curr_str::String = raw""
        if do_sigma
            curr_str *= raw"\hat{\sigma}_{" * ops_str[ind] * "}"
            curr_str *= raw"^{(" * sym * ")}"
        else
            curr_str *= raw"\hat{" * ops[ind] * "}"
            curr_str *= raw"_{" * sym * "}"
        end
        return curr_str
    end
    return OperatorSet("PM Qubit", true, 1, 4, base_pm, non_base_ops, ops, pm_product, pm_dag, pmstr2ind, pm2str, pm2latex)
end
# Test 
#q = QubitPM()
#q.strs2ind("p")
#display(latexstring(q.op2latex(1, "i")))
#q.op_product(1, 2)