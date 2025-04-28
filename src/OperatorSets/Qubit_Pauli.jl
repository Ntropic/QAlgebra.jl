export QubitPauli

#### Cayley tables for Pauli Operators
const PAULI_TRANSFORM_TABLE = [
    4 3 2 1;
    3 4 1 2;
    2 1 4 3;
    1 2 3 4]
const PAULI_COEFF_TABLE = [
    1 im -im 1;
    -im 1 im 1;
    im -im 1 1;
    1 1 1 1]

"QubitPauli returns an OperatorSet for qubit (Pauli) operators."
"""
    QubitPauli() -> OperatorSet

Creates the OperatorSet for a qubit using Pauli operators ($\sigma_x$, $\sigma_y$, $\sigma_z$, $\sigma_I$).
"""
function QubitPauli()
    ops = ["x", "y", "z", "I"]
    base_pauli = [1, 2, 3]
    # Define the transformation function using the PAULI_* tables.
    function pauli_product(op1::Int, op2::Int)::Vector{Tuple{Complex,Int}}
        # Look up the coefficient and new operator index.
        c = PAULI_COEFF_TABLE[op1, op2]
        new_index = PAULI_TRANSFORM_TABLE[op1, op2]
        # Return as a one-term sum.
        return [(c, new_index)]
    end
    function pauli_dag(op::Int)::Vector{Tuple{Complex,Int}}
        return [(1.0, op)]
    end
    function paulistr2ind(str::String)::Vector{Tuple{Complex,Int}}
        if length(str) == 0
            return [(1.0 + 0im, 4)]
        end
        res = findfirst(==(str), ops)
        if res == nothing
            error("Invalid Pauli string: $str, must be one of $ops, or empty")
        end
        return [(1.0 + 0im, res)]
    end
    function paulistr2ind(str::Vector{String})::Vector{Tuple{Complex,Int}}
        # iteratively take paulistr2ind for each element and unify them iteratively using pauli_product 
        if length(str) == 0
            return [(1.0 + 0im, 4)]
        end
        inds = [paulistr2ind(s)[1][2] for s in str]
        coeff = 1.0 + 0im
        while length(inds) > 1
            new_coeff, new_ind = pauli_product(inds[end-1], inds[end])[1]
            coeff *= new_coeff
            inds[end-1] = new_ind
            deleteat!(inds, length(inds))
        end
        return [(coeff, inds[1])]
    end
    function pauli2str(ind::Int, sym::String="")::String
        # create underscored string representation of sym using subscript_indexes
        return ops[ind] * str2sub(sym)
    end
    function pauli2latex(ind::Int, sym::String, do_sigma::Bool=false)::String
        # create underscored string representation of sym using subscript_indexes
        curr_str::String = raw""
        if do_sigma
            curr_str *= raw"\hat{\sigma}_{" * ops[ind] * "}"
            curr_str *= raw"^{(" * sym * ")}"
        else
            curr_str *= raw"\hat{" * ops[ind] * "}"
            curr_str *= raw"_{" * sym * "}"
        end
        return curr_str
    end
    return OperatorSet("Pauli Qubit", true, 1, 4, base_pauli, ops, pauli_product, pauli_dag, paulistr2ind, pauli2str, pauli2latex)
end
## Test 
#q = QubitPauli()
#q.strs2ind("x")
#display(latexstring(q.op2latex(1, "i")))
#q.op_product(1, 2)
#q