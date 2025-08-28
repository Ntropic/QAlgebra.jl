module QSpace

using ComplexRationals
using ..CFunctions
using ..StringUtils

export OperatorSet
export SubSpace, SubSpaceDefinitions, SubSpaceInfo, SubSpaceIndex, outer, inner, expanded, Index2Symbol, Index2String
export OperatorType, OperatorTypeInfo, OperatorDefinitions
export Parameter, ParameterInfo, ParameterDefinitions, map_by_subspace, map_by_tindex
export StateSpace

Is = Vector{Int}
"""
    OperatorSet(name::String, fermion::Bool, len::Int, neutral_element::Union{Int,Vector{Int}}, base_ops::Union{Vector{Int},Vector{Vector{Int}}}, ops::Vector{String}, op_product::Function, op_dag::Function, strs2ind::Function, op2str::Function, op2latex::Function)

OperatorSets define the algebraic structure of a quantum system, defining ways to multiply and conjugate operators within the space, how to print them (both for plain and latex formatting), how to extract operators from strings.
We provide a few standard operator sets, such as QubitPauli, QubitPM and Ladder.
"""
struct OperatorSet
    name::String
    particle_type::String   # fermion, boson, anyon...
    len::Int                # length of indexes describing operator
    neutral_element::Vector{Int}   # neutral element of the operator set
    base_ops::Vector{Vector{Int}}
    ops::Vector{String}     # operator symbols
    op_product::Function    # takes operator indexes of two operators of this set and outputs a vector of tuples of coefficients and associated indexes for the resulting operators in this set
    op_dag::Function        # Create Complex Transpoose Conjugate
    op2str::Function        # transforms an operator index into a string for console printing
    op2latex::Function      # transforms an operator index into a LaTeX string for formatted LaTeXStrings
    commutes::Function
    function OperatorSet(name::String, particle_type::String, len::Int, neutral_element::Vector{Int}, base_ops::Vector{Vector{Int}}, ops::Vector{String}, op_product::Function, op_dag::Function, op2str::Function, op2latex::Function, commutes::Function)
        return new(name, particle_type, len, neutral_element, base_ops, ops, op_product, op_dag, op2str, op2latex, commutes)
    end
    function OperatorSet(name::String, particle_type::String, len::Int, neutral_element::Vector{Int}, base_ops::Vector{Vector{Int}}, ops::Vector{String}, op_product::Function, op_dag::Function, op2str::Function, op2latex::Function)
        function commutes(op1::Vector{Int}, op2::Vector{Int}) # multiply to test commute => probably much slower than a custom implementation
            if op1 == op2 || op1 == neutral_element || op2 == neutral_element
                return true
            end
            prod_1 = op_product(op1, op2) # isa Vector{Tuple{ComplexRational,Vector{Int}}}
            prod_2 = op_product(op2, op1) # isa Vector{Tuple{ComplexRational,Vector{Int}}}
            # sort prod1 and prod2
            if length(prod_1) != length(prod_2)
                return false
            end
            sort!(prod_1, by=x -> x[2])
            sort!(prod_2, by=x -> x[2])
            for k in eachindex(prod_1)
                if prod_1[k][2] != prod_2[k][2] || prod_1[k][1] != -prod_2[k][1]
                    return false
                end
            end
            return true
        end
        return new(name, particle_type, len, neutral_element, base_ops, ops, op_product, op_dag, op2str, op2latex, commutes)
    end
    function OperatorSet() # Dummy Operator Set 
        dummy_fun(args...; kwargs...) = error("OperatorSet not initialized")
        return OperatorSet("Unspecified", "none", 1, Int[0],
                        Vector{Vector{Int}}(),
                        Dict{String, Vector{Tuple{ComplexRational, Vector{Int}}}}(),
                        String[],
                        dummy_fun, dummy_fun, dummy_fun, dummy_fun, dummy_fun)
    end
end
#function OperatorSet
function Base.show(io::IO, os::OperatorSet)
    op_str = ""
    for i in 1:length(os.ops)
        for j in 1:os.len
            # place the i in j'th position 
            z = zeros(Int, os.len)
            z[j] = i
            op_str *= os.op2str(z, "p")
        end
        if i == os.neutral_element
            op_str *= " (identity)"
        end
        if i < length(os.ops)
            op_str *= ", "
        end
    end
    print(io, os.name, " (", os.particle_type, "):  " * op_str)
end
include("OperatorSets/Qubit_Pauli.jl")
include("OperatorSets/Qubit_PM.jl")
include("OperatorSets/Ladder.jl")

include("QSpaceOps/QSpace_subspaces.jl")
include("QSpaceOps/QSpace_abstract.jl")
include("QSpaceOps/QSpace_parameters.jl")

"""
    StateSpace(subspace_def::SubSpaceDefinitions, op_def::OperatorDefinitions, param_def::ParameterDefinitions; max_t_ind::Int=0) -> StateSpace

Constructs a combined Hilbert and Parameter space. The Hilbert space consists of different subspaces, themselves composed of different operator sets.
The Parameter space also defines the variables, that are needed to describe equations on the Hilbert space and abstract operators, that are not yet specified. 
Optionally you can also allow for multiple time dimensions, which can be useful for solving nested integrals over different time parameters.
"""
struct StateSpace
    # Subspace definitions:
    subspaces::Vector{SubSpace}
    subspace_info::SubSpaceInfo    # Info object containing references to all the indexing of outer and inner subspaces

    # Abstract operators
    operatortypes::Vector{OperatorType}
    operatortype_info::OperatorTypeInfo

    # Parameter fields:
    vars::Vector{Parameter}
    param_info::ParameterInfo

    I_op::Vector{Is}               # Neutral Vector of all expanded subspaces
    I_ensemble_op::Vector{Vector{Is}}      # Neutral Vector of all expanded ensemble subspaces
    c_one::CAtom                            # onelike function in CFunctions 
    c_zero::CAtom                           # zerolike function in CFunctions 

    function StateSpace(subspace_def::SubSpaceDefinitions, op_def::OperatorDefinitions, param_def::ParameterDefinitions; max_t_ind::Int=0)
        # ==========> 1st Subspaces <==========
        subspaces = subspace_def.subspaces
        subspace_info = SubSpaceInfo(subspaces)
        used_symbols = subspace_def.used_symbols
        I_op, I_ensemble_op = subspace_def.I_op, subspace_def.I_ensemble_op  # These are

        # ==========> 2nd Abstract Operators <==========
        operatortypes = OperatorDefinitions2OperatorType(op_def, subspace_def)
        operatortype_info = OperatorTypeInfo(operatortypes, commute_fun=op_def.commute_fun, check_n=op_def.check_n) 

        # ==========> 3rd Parameters <==========
        vars, param_info = ParameterDefinitions2Parameters(param_def, subspace_info, used_symbols, max_t_ind)
    
        # Generate the string representations
        c_one = CAtom(ComplexRational(1,0,1), zeros(Int, length(vars)))
        c_zero = CAtom(ComplexRational(0,0,1), zeros(Int, length(vars)))
        qss = new( subspaces, subspace_info,                                      # Subspaces
                operatortypes, operatortype_info,                                 # Abstract Operators 
                vars, param_info,                                                 # Variables / Parameters
                I_op, I_ensemble_op, c_one, c_zero)                         # Pecomputed operator blueprints 
        return qss
    end
end
# Define the custom show for StateSpace.
function Base.show(io::IO, statespace::StateSpace)
    # First line: StateSpace and its variables.
    var_str = join([p.var_str for p in statespace.vars], ", ")
    println(io, "StateSpace: [" * var_str * "]")
    # Then print each subspace on its own line.
    for ss in statespace.subspaces
        println(io, "   - ", string(ss))
    end
    for op in statespace.operatortypes
        println(io, "   - ", string(op))
    end
end

## Test 
#xi, yi, zi = base_operators("i", qs)
#I = base_operators("I", qs)
#alpha, beta = base_operators("vars", qs)
function cleanup_terms(terms::Vector{Tuple{T,S}})::Vector{Tuple{T,S}} where {T<:Number,S}
    # 1) sort once by index
    sort!(terms, by = x -> x[2])
    # 2) prealloc output to worst‑case length and scan in one pass
    n = length(terms)
    T0 = typeof(terms[1][1])
    S0 = typeof(terms[1][2])
    cleaned = Vector{Tuple{T0,S0}}(undef, n)
    cnt = 0
    i = 1
    @inbounds while i ≤ n
        sumc, idx = terms[i]           # destructure once
        j = i + 1
        # inner loop: accumulate identical idx
        @inbounds while j ≤ n && terms[j][2] == idx
            sumc += terms[j][1]
            j += 1
        end

        # push nonzero
        if sumc != zero(T0)
            cnt += 1
            cleaned[cnt] = (sumc, idx)
        end

        i = j
    end
    resize!(cleaned, cnt)                # trim unused slots
    return cleaned
end


end # module QSpace