module QSpace

using ComplexRationals
using ..CFunctions
using ..StringUtils

export OperatorSet, SubSpace, Parameter, OperatorType, StateSpace, string2operator_type

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
    non_base_ops::Dict{String, Vector{Tuple{ComplexRational, Is}}}
    ops::Vector{String}     # operator symbols
    op_product::Function    # takes operator indexes of two operators of this set and outputs a vector of tuples of coefficients and associated indexes for the resulting operators in this set
    op_dag::Function        # Create Complex Transpoose Conjugate
    op2str::Function        # transforms an operator index into a string for console printing
    op2latex::Function      # transforms an operator index into a LaTeX string for formatted LaTeXStrings
    commutes::Function
    function OperatorSet(name::String, particle_type::String, len::Int, neutral_element::Vector{Int}, base_ops::Vector{Vector{Int}}, non_base_ops::Dict{String, Vector{Tuple{ComplexRational, Is}}}, ops::Vector{String}, op_product::Function, op_dag::Function, op2str::Function, op2latex::Function, commutes::Function)
        return new(name, particle_type, len, neutral_element, base_ops, non_base_ops, ops, op_product, op_dag, op2str, op2latex, commutes)
    end
    function OperatorSet(name::String, particle_type::String, len::Int, neutral_element::Vector{Int}, base_ops::Vector{Vector{Int}}, non_base_ops::Dict{String, Vector{Tuple{ComplexRational, Is}}}, ops::Vector{String}, op_product::Function, op_dag::Function, op2str::Function, op2latex::Function)
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
        return new(name, particle_type, len, neutral_element, base_ops, non_base_ops, ops, op_product, op_dag, op2str, op2latex, commutes)
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
    print(io, os.name, " ($os.particle_type):  " * op_str)
end
include("OperatorSets/Qubit_Pauli.jl")
include("OperatorSets/Qubit_PM.jl")
include("OperatorSets/Ladder.jl")

include("QSpaceOps/QSpace_subspaces.jl")
include("QSpaceOps/QSpace_abstract.jl")
include("QSpaceOps/QSpace_parameters.jl")


struct VarDefinitions
    args::Vector{String}
    args_symbols::Vector{Symbol} = []  
    ts::Vector{Int} # time dependent variables
    function VarDefinitions(args...; ts::Vector{Int} = Int[]) 
        new_args = String[]
        args_symbols::Vector{Symbol} = []  
        for arg in args 
            if arg isa Symbol 
                push!(new_args, string(arg))
                push!(args_symbols, arg)  # Store symbols
            else
                push!(new_args, string(arg)) 
                push!(args_symbols, Symbol(arg))  # Convert
            end
        end
        return new(new_args, args_symbols, ts)  
    end 
end

"""
    StateSpace(args...; kwargs...)

Constructs a combined Hilbert and Parameter space. The Hilbert space consists of different subspaces, themselves composed of different operator sets. The Parameter space defines the variables, that are needed to describe equations on the Hilbert space.
    - **args**: A variable number of symbols or strings representing the state variables. Can refer to indexes of subsystems via for underscore notation, i.e., "alpha_i" or declare time dependence via for example "alpha(t)".
    - **kwargs**: Each keyword is interpreted as a subspace label. The values are either Operator Sets or Tuples with an integer and an OperatorSet. The integer is the number of indexes generated for the subspace.
    - **operators**: Specifies the types of abstract Operators from Strings such as "A(U,H)", which would result in an operator with name "A" and properties Hermitian (H) and Unitary (U).
                    You can specify the subspaces on which the operator acts non trivially via the subspace keys (kwarg keys). Alternatively, you can specify the subspaces on which it doesn#t act non trivially via the negatved subspace keys (!kwarg keys).
                    If no subspace is specified, it defaults to all subspaces.
                    Example: "A(U,H,i)" constructs an abstract operator "A" that is both hermitian and unitary and acts non trivially on the subspace "i".
                    Example: "B(H,!i)" constructs an abstract operator "B" that is hermitian and acts trivially on the subspace "i".
                    Example: "C(U)" constructs an abstract operator "C" that is unitary and acts non-trivially on all subspaces.
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
    how_many_by_ensemble::Dict{Int,Int}   # for ensemble subspaces, how many variables do we have 
    where_by_ensemble::Dict{Int,Vector{Vector{Int}}}   # for ensemble subspaces, where are our variables  
    where_by_ensemble_var::Vector{Vector{Vector{Int}}}
    where_by_time::Vector{Int}
    where_const::Vector{Int}
    
    subspace_by_ind::Vector{Int}
    c_one::CAtom

    function StateSpace(subspace_def::SubSpaceDefinitions, op_def::OpDefinitions) 
        # ==========> 1st Subspaces <==========
        subspaces = subspace_def.subspaces
        subspace_info = SubSpaceInfo(subspaces)
        used_symbols = subspace_def.used_symbols

        # ==========> 2nd Abstract Operators <==========
        operatortypes = OpDefinitions2OperatorType(op_def, subspace_def)
        operatortype_info = OperatorTypeInfo(operatortypes, op_def.commute_fun, op_def.check_n) 

        # ==========> 3rd Variables <==========
        vars::Vector{Parameter} = []
        vars_cont::Vector{Tuple{Vector{Int},Tuple}} = []
        where_by_ensemble::Dict{Int,Vector{Vector{Int}}} = Dict()
        how_many_by_ensemble::Dict{Int,Int} = Dict()
        where_by_time::Vector{Int} = []
        where_const::Vector{Int} = []
        for i in 1:length(subspaces)
            if subspaces[i].ensemble
                where_by_ensemble[i] = Vector[]
                how_many_by_ensemble[i] = 0
            end
        end
        # Convert state variables (positional arguments) into strings.
        # Convert state variables (positional arguments) into Parameters.
        #args = sort(collect(args), by=length, rev=true)  # first sort args by length -> this helps avoid conflicts in potential string parsing, so that alpha' and alpha are both deteccted properly for example, since we search first for alpha'...
        for arg in args
            if !isa(arg, String)
                error("Parameters must be declared via Strings, problem with $arg.")
            end

            var_of_t = occursin("(t)", arg)
            raw_arg = replace(arg, "(t)" => "")
            var_of_ensemble = occursin("_", raw_arg)

            pre, sub = "", ""
            var_ensemble_index = 0

            if var_of_ensemble
                parts = split(raw_arg, "_")
                if length(parts) != 2
                    error("Malformed ensemble variable: $arg")
                end
                pre, sub = parts
                if length(sub) > 1
                    error("ensemble index ($sub) in ($arg) must be a single character.")
                end

                if pre in used_symbols
                    error("Duplicate variable name ($pre) found.")
                end
                push!(used_symbols, pre)

                curr_subs::Vector{String} = []
                for (index, (key, val)) in enumerate(kwargs)
                    key_str = String(key)
                    if sub[1] == key_str[1]
                        var_of_ensemble = true
                        var_ensemble_index = index
                        curr_subs = copy(subspaces[var_ensemble_index].keys)
                        break
                    end
                end
                if var_ensemble_index == 0
                    error("Couldn't find a ensemble with key ($sub) for variable ($arg)")
                end

                how_many_by_ensemble[var_ensemble_index] += 1
                i_vec::Vector{Int} = []
                curr_where = Int[]
                for sublabel in curr_subs
                    p = Parameter(string(pre), var_of_t, true, var_ensemble_index; var_suffix=sublabel)
                    push!(vars, p)
                    push!(i_vec, length(vars))
                    push!(curr_where, length(vars))
                    if var_of_t
                        push!(where_by_time, length(vars))
                    else
                        push!(where_const, length(vars))
                    end
                end
                push!(where_by_ensemble[var_ensemble_index], curr_where)
                push!(vars_cont, (i_vec, ()))

            else
                pre = raw_arg
                if pre in used_symbols
                    error("Duplicate variable name ($pre) found.")
                end
                push!(used_symbols, pre)
                p = Parameter(pre, var_of_t, false, 0)
                push!(vars, p)
                if var_of_t
                    push!(where_by_time, length(vars))
                else
                    push!(where_const, length(vars))
                end
            end
        end

        # Generate the string representations
        c_one = CAtom(ComplexRational(1,0,1), zeros(Int, length(vars)))
        
        where_by_ensemble_var::Vector{Vector{Vector{Int}}} = [where_by_ensemble[i] for i in where_ensemble]
        subspace_by_ind::Vector{Int} = zeros(Int, length(I_op))
        for (i, sub) in enumerate(subspaces)
            inds = sub.ss_inner_ind
            subspace_by_ind[inds] .= i
        end
        qss = new( subspaces, subspace_info,                                      # Subspaces
                operatortypes, operatortype_info,                                 # Abstract Operators 
                vars, how_many_by_ensemble, where_by_ensemble, where_by_ensemble_var, where_by_time, where_const, subspace_by_ind, c_one)   # Variables 
        return qss
    end
end
# Define the custom show for StateSpace.
function Base.show(io::IO, statespace::StateSpace)
    # First line: StateSpace and its variables.
    var_str = join([p.var_str_fun() for p in statespace.vars], ", ")
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