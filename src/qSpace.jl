module qSpace

using ComplexRationals
using ..FFunctions
using ..StringUtils

export OperatorSet, SubSpace, Parameter, OperatorType, StateSpace, string2operator_type, GLOBAL_STATE_SPACE

Is = Union{Int,Vector{Int}}
"""
    OperatorSet(name::String, fermion::Bool, len::Int, neutral_element::Union{Int,Vector{Int}}, base_ops::Union{Vector{Int},Vector{Vector{Int}}}, ops::Vector{String}, op_product::Function, op_dag::Function, strs2ind::Function, op2str::Function, op2latex::Function)

OperatorSets define the algebraic structure of a quantum system, defining ways to multiply and conjugate operators within the space, how to print them (both for plain and latex formatting), how to extract operators from strings.
We provide a few standard operator sets, such as QubitPauli, QubitPM and Ladder.
"""
struct OperatorSet
    name::String
    fermion::Bool           # fermions have symbol aftter operator name, bosons before
    len::Int                # length of indexes describing operator
    neutral_element::Union{Int,Vector{Int}}   # neutral element of the operator set
    base_ops::Union{Vector{Int},Vector{Vector{Int}}}
    non_base_ops::Dict{String, Vector{Tuple{ComplexRational, Is}}}
    ops::Vector{String}     # operator symbols
    op_product::Function    # takes operator indexes of two operators of this set and outputs a vector of tuples of coefficients and associated indexes for the resulting operators in this set
    op_dag::Function        # Create Complex Transpoose Conjugate
    op2str::Function        # transforms an operator index into a string for console printing
    op2latex::Function      # transforms an operator index into a LaTeX string for formatted LaTeXStrings
end
function Base.show(io::IO, os::OperatorSet)
    type_str = os.fermion ? "Fermionic" : "Bosonic"
    op_str = ""
    for i in 1:length(os.ops)
        op_str *= os.op2str(i, "p")
        if i == os.neutral_element
            op_str *= " (identity)"
        end
        if i < length(os.ops)
            op_str *= ", "
        end
    end
    print(io, os.name, " ($type_str):  " * op_str)
end

# An auxiliary struct representing one subspace of the Hilbert space.
""" 
    SubSpace(key::String, keys::Vector{String}, statespace_main_ind::Int, statespace_inds::Vector{Int}, op_set::OperatorSet, continuum::Bool, fermion::Bool)

SubSpace defines a subspace of a Hilbert space. It contains an operator set, aswell as additional information to reference and work with a subspace. 
Fermionic subspaces support multiple copies of the same subspace, so as to support continuous generalisations of the subspace. 
"""
struct SubSpace
    key::String                     # Original input key
    keys::Vector{String}            # Allowed keys for this subspace 
    statespace_main_ind::Int        # Which Vector to use for statespace_inds  (this is for accessing the string elements)
    statespace_inds::Vector{Int}    # Indices to access operator values in the corresponding statespace main ind  (this is for accessing the string elements)
    op_index_inds::Vector{Int}
    op_set::OperatorSet             # The operator set for this subspace.
    continuum::Bool
    fermion::Bool
end
# Define the custom show for SubSpace.
function Base.show(io::IO, statespace::SubSpace)
    # Print the subspace key and allowed keys.
    print(io, "SubSpace ", statespace.keys, ": ")
    # Use the OperatorSet's show for the op_set field.
    show(io, statespace.op_set)
end

""" 
    Parameter(var_name::String, var_of_t::Bool, var_of_continuum::Bool, var_continuum_index::Int=0; var_val::Union{Nothing,Number,Vector{Number},Function}=nothing, var_suffix::String="")

Parameter is a struct that represents a variable in the state space, and information of how to access and print it.
"""
mutable struct Parameter
    var_name::String
    var_str::String
    var_latex::String
    var_val::Union{Nothing,Number,Vector{Number},Function}
    var_of_t::Bool
    var_of_continuum::Bool
    var_continuum_index::Int
    var_subspace_indexes::Vector{Int}

    function Parameter(var_name::String, var_of_t::Bool, var_of_continuum::Bool, var_continuum_index::Int=0; var_val::Union{Nothing,Number,Vector{Number},Function}=nothing, var_suffix::String="", var_subspace_indexes::Vector{Int}=Int[])
        var_suffix_str = str2sub(var_suffix)
        var_str = ""
        if haskey(var_substitution, var_name)
            var_str *= var_substitution[var_name] * var_suffix_str
        else
            var_str *= var_name * var_suffix_str
        end
        var_suffix_latex = ""
        if length(var_suffix) > 0
            var_suffix_latex *= "_{" * var_suffix * "}"
        end
        var_latex = ""
        if haskey(var_substitution_latex, var_name)
            var_latex = var_substitution_latex[var_name] * var_suffix_latex
        else
            var_latex = var_name * var_suffix_latex
        end
        if var_of_t
            var_str *= "(t)"
            var_latex *= "(t)"
        end
        var_name_mod = var_name
        if length(var_suffix) > 0
            var_name_mod *= "_" * var_suffix
        end
        new(var_name_mod, var_str, var_latex, var_val, var_of_t, var_of_continuum, var_continuum_index, var_subspace_indexes)
    end
end

function Base.show(io::IO, p::Parameter)
    print(io, "par: ", p.var_str)
end

""" 
    OperatorType(name::String, hermitian::Bool=false, unitary::Bool=false, subspaces::Union{Nothing, Vector{Bool}}=nothing)


Define an operator Type, by declaring its name and properties such as hermitian and unitary. 
"""
struct OperatorType
    name::String 
    hermitian::Bool 
    unitary::Bool 
    subspaces::Vector{Bool}   # subspace groups
    non_subspaces::Vector{Bool} 
    non_trivial_op_indices::Vector{Bool} # indices of non-trivial components in the op_indices pciture for bath subsystems
    function OperatorType(name::String; hermitian::Bool=false, unitary::Bool=false, subspaces::Vector{Bool}, non_trivial_op_indices::Vector{Bool})
        if length(name) == 0
            error("OperatorType name cannot be empty.")
        end
        non_subspaces = .!subspaces
        new(name, hermitian, unitary, subspaces, non_subspaces, non_trivial_op_indices)
    end
end
function string2operator_type(s::String, subspace_keys::Vector{String}, dim::Int, subspace_vec::Vector{SubSpace})
    # use name(U,H) notation, separate the name from the brace. and parse the , separated elements within the brace to determine which properties are true (by default all are false) 
    hermitian = false 
    unitary = false 
    subspaces::Vector{Bool} = [true for _ in subspace_keys]
    name::String = ""
    modified_subspace_keys::Vector{String} = ["!"*key for key in subspace_keys]
    if occursin("(", s)  # Check if there is a brace in the string
        name, brace = split(s, "(")  # Output: ["As", "U,H"]
        if length(brace) > 0
            brace = split(brace, ")")[1]  # Output: "U,H"
            tokens = split(brace, ",")  # Output: ["U", "H"]
            tokens = [strip(t) for t in tokens]
            hermitian = "H" in tokens
            unitary = "U" in tokens 
            token_list = vcat(["H", "U"], subspace_keys, modified_subspace_keys)
            # if any token not in token_list, return error
            for token in tokens
                if !(token in token_list)
                    throw(ArgumentError("Invalid token: $token"))
                end
            end
            any_subspace_keys = false
            for s in subspace_keys
                if s in tokens
                    any_subspace_keys = true
                    break 
                end
            end
            if any_subspace_keys
                # check if any modified_subspace_keys in tokens
                for s in modified_subspace_keys
                    if s in tokens
                        throw(ArgumentError("Cannot have negatved subspace and subspace keys at the same time"))
                    end
                end
                for (i, s) in enumerate(subspace_keys)
                    subspaces[i] = s in tokens
                end
            else
                subspaces = fill(true, length(subspace_keys))
                for (i, s) in enumerate(modified_subspace_keys)
                    if s in tokens
                        subspaces[i] = false
                    end
                end
            end
        else
            throw(ArgumentError("Invalid Operator string: $s"))
        end
    else
        name = s
    end
    non_trivial_op_indices::Vector{Bool} = fill(false, dim)
    for (s_bool, subspace) in zip(subspaces, subspace_vec)
        if s_bool == true
            for ind in subspace.op_index_inds
                non_trivial_op_indices[ind] = true
            end
        end
    end
    return OperatorType(name, hermitian=hermitian, unitary=unitary, subspaces=subspaces, non_trivial_op_indices=non_trivial_op_indices)
end
function operator_type2string(p::OperatorType)
    curr_str = p.name 
    if any([p.hermitian, p.unitary])
        curr_str *= "("
        if p.hermitian
            curr_str *= "H,"
        end
        if p.unitary 
            curr_str *= "U,"
        end
        curr_str = curr_str[1:end-1] * ")"
    end
    return curr_str
end
function Base.show(io::IO, p::OperatorType)
    print(io, "Op: ", operator_type2string(p))
end

function operatertypes2commutator_matrix(optypes::Vector{OperatorType})
    # matrix of operatortypes commutation, true ==> the two commutators commute
    mat = Matrix{Bool}(undef, length(optypes), length(optypes))
    for i in 1:length(optypes)
        mat[i,i] = true
        for j in i+1:length(optypes)
            # check for any collisions of the operator_type subspaces 
            mat[i,j] = !any(optypes[i].subspaces .& optypes[j].subspaces)
            mat[j,i] = mat[i,j]
        end
    end
    return mat
end

"""
    GLOBAL_STATE_SPACE

Can define a globally accessible StateSpace
"""
GLOBAL_STATE_SPACE = nothing

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
    # Parameter fields:
    vars::Vector{Parameter}
    vars_str::Vector{String}
    vars_cont::Vector{Tuple{Vector{Int},Tuple}}    # For continuum variables: (subspace indices, standardized tuple for distribution) -> distribution is held here 
    how_many_by_continuum::Dict{Int,Int}   # for continuum subspaces, how many variables do we have 
    where_by_continuum::Dict{Int,Vector{Vector{Int}}}   # for continuum subspaces, where are our variables  
    where_by_continuum_var::Vector{Vector{Vector{Int}}}
    where_by_time::Vector{Int}
    where_const::Vector{Int}
    # Subspace definitions:
    subspaces::Vector{SubSpace}
    where_continuum::Vector{Int}
    continuum_indexes::Vector{Vector{Int}}
    fermionic_keys::Vector{String}
    bosonic_keys::Vector{String}
    neutral_op::Vector{Is}
    neutral_continuum_op::Vector{Vector{Is}}
    subspace_by_ind::Vector{Int}
    fone::FAtom
    # Abstract operators
    operatortypes::Vector{OperatorType}
    operator_names::Vector{String}
    operatortypes_commutator_mat::Matrix{Bool}
    function StateSpace(args...; operators::Union{String, Vector{String}}="A", make_global::Bool=true, kwargs...)
        subspaces = Vector{SubSpace}()
        fermionic_keys = String[]
        bosonic_keys = String[]
        used_strings = Set{String}()
        op_index_ind_max = 1
        for (key, val) in kwargs
            key_str = String(key)
            n = 1
            if isa(val, Tuple)
                n = val[1]
                op_set = val[2]
                if !(n isa Integer) || n ≤ 0
                    error("For subspace pattern key \"$key_str\", n must be a positive integer")
                end
            else
                n = 1
                op_set = val
            end

            if !isa(op_set, OperatorSet)
                error("Value for subspace key \"$key_str\" must be an OperatorSet")
            end
            continuum = false
            statespace_main_ind = 2 + !op_set.fermion
            # Check for a pattern key: all characters are consecutive in the alphabet.
            if n > 1
                continuum = true
                # must be fermionic 
                #if !op_set.fermion       # Outdated? Check if this works!!!
                #    error("Baths (n > 1) must be fermionic operators.")
                #end
                # Generate indices starting from the first character of key_str.
                if length(key_str) > 1
                    error("Pattern key ($key_str) must be a single character for spin baths.")
                end
                if key_str == "t"
                    error("Letter t is reserved for time, cannot be used for baths. ")
                end
                start_char = key_str[1]
                keys = [string(Char(start_char) + i) for i in 0:(n-1)]
                indices = [length(fermionic_keys) + i for i in 1:n]
                # Ensure no letter is reused across subspaces.
            else
                # Here, the index list is just the key itself.
                if op_set.fermion
                    index = length(fermionic_keys) + 1
                else
                    index = length(bosonic_keys) + 1
                end
                keys = [key_str]
                indices = [index]
            end
            if continuum && !op_set.fermion
                error("Cannot Create non fermionic Continuum Subspaces!")
            end
            op_index_inds = collect(op_index_ind_max:op_index_ind_max+n-1)
            op_index_ind_max += n
            push!(subspaces, SubSpace(key_str, keys, statespace_main_ind, indices, op_index_inds, op_set, continuum, op_set.fermion))
            if op_set.fermion
                append!(fermionic_keys, keys)
            else
                append!(bosonic_keys, keys)
            end
            for c in keys
                if c in used_strings
                    error("Duplicate string ($c) found in subspace definitions")
                end
                push!(used_strings, c)
            end
        end

        vars::Vector{Parameter} = []
        vars_cont::Vector{Tuple{Vector{Int},Tuple}} = []
        where_by_continuum::Dict{Int,Vector{Vector{Int}}} = Dict()
        how_many_by_continuum::Dict{Int,Int} = Dict()
        where_by_time::Vector{Int} = []
        where_const::Vector{Int} = []
        for i in 1:length(subspaces)
            if subspaces[i].continuum
                where_by_continuum[i] = Vector[]
                how_many_by_continuum[i] = 0
            end
        end
        # Convert state variables (positional arguments) into strings.
        # Convert state variables (positional arguments) into Parameters.
        # first sort args by length -> this helps avoid conflicts in potential string parsing, so that alpha' and alpha are both deteccted properly for example, since we search first for alpha'...
        args = sort(collect(args), by=length, rev=true)
        for arg in args
            if !isa(arg, String)
                error("Parameters must be declared via Strings, problem with $arg.")
            end

            var_of_t = occursin("(t)", arg)
            raw_arg = replace(arg, "(t)" => "")
            var_of_continuum = occursin("_", raw_arg)

            pre, sub = "", ""
            var_continuum_index = 0

            if var_of_continuum
                parts = split(raw_arg, "_")
                if length(parts) != 2
                    error("Malformed continuum variable: $arg")
                end
                pre, sub = parts
                if length(sub) > 1
                    error("Continuum index ($sub) in ($arg) must be a single character.")
                end

                if pre in used_strings
                    error("Duplicate variable name ($pre) found.")
                end
                push!(used_strings, pre)

                curr_subs::Vector{String} = []
                for (index, (key, val)) in enumerate(kwargs)
                    key_str = String(key)
                    if sub[1] == key_str[1]
                        var_of_continuum = true
                        var_continuum_index = index
                        curr_subs = copy(subspaces[var_continuum_index].keys)
                        break
                    end
                end
                if var_continuum_index == 0
                    error("Couldn't find a continuum with key ($sub) for variable ($arg)")
                end

                how_many_by_continuum[var_continuum_index] += 1
                i_vec::Vector{Int} = []
                curr_where = Int[]
                for sublabel in curr_subs
                    p = Parameter(string(pre), var_of_t, true, var_continuum_index; var_suffix=sublabel)
                    push!(vars, p)
                    push!(i_vec, length(vars))
                    push!(curr_where, length(vars))
                    if var_of_t
                        push!(where_by_time, length(vars))
                    else
                        push!(where_const, length(vars))
                    end
                end
                push!(where_by_continuum[var_continuum_index], curr_where)
                push!(vars_cont, (i_vec, ()))

            else
                pre = raw_arg
                if pre in used_strings
                    error("Duplicate variable name ($pre) found.")
                end
                push!(used_strings, pre)
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
        vars_str::Vector{String} = [p.var_str for p in vars]

        neutral_op = [s.op_set.neutral_element for s in subspaces for key in s.keys]
        fone = FAtom(ComplexRational(1,0,1), zeros(Int, length(vars)))
        if !isa(operators, Vector) 
            operators = [operators]
        end
        operatortypes::Vector{OperatorType} = []
        subspace_keys::Vector{String} = [s.key for s in subspaces]
        for s in operators 
            push!(operatortypes, string2operator_type(s, subspace_keys, length(neutral_op), subspaces))
        end
        operator_names::Vector{String} = [ot.name for ot in operatortypes]
        operatortypes_commutator_mat = operatertypes2commutator_matrix(operatortypes)

        where_continuum::Vector{Int} = sort(collect(keys(how_many_by_continuum))) 
        continuum_indexes::Vector{Vector{Int}} = [copy(subspaces[i].statespace_inds) for i in where_continuum]
        neutral_continuum_op::Vector{Vector{Is}} = [[neutral_op[i] for i in cont] for cont in continuum_indexes]
        where_by_continuum_var::Vector{Vector{Vector{Int}}} = [where_by_continuum[i] for i in where_continuum]
        subspace_by_ind::Vector{Int} = zeros(Int, length(neutral_op))
        for (i, s) in zip(where_continuum, continuum_indexes)
            subspace_by_ind[s] .= i
        end
        qss = new(vars, vars_str, vars_cont, how_many_by_continuum, where_by_continuum, where_by_continuum_var, where_by_time, where_const, subspaces, where_continuum, continuum_indexes, fermionic_keys, bosonic_keys, neutral_op, neutral_continuum_op, subspace_by_ind, fone, operatortypes, operator_names, operatortypes_commutator_mat)
        if make_global
            global GLOBAL_STATE_SPACE = qss
        end
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

include("OperatorSets/Qubit_Pauli.jl")
include("OperatorSets/Qubit_PM.jl")
include("OperatorSets/Ladder.jl")

end # module qSpace