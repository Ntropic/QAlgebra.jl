module qSpace

using ..StringUtils

export OperatorSet, SubSpace, Parameter, StateSpace

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
    ops::Vector{String}     # operator symbols
    op_product::Function    # takes operator indexes of two operators of this set and outputs a vector of tuples of coefficients and associated indexes for the resulting operators in this set
    op_dag::Function        # Create Complex Transpoose Conjugate
    strs2ind::Function      # transform a vector of strings into the index representation
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
    statespace_main_ind::Int        # Which Vector to use for statespace_inds
    statespace_inds::Vector{Int}    # Indices to access operator values 
    op_set::OperatorSet             # The operator set for this subspace.
    continuum::Bool
    fermion::Bool
end
# Define the custom show for SubSpace.
function Base.show(io::IO, qspace::SubSpace)
    # Print the subspace key and allowed keys.
    print(io, "SubSpace ", qspace.keys, ": ")
    # Use the OperatorSet's show for the op_set field.
    show(io, qspace.op_set)
end

""" 
    Parameter(var_name::String, var_of_t::Bool, var_of_continuum::Bool, var_continuum_index::Int=0; var_val::Union{Nothing,Number,Vector{Number},Function}=nothing, var_suffix::String="")

Parameter is a struct that represents a variable in the state space, and information of how to access and print it.
"""
struct Parameter
    var_name::String
    var_str::String
    var_latex::String
    var_val::Union{Nothing,Number,Vector{Number},Function}
    var_of_t::Bool
    var_of_continuum::Bool
    var_continuum_index::Int

    function Parameter(var_name::String, var_of_t::Bool, var_of_continuum::Bool, var_continuum_index::Int=0; var_val::Union{Nothing,Number,Vector{Number},Function}=nothing, var_suffix::String="")
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
        new(var_name_mod, var_str, var_latex, var_val, var_of_t, var_of_continuum, var_continuum_index)
    end
end

function Base.show(io::IO, p::Parameter)
    print(io, "par: ", p.var_str)
end

global GLOBAL_STATE_SPACE = nothing

"""
    StateSpace(args...; kwargs...)

Constructs a combined Hilbert and Parameter space. The Hilbert space consists of different subspaces, themselves composed of different operator sets. The Parameter space defines the variables, that are needed to describe equations on the Hilbert space.
    - **args**: A variable number of symbols or strings representing the state variables. Can refer to indexes of subsystems via for underscore notation, i.e., "alpha_i" or declare time dependence via for example "alpha(t)".
    - **kwargs**: Each keyword is interpreted as a subspace label. The values are either Operator Sets or Tuples with an integer and an OperatorSet. The integer is the number of indexes generated for the subspace.
"""
struct StateSpace
    # Parameter fields:
    vars::Vector{Parameter}
    vars_str::Vector{String}
    vars_cont::Vector{Tuple{Vector{Int},Tuple}}    # For continuum variables: (subspace indices, standardized tuple)
    how_many_by_continuum::Dict{Int,Int}   # for continuum subspaces, how many variables do we have 
    where_by_continuum::Dict{Int,Vector{Vector{Int}}}   # for continuum subspaces, how many variables do we have 
    where_by_time::Vector{Int}
    where_const::Vector{Int}
    # Subspace definitions:
    subspaces::Vector{SubSpace}
    fermionic_keys::Vector{String}
    bosonic_keys::Vector{String}
    neutral_op::Vector{Is}

    function StateSpace(args...; make_global::Bool=true, kwargs...)
        subspaces = Vector{SubSpace}()
        fermionic_keys = String[]
        bosonic_keys = String[]
        used_strings = Set{String}()

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
                if !op_set.fermion
                    error("Baths (n > 1) must be fermionic operators.")
                end
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
            push!(subspaces, SubSpace(key_str, keys, statespace_main_ind, indices, op_set, continuum, op_set.fermion))
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
        for arg in args
            # Create a Parameter for each arg 
            if !isa(arg, String)
                error("Parameters are declared via Strings, problem for $arg.")
            end
            var_of_continuum = false
            var_of_t = false
            var_continuum_index = 0
            if occursin("(t)", arg)
                var_of_t = true
                pre, sub = split(arg, "(t)")
                pre = string(pre)
                push!(vars, Parameter(pre, var_of_t, var_of_continuum, var_continuum_index))
                push!(where_by_time, length(vars))
            elseif occursin("_", arg)
                pre, sub = split(arg, "_")
                pre = string(pre)
                sub = string(sub)
                if length(sub) > 1
                    error("Continuum index ($sub) in ($arg) must be a single character.")
                end
                if pre in used_strings
                    error("Duplicate string ($pre) found in variable definitions")
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
                how_many_by_continuum[var_continuum_index] += 1
                # add Parameter for each sub in curr_subs
                if var_continuum_index == 0
                    error("Couldn't find a continuum with key ($pre) for original input ($arg)")
                end
                i_vec::Vector{Int} = []
                curr_where = Int[]
                for (i, sub) in enumerate(curr_subs)
                    push!(vars, Parameter(pre, var_of_t, var_of_continuum, var_continuum_index; var_suffix=sub))
                    push!(curr_where, length(vars))
                    push!(i_vec, length(vars))
                end
                push!(where_by_continuum[var_continuum_index], curr_where)
                push!(vars_cont, (i_vec, ()))
            else
                push!(vars, Parameter(arg, var_of_t, var_of_continuum, var_continuum_index))
                push!(where_const, length(vars))
                if arg in used_strings
                    error("Duplicate string ($arg) found in variable definitions")
                end
                push!(used_strings, arg)
            end
        end
        vars_str::Vector{String} = []
        for p in vars
            push!(vars_str, p.var_name)
        end
        neutral_op = [s.op_set.neutral_element for s in subspaces for key in s.keys]
        qss = new(vars, vars_str, vars_cont, how_many_by_continuum, where_by_continuum, where_by_time, where_const, subspaces, fermionic_keys, bosonic_keys, neutral_op)
        if make_global
            global GLOBAL_STATE_SPACE = qss
        end
        return qss
    end
end
# Define the custom show for StateSpace.
function Base.show(io::IO, qspace::StateSpace)
    # First line: StateSpace and its variables.
    var_str = join([p.var_str for p in qspace.vars], ", ")
    println(io, "StateSpace: [" * var_str * "]")
    # Then print each subspace on its own line.
    for ss in qspace.subspaces
        println(io, "   - ", string(ss))
    end
end

## Test 
#xi, yi, zi = base_operators("i", qs)
#I = base_operators("I", qs)
#alpha, beta = base_operators("vars", qs)

Is = Union{Int,Vector{Int}}
function cleanup_terms(terms::Vector{Tuple{T,S}}; tol::Real=1e-12)::Vector{Tuple{T,S}} where {T<:Number,S}

    sort!(terms, by=x -> x[2])  # Sort by the index vector

    cleaned = Tuple{eltype(terms[1][1]),Vector{Is}}[]  # Result container
    i = 1
    while i <= length(terms)
        coeff_sum = terms[i][1]
        idx_vec = terms[i][2]
        j = i + 1
        while j <= length(terms) && terms[j][2] == idx_vec
            coeff_sum += terms[j][1]
            j += 1
        end
        if abs(coeff_sum) > tol
            push!(cleaned, (coeff_sum, idx_vec))
        end
        i = j
    end

    return cleaned
end


include("OperatorSets/Qubit_Pauli.jl")
include("OperatorSets/Qubit_PM.jl")
include("OperatorSets/Ladder.jl")

end # module 