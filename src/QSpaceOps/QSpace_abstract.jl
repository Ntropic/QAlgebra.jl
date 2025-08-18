""" 
    OperatorType(name::String, hermitian::Bool=false, unitary::Bool=false, acting_ss::Union{Nothing, Vector{Bool}}=nothing)


Define an operator Type, by declaring its name and properties such as hermitian and unitary. 
"""
struct OperatorType
    name::String 
    operator_index::Int
    hermitian::Bool 
    unitary::Bool 
    of_time::Bool
    acting_ss::Vector{Bool}   # formerly acting_ss
    expanded_ss_acting::Vector{Bool} # indices of non-trivial components in the op_indices picture for bath subsystems
    function OperatorType(name::String; operator_index::Int, hermitian::Bool=false, unitary::Bool=false, of_time::Bool, acting_ss::Vector{Bool}, expanded_ss_acting::Vector{Bool})
        if length(name) == 0
            error("OperatorType name cannot be empty.")
        end
        new(name, operator_index, hermitian, unitary, of_time, subsystems,acting_ss, expanded_ss_acting)
    end
end
function string2operator_type(operator_index::Int, s::String, subspace_keys::Vector{String}, dim::Int, subspaces::Vector{SubSpace})
    # use name(U,H) notation, separate the name from the brace. and parse the , separated elements within the brace to determine which properties are true (by default all are false) 
    hermitian = false 
    unitary = false 
    of_time = false
    acting_ss::Vector{Bool} = [true for _ in subspace_keys]
    name::String = ""
    modified_subspace_keys::Vector{String} = ["!"*key for key in subspace_keys]
    if occursin("(", s)  # Check if there is a brace in the string
        name, brace = split(s, "(")  # Output: ["As", "U,H"]
        if length(brace) > 0
            brace = split(brace, ")")[1]  # Output: "U,H"
            tokens = split(brace, ",")  # Output: ["U", "H"]
            tokens = [strip(t) for t in tokens]
            of_time = "t" in tokens
            hermitian = "H" in tokens
            unitary = "U" in tokens 
            token_list = vcat(["t", "H", "U"], subspace_keys, modified_subspace_keys)
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
                    acting_ss[i] = s in tokens
                end
            else
                acting_ss = fill(true, length(subspace_keys))
                for (i, s) in enumerate(modified_subspace_keys)
                    if s in tokens
                        acting_ss[i] = false
                    end
                end
            end
        else
            throw(ArgumentError("Invalid Operator string: $s"))
        end
    else
        name = s
    end
    expanded_ss_acting::Vector{Bool} = fill(false, dim)
    for (s_bool, subspace) in zip(acting_ss, subspaces)
        if s_bool == true
            for ind in subspace.op_index_inds
                expanded_ss_acting[ind] = true
            end
        end
    end
    return OperatorType(name, operator_index, hermitian=hermitian, unitary=unitary, of_time=of_time, acting_ss=acting_ss, expanded_ss_acting=expanded_ss_acting)
end
function operator_type2string(p::OperatorType, time_index::Int=0)
    curr_str = p.name 
    if any([p.hermitian, p.unitary])
        curr_str *= "("
        if p.of_time 
            curr_str *= "t"
            if time_index > 0
                curr_str *= str2sub(String(time_index))
            end
            curr_str *= ","
        end
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

function operatertypes2commutator_matrix(optypes::Vector{OperatorType})::Matrix{Bool}
    # matrix of operatortypes commutation, true ==> the two commutators commute
    mat = Matrix{Bool}(undef, length(optypes), length(optypes))
    for i in 1:length(optypes)
        mat[i,i] = true
        for j in i+1:length(optypes)
            # check for any collisions of the operator_type acting_ss 
            mat[i,j] = !any(optypes[i].acting_ss .& optypes[j].acting_ss)
            mat[j,i] = mat[i,j]
        end
    end
    return mat
end

function accepts_n_args(f, n::Int)
    any(m -> length(m.sig.parameters) - 1 == n, methods(f))
end
# Generates function that takes Operator indexes and returns commutations
function gen_commutes_function(optypes::Vector{OperatorType}, fun::Union{Function, Nothing} = nothing; check_n::Int=1)::Tuple{Matrix{Bool}, Function}
    commutation_matrix = operatertypes2commutator_matrix(optypes) # if sth commutes for commutation matrix it must also commute for fun
    # but if sth doesn't commute in commutation matrix it might still commute for fun!
    if isnothing(fun)
        # assume this function structure 
        fun = (a_key, a_sub, a_dag, b_key, b_sub, b_dag) -> commutation_matrix[a_key, b_key]
        return commutation_matrix, fun
    else
        if accepts_n_args(fun, 6)
            # function that takes 6 arguments
            return commutation_matrix, fun
        end 
        for optype in optypes
            if !optype.hermitian && !optype.unitary
                error("Commutation functions with 2 or 4 arguments require that the operators are either hermitian or unitary, to generalize to Daggered variants of commutation. ")
            end
        end
        n = length(optypes) 
        if accepts_n_args(fun, 2)
            # check every combination of a_key and b_key against commutation matrix 
            for i in 1:n, j in 1:i-1
                if commutation_matrix[i,j]
                    if !fun(i,j) # must output bool
                        op_A_str = operator_type2string(optypes[i])
                        op_B_str = operator_type2string(optypes[j])  
                        error("For [$op_A_str, $op_B_str] the provided commutation function (fun) must return true for arguments ($i, $j), because the operators are acting on non overlapping subspaces!")
                    end
                end
            end
            return commutation_matrix, fun
        elseif accepts_n_args(fun, 4)
            for i in 1:n, j in 1:i-1
                if commutation_matrix[i,j]
                    # create check_n random integers picked in the interval 0 to 4*check_n
                    sub_inds_i = rand(0:4*check_n, check_n)
                    sub_inds_j = rand(0:4*check_n, check_n)
                    for (sub_i, sub_j) in zip(sub_inds_i, sub_inds_j)
                    if !fun(i, sub_i, j, sub_j) # must output bool
                        op_A_str = operator_type2string(optypes[i])
                        op_B_str = operator_type2string(optypes[j])  
                        error("For [$op_A_str, $op_B_str] the provided commutation function (fun) must return true for arguments ($i, sub index 1, $j, sub index 2) with arbitrary sub indexes, because the operators are acting on non overlapping subspaces! Returned false for sub indexes ($sub_i, $sub_j). ")
                    end
                end
            end
            return commutation_matrix, fun 
        else 
            error("The commutation function (fun) must accept 2, 4 or 6 arguments!")
        end
    end
end

