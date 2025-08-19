struct OpDefinitions
    operators::Vector{Tuple{String, Vector{String}}}
    commute_fun::Function 
    check_n::Int
    function OpDefinitions(operators...; commute_fun::Union{Nothing,Function}=nothing, check_n=5) 
        new_operators::Vector{Tuple{String, Vector{String}}} = []
        for op in operators
            if isa(op, String)
                push!(new_operators, brace_separate(s))
            elseif isa(op, Symbol)
                push!(new_operators, (string(op), String[]))
            elseif isa(op, Tuple)
                push!(new_operators, op) 
            else 
                error("Invalid operator type: $op, must be String, Tuple{String, Vector{String}} or Symbol!")
            end
        end
        return new(new_operators, commute_fun, check_n)  
    end
end

struct OperatorTypeInfo
    operator_names::Vector{String} # names of the operators (e
    operator_types::Vector{OperatorType}
    commutation_matrix::Matrix{Bool} # if true, then the operators commute. If false, they probably don't 
    commute_fun::Function 
    function OperatorTypeInfo(operator_types::Vector{OperatorType}; commute_fun::Union{Nothing,Function}=nothing, check_n=5) 
        operator_names::Vector{String} = [ot.name for ot in operator_types]
        commutation_matrix, fun = generate_commutation_matrix(operator_types, commute_fun, check_n)
        return new(operator_names, operator_types, commutation_matrix, fun)
    end
end

""" 
    OperatorType(name::String, hermitian::Bool=false, unitary::Bool=false, acting_ss::Union{Nothing, Vector{Bool}}=nothing)


Define an operator Type, by declaring its name and properties such as hermitian and unitary. 
"""
struct OperatorType
    name::String 
    name_sym::Symbol
    operator_index::Int
    hermitian::Bool 
    unitary::Bool 
    of_time::Bool
    acting_ss::Vector{Bool}   # formerly acting_ss
    expanded_ss_acting::Vector{Bool} # indices of non-trivial components in the op_indices picture for bath subsystems
end

function OpDefinitions2OperatorType(opdefs::OpDefinitions, subspace_definitions::SubSpaceDefinitions)::Vector{OperatorType}  
     return [SingleOpDefinition2OperatorType(op_def, i, cond, subspace_definitions) for (i, (op_def, cond) in enumerate(opdefs.operators))]
end
function SingleOpDefinition2OperatorType(op_str::String, operator_index::Int, conditions::Vector{String}, subspace_definitions::SubSpaceDefinitions)::OperatorType
    used_symbols::Vector{Symbol} = subspace_definitions.used_symbols
    op_sym::Symbol = Symbol(op_str)
    conditions_sym::Vector{Symbol} = Symbol.(conditions)
    if (op_sym in used_symbols)
        error("Operator $op_str already defined for either an Subspace definition or another Abstract Operator!") 
    end

    key_symbols::Vector{Symbol} = Symbol[]  # Add all
    for subspace in subspace_definitions.subspaces
        push!(key_symbols, subspace.key_symbol)
    end
    if length(conditions_sym) != length(unique(conditions_sym))
        throw(ArgumentError("Duplicate conditions in operator definition $(conditions_sym)")
    end

    of_time, unitary, hermitian = false, false, false 
    conditions_str::Vector{String} = []
    contains_negation::Vector{Bool} = []
    reduced_conditions_sym::Vector{Symbol} = []
    for condition_sym in conditions_sym
        if condition_sym == :t
            of_time = true
        elseif condition_sym == :U || condition_sym == :unitary
            unitary = true
        elseif condition_sym == :H || condition_sym == :hermitian
            hermitian = true
        else 
            condition_str = string(condition_sym)
            if startswith(condition_str, "!") 
                push!(contains_negation, true)
                push!(conditions_str, condition_str[2:end])
                push!(reduced_conditions_sym, Symbol(condition_str[2:end]))
            else
                push!(contains_negation, false)
                push!(conditions_str, condition_str)
                push!(reduced_conditions_sym, condition_sym)
            end
        end
    end
    for reduced_conditions_sym in conditions_sym
        if !(condition_sym in key_symbols)
            throw(ArgumentError("Invalid condition: $condition_sym, must be one of: $(accepted_conditions)")) 
        end 
    end
    # check if all negations are the same 
    if !(all(contains_negation)  || all(.!contains_negation))
        error("Mixed negations in conditions are not allowed, got for $op_str the negation pattern $contains_negation.")
    end
    negation::Bool = contains_negation[1]
    acting_ss::Vector{Bool} = [negation for _ in subspace_keys]

    for (i, key_symbol) in enumerate(key_symbols)  
        if key_symbol in redduced_conditions_sym 
            acting_ss[i] = !acting_ss[i]
        end
    end
    expanded_ss_acting::Vector{Bool} = fill(false, dim)
    for (s_bool, subspace) in zip(acting_ss, subspaces)
        if s_bool == true
            for ind in subspace.op_index_inds
                expanded_ss_acting[ind] = true
            end
        end
    end
    return OperatorType(op_str, op_sym, operator_index, hermitian, unitary, of_time, acting_ss, expanded_ss_acting)
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

struct OperatorTypeInfo
    operator_names::Vector{String} # names of the operators (e
    operator_types::Vector{OperatorType}
    commutation_matrix::Matrix{Bool} # if true, then the operators commute. If false, they probably don't 
    commute_fun::Function 
    function OperatorTypeInfo(operator_types::Vector{OperatorType}; commute_fun::Union{Nothing,Function}=nothing, check_n=5) 
        operator_names::Vector{String} = [ot.name for ot in operator_types]
        commutation_matrix, fun = generate_commutation_matrix(operator_types, commute_fun, check_n)
        return new(operator_names, operator_types, commutation_matrix, fun)
    end
end