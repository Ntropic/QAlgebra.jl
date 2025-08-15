export term

# ==========================================================================================================================================================
# --------> Terms from Strings <----------------------------------------------------------------------------------------------------------------------------
# ==========================================================================================================================================================
function get_op_inds(space::SubSpace, res_strings::Tuple{Vector{Vector{String}},Vector{Vector{String}},Vector{Vector{String}}})::Vector{Vector{Tuple{Number,Is}}}
    str_elements::Vector{Vector{String}} = res_strings[space.ss_outer_ind]
    all_results::Vector{Vector{Tuple{Number,Is}}} = []
    for sub_ind in space.ss_inner_ind
        curr_elements::Vector{String} = str_elements[sub_ind]
        curr_results = space.op_set.strs2ind(curr_elements)
        push!(all_results, curr_results)
    end
    return all_results
end
function string2qabstract(statespace::StateSpace, operator_str::String)::QAbstract
    # find which statespace.operator_names fits the beginning of the string 
    name = nothing
    rest::String = ""
    key_index = -1
    for (i, op_name) in enumerate(statespace.operator_names)
        if startswith(operator_str, op_name)
            name = op_name
            key_index = i 
            rest = operator_str[length(op_name)+1:end]
            break
        end
    end
    if isnothing(name)
        error("Invalid string: $operator_str, must be one of $(statespace.operator_names)")
    end
    # process rest separate by separating into pre and post ^ 
    expstr, exp = expstr_separate(rest)
    # if expstr ends with ' conjugate = true 
    conjugate = occursin("'", expstr) 
    subindex = -1
    if conjugate
        a,b = split(expstr, "'")
        if length(a) > 0
            subindex = parse(Int, a)
        end
        # find the operator in statespace 
        if length(b) > 0
            error("For an abstract operator, the Dagger symbol (') can only be followed up by an exponential, nothing else.")
        end
    else 
        if length(expstr) > 0
            subindex = parse(Int, expstr)
        end
    end
    operator_type = statespace.operatortypes[key_index]
    return QAbstract(operator_type, key_index, subindex, exp, conjugate) # subindex only printed if not -1 
end
