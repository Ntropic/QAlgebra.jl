
"""
    base_operators(statespace::StateSpace, letter::String; do_fun::Bool=false) -> Union{Int,Vector{Int}}
    base_operators(ss:StateSpace; do_fun::Bool=false) -> Tuple{Dict{String,qExpr},Dict{String,qExpr},Dict{String,Function}}

Returns variables and/or operators in the state space `ss`.
Specifc variables/operators can be selected by passing a string `letter`.
If no `letter` is passed, the function returns a tuple of 3 dictionaries:
- The first dictionary contains the variables in the state space, with their corresponding qExpr objects.
- The second dictionary contains the operators in the state space, with their corresponding qExpr objects.
- The third dictionary contains the abstract operators in the state space either as a callable function to specify the subtype or as a standard qExpr.
If you pass "vars", it will return a tuple with elements for each variable
    - do_fun specifies if abstract operators are returned as functions, that can be called with no arguments or with an integer to specify the subindex of the abstract operator. 
"""
function base_operators(statespace::StateSpace, letter::String; do_fun::Bool=false) 
    my_ops::Vector{qExpr} = []
    var_exponents = zeros(Int, length(statespace.vars))
    neutral_operator = [s.op_set.neutral_element for s in statespace.subspaces for key in s.keys]
    index = 1
    if letter == "I"
        return qExpr(statespace, qTerm[])
    end
    if letter == "vars"
        for i in 1:length(statespace.vars)
            var_exponents[i] += 1
            push!(my_ops, qExpr(qAtomProduct(statespace, 1, var_exponents, qTerm[]), statespace))
            var_exponents[i] -= 1
        end
        if length(my_ops) > 1
            return tuple(my_ops...)
        elseif length(my_ops) == 1
            return my_ops[1]
        else
            return nothing
        end
    end
    for (i, var) in enumerate(statespace.vars)
        if letter == var.var_name || letter == var.var_str
            var_exponents[i] += 1
            return qExpr(statespace, qAtomProduct(statespace, 1, var_exponents, qTerm[]))
        end
    end
    # Option 2 for vars: gather all the ones for which it occurs in the name, colelct those return if vector is not empty
    curr_ops = qExpr[]
    for (i, var) in enumerate(statespace.vars)
        if occursin(letter, var.var_name) || occursin(letter, var.var_str)
            var_exponents[i] += 1
            push!(curr_ops, qExpr(statespace, qAtomProduct(statespace, 1, var_exponents, qTerm[])))
            var_exponents[i] -= 1
        end
    end
    if length(curr_ops) > 0 
        return curr_ops
    end

    # check subspaces
    for sub in statespace.subspaces
        for key in sub.keys
            if key == letter
                keys = sub.keys
                op_set = sub.op_set
                base_ops = op_set.base_ops
                for base_op in base_ops
                    curr_operator = copy(neutral_operator)
                    curr_operator[index] = base_op
                    push!(my_ops, qExpr(statespace, qTerm(copy(curr_operator))))
                end

                # non base ops # in a Dict 
                non_base_ops = op_set.non_base_ops
                for (key, ops) in non_base_ops
                    curr_terms::Vector{qAtomProduct} = qAtomProduct[]
                    for op in ops
                        curr_operator = copy(neutral_operator)
                        curr_operator[index] = op[2]
                        coeff = op[1]
                        curr_prod = qAtomProduct(statespace,coeff, var_exponents, qTerm(copy(curr_operator)))
                        push!(curr_terms, curr_prod)
                    end
                    push!(my_ops, qExpr(statespace, curr_terms))
                end

                if length(my_ops) > 1
                    return tuple(my_ops...)
                else
                    return my_ops[1]
                end
            end
            index += 1
        end
    end
    # check for abstract operators
    for (key_index, (name, operatortype)) in enumerate(zip(statespace.operator_names, statespace.operatortypes))
        if name == letter
            if do_fun 
                return (subindex=-1) -> qExpr(statespace, qAbstract(operatortype, statespace.key_index, subindex))
            else
                return qExpr(statespace, qAbstract(operatortype, key_index))
            end
        elseif occursin(name, letter)
            return qExpr(statespace, string2qabstract(statespace, replace(letter, "_" => "")))
        end
    end
    error("No variable, subspace component or abstract operator with key starting with '$letter' found in the state space.")
end
function base_operators(statespace::StateSpace; do_fun::Bool=false)::Tuple{Dict{String,qExpr},Dict{String,qExpr}, Dict{String, Union{Function, qExpr}}}
    # return 2 dicctionaries, one with the vars and one with the operators 
    var_dict::Dict{String,qExpr} = Dict()
    op_dict::Dict{String,qExpr} = Dict()
    CRone = one(ComplexRational)
    var_exponents = zeros(Int, length(statespace.vars))
    #abstract_dict::Dict{String,qExpr} = Dict()
    neutral_operator = statespace.neutral_op
    for (i, vars_str) in enumerate(statespace.vars_str)
        var_exponents[i] += 1
        var_dict[vars_str] = qExpr(statespace, qAtomProduct(statespace, CRone, var_exponents, qTerm[]))
        var_exponents[i] -= 1
    end
    index = 1
    for sub in statespace.subspaces
        op_set = sub.op_set
        base_ops = op_set.base_ops
        for key in sub.keys
            for base_op in base_ops
                curr_operator = copy(neutral_operator)
                curr_operator[index] = base_op
                term = qTerm(curr_operator)
                curr_name = op_set.op2str(base_op, key)
                op_dict[curr_name] = qExpr(statespace, term)
            end

            # non base ops # in a Dict 
            non_base_ops = op_set.non_base_ops
            for (key, ops) in non_base_ops
                curr_terms::Vector{qAtomProduct} = qAtomProduct[]
                for op in ops
                    curr_operator = copy(neutral_operator)
                    curr_operator[index] = op[2]
                    coeff = op[1]
                    curr_prod = qAtomProduct(statespace, coeff, var_exponents, qTerm(copy(curr_operator)))
                    push!(curr_terms, curr_prod)
                end
                op_dict[key] = qExpr(statespace, curr_terms)
            end
            index += 1
        end
    end
    op_dict["I"] = qExpr(statespace, qTerm[])
    # now abstract operators 
    abstract_dict::Dict{String, Union{Function, qExpr}} = Dict()
    for (key_index, (name, operatortype)) in enumerate(zip(statespace.operator_names, statespace.operatortypes))
        if do_fun
            abstract_dict[name] = (subindex=-1) -> qExpr(statespace, qAbstract(operatortype, key_index, subindex))
        else
            abstract_dict[name] = qExpr(statespace, qAbstract(operatortype, key_index))
        end
    end
    return var_dict, op_dict, abstract_dict
end

