import ..CFunctions: reorder
export reorder

function where_defined_to_index_order(statespace::StateSpace, where_defined::Vector{Vector{Bool}})::Tuple{Vector{Int}, Vector{Int}}
    # takes where_defined and the ensemble indexes to determine the new order for both operators and variables 
    # for each element of where_defined, we shift all the true elements to the left, and all false elements to the right, we want to get the indexes of the permutation that achieves that 
    n_vars = length(statespace.vars_str)
    n_ops = length(statespace.I_op)
    op_inds = collect(1:n_ops)
    var_inds = collect(1:n_vars)
    ensemble_indexes = statespace.subspaceinfo.ensemble_indexes
    variable_indexes = statespace.where_by_ensemble_var
    for (w, c, vs)  in zip(where_defined, ensemble_indexes, variable_indexes)
        w_order = sortperm(w, rev=true)
        for v in vs
            var_inds[v] = var_inds[v][w_order]
        end
        op_inds[c] = op_inds[c][w_order]
    end
    return op_inds, var_inds
end


# QObj
function reorder(q::QTerm, index_order::Vector{Int})::QTerm
    op_indices = q.op_indices[index_order]
    return QTerm(op_indices)
end
function reorder(q::QAtomProduct, add_at_sum::Bool,  where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_index_order::Vector{Int})::QAtomProduct
    return modify_coeff_expr(q, reorder(q.coeff_fun, var_index_order), QAtom[reorder(x, index_order) for x in q.expr])
end
function reorder(q::QExpr, add_at_sum::Bool, where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_index_order::Vector{Int})::QExpr
    return QExpr(q.statespace, [reorder(qq, add_at_sum, where_defined, index_order, var_index_order) for qq in q.terms])
end
function reorder(q::T, add_at_sum::Bool, where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_index_order::Vector{Int})::T where T <: QComposite
    return modify_expr(q, reorder(q.expr, add_at_sum, where_defined, index_order, var_index_order))
end
function reorder(q::T, add_at_sum::Bool, where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_index_order::Vector{Int})::T where T <: QMultiComposite
    modify_expr(q, [reorder(qq, add_at_sum, where_defined, index_order, var_index_order) for qq in q.expr])
end
function reorder(q::QSum, add_at_sum::Bool, where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_index_order::Vector{Int})::QSum
    # define improved index_order and var_index_order
    qspace = q.statespace
    if add_at_sum
        subsystem = q.subsystem_index 
        element_indexes = q.element_indexes
        # check where subsystem is among qspace.where_ensemble
        outer_ind = findfirst(x -> x == subsystem, q.statespace.where_ensemble)
        if isnothing(outer_ind)
            error("Subsystem $subsystem not found in qspace.subspace_by_ind!")
        end
        if !all(where_defined[outer_ind][element_indexes] .== false)
            error("Summation indexes already defined, cannot sum over defined indexes!")
        end
        new_where_defined = copy.(where_defined)
        new_where_defined[outer_ind][element_indexes] .= true
        op_ind, var_inds = where_defined_to_index_order(q.statespace, new_where_defined)

        # we need to change the other parameters of sum aswell determining what is summed over 
        curr_subspace = qspace.subspaces[subsystem]
        new_element_indexes::Vector{Int} = []
        new_indexes::Vector{String} = []
        for (i, element_index) in enumerate(element_indexes)
            prev_index = curr_subspace.op_index_inds[element_index]
            new_ind = findfirst(x -> x == prev_index, op_ind)
            if isnothing(new_ind) 
                error("Not all element_indexes found in new indexing!")
            end
            # find in subspace 
            new_subind = findfirst(x -> x == new_ind, curr_subspace.op_index_inds)
            if isnothing(new_subind) 
                error("Not all element_indexes found in new indexing!")
            end
            push!(new_element_indexes, new_subind)
            push!(new_indexes, curr_subspace.keys[new_subind])
        end
        return modify_expr_indexes(q, reorder(q.expr, add_at_sum, new_where_defined, op_ind, var_inds), new_indexes, q.subsystem_index, new_element_indexes)
    else
        return modify_expr(q, reorder(q.expr, add_at_sum, where_defined, index_order, var_index_order))
    end
end

"""^
    reorder!(q::diff_QEq) -> diff_QEq

Reorders the indexes of ensemble-subspaces to the left, so that present indexes are i,j,k and not i,k,m.
This allows simplify to further simplify expressions, by removing 
"""
function reorder(q::diff_QEq)::diff_QEq
    # check index order on left side 
    where_defined_lhs = which_ensemble_acting(q.left_hand_side)
    op_inds, var_inds = where_defined_to_index_order(q.statespace, where_defined_lhs)
    # check if op_inds is not sorted (i.e. not equal to 1:length(op_inds))
    if op_inds != 1:length(op_inds) 
        # first we need to sort without changing at sums 
        left_hand_side = reorder(q.left_hand_side, false, where_defined_lhs, op_inds, var_inds)
        # then expr 
        expr = reorder(q.expr, false, where_defined_lhs, op_inds, var_inds)
        where_defined_lhs = which_ensemble_acting(left_hand_side)
        op_inds = collect(1:length(op_inds))
        var_inds = collect(1:length(var_inds))
        q = diff_QEq(q.statespace, left_hand_side, expr, q.braket)
    end
    expr = reorder(q.expr, true, where_defined_lhs, op_inds, var_inds)
    return diff_QEq(q.statespace, q.left_hand_side, expr, q.braket)
end
