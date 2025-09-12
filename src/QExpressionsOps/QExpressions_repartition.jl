import ..CFunctions: repartition
export repartition

# Fix this function vor variables 
function where_defined_to_index_order(statespace::StateSpace, where_defined::Vector{Vector{Bool}})::Tuple{Vector{Int}, Vector{Tuple{Int, Int}}}
    # takes where_defined and the ensemble indexes to determine the new order for both operators and variables 
    # for each element of where_defined, we shift all the true elements to the left, and all false elements to the right, we want to get the indexes of the permutation that achieves that 
    function permutation_moves(p::Vector{Int})::Vector{Tuple{Int, Int}}  # helper to extract the moves 
        swaps::Vector{Tuple{Int, Int}} = []
        for (i, pi) in enumerate(p)
            if pi > i   # only count once
                push!(swaps, (i, pi))
            end
        end
        return swaps
    end

    n_ops = length(statespace.I_op)
    n_vars = length(statespace.params)

    param_info = statespace.param_info 
    subspace_info = statespace.subspace_info
    
    op_inds = collect(1:n_ops)
    var_inds = collect(1:n_vars)

    ensemble_indexes = subspace_info.ensemble_indexes
    where_ensembles = subspace_info.where_ensembles
    for (outer, w, c)  in zip(where_ensembles, where_defined, ensemble_indexes) # iterate over ensemble subspaces
        w_order = sortperm(w, rev=true)
        op_inds[c] = op_inds[c][w_order]
        inner_perms = permutation_moves(w_order)
        for inner_perm in inner_perms  
            curr_perm_params = map_by_subspace(SubSpaceIndex(outer, inner_perm[1], c[inner_perm[1]]), SubSpaceIndex(outer, inner_perm[2], c[inner_perm[2]]), param_info)
            var_inds = var_inds[curr_perm_params]
        end
    end
    var_tuples::Vector{Tuple{Int, Int}} = []
    @inbounds for i in 1:n_vars 
        if var_inds[i] != i 
            push!(var_tuples, (i, var_inds[i]))
        end
    end
    return op_inds, var_tuples
end


# QObj
"""
    repartition!(q::diff_QEq) -> diff_QEq

Reorders the indexes of ensemble-subspaces to the left, so that present indexes are for example  i,j,k and not i,_,k,m .
This allows simplify to further simplify expressions, by removing differences only in indexing parameters
"""
function repartition(q::QTerm, index_order::Vector{Int})::QTerm
    op_indices = q.op_indices[index_order]
    return QTerm(op_indices)
end
function repartition(q::QAtomProduct, add_at_sum::Bool,  where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_tuples::Vector{Tuple{Int, Int}})::QAtomProduct
    return modify_coeff_expr(q, repartition(q.coeff_fun, var_tuples), QAtom[repartition(x, index_order) for x in q.expr])
end
function repartition(q::QExpr, add_at_sum::Bool, where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_tuples::Vector{Tuple{Int, Int}})::QExpr
    return QExpr(q.statespace, [repartition(qq, add_at_sum, where_defined, index_order, var_tuples) for qq in q.terms])
end
function repartition(q::T, add_at_sum::Bool, where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_tuples::Vector{Tuple{Int, Int}})::T where T <: QComposite
    return modify_coeff_expr(q, repartition(q.coeff_fun, var_tuples), repartition(q.expr, add_at_sum, where_defined, index_order, var_tuples))
end
function repartition(q::T, add_at_sum::Bool, where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_tuples::Vector{Tuple{Int, Int}})::T where T <: QMultiComposite
    return modify_coeff_expr(q, repartition(q.coeff_fun, var_tuples), [repartition(qq, add_at_sum, where_defined, index_order, var_tuples) for qq in q.expr])
end
function repartition(q::QSum, add_at_sum::Bool, where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_tuples::Vector{Tuple{Int, Int}})::QSum
    # define improved index_order and var_index_order
    qspace = q.statespace
    info =  qspace.subspace_info
    if add_at_sum
        indexes = q.indexes

        new_where_defined = copy.(where_defined)
        for index in indexes 
            ensemble = Index2Ensemble(index, info)
            if new_where_defined[ensemble][index.inner] == true
                index_str = Index2String(index, info)
                error("Summation index $index_str already defined, cannot sum over defined indexes!")
            else
                new_where_defined[ensemble][index.inner] = true 
            end
        end
        op_ind, var_tuples_new = where_defined_to_index_order(qspace, new_where_defined) # permutation vectors 

        # we need to change the other parameters of sum aswell determining what is summed over 
        new_indexes::Vector{SubSpaceIndex} = []
        for index in indexes
            # find index in new order 
            new_ind_expanded = findfirst(==(index.expanded), op_ind)  # theres probably a better way to do this than findfirst 
            diff = new_ind_expanded - index.expanded 
            new_inner = index.inner + diff 
            if new_inner < 1 
                error("New index no longer in the same ensemble subspace.")
            end 
            # find in subspace 
            push!(new_indexes, SubSpaceIndex(index.outer, new_inner, new_ind_expanded))
        end
        return modify_expr_indexes(q, repartition(q.expr, add_at_sum, new_where_defined, op_ind, var_tuples_new), new_indexes)
    else
        return modify_expr(q, repartition(q.expr, add_at_sum, where_defined, index_order, var_tuples))
    end
end

function repartition(q::diff_QEq)::diff_QEq
    # check index order on left side 
    where_defined_lhs = which_ensemble_acting(q.left_hand_side)
    op_inds, var_inds = where_defined_to_index_order(q.statespace, where_defined_lhs)
    # check if op_inds is not sorted (i.e. not equal to 1:length(op_inds))
    if op_inds != 1:length(op_inds) 
        # first we need to sort without changing at sums 
        left_hand_side = repartition(q.left_hand_side, false, where_defined_lhs, op_inds, var_inds)
        # then expr 
        expr = repartition(q.expr, false, where_defined_lhs, op_inds, var_inds)
        where_defined_lhs = which_ensemble_acting(left_hand_side)
        op_inds = collect(1:length(op_inds))
        var_inds = collect(1:length(var_inds))
        return diff_QEq(q.statespace, left_hand_side, expr, q.do_braket)
    end
    expr = repartition(q.expr, true, where_defined_lhs, op_inds, var_inds)
    return diff_QEq(q.statespace, q.left_hand_side, expr, Val(:nosimp), do_braket=q.do_braket)
end
