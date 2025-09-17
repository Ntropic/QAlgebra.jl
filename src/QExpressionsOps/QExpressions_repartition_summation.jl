# Repartitioning for collision free QSum indexing: in a  QSum:
# 0th initialize empty where_defined
# 1st create subsitutions in case of collision of current indexes with where_defined
# 2nd update where_defined 
# 3rd apply repartition to terms.


# Detect summation index collisions, and find the first ensemble index that is free, to switch to, returns collision::Bool, the (potentially) updated where_acting, and a potentially new summation index, 
function collision_find_first_free(where_acting::Vector{BitVector}, subspace_ind::SubSpaceIndex, info::SubSpaceInfo)::Tuple{Bool, Vector{BitVector}, SubSpaceIndex}
    ensemble, summation = Index2Ensemble_and_Summation(subspace_ind, info)
    if where_acting[ensemble][summation] # found a collision 
        # find first free 
        free_sum_ind = findfirst(!, where_acting[ensemble])
        where_acting[ensemble][free_sum_ind] = true # QSums need to copy this! to not everright each other 
        return true, where_acting, SummationIndex2SubSpaceIndex(subspace_ind.outer, ensemble, free_sum_ind, info)
    else
        return false, where_acting, subspace_ind
    end
end

# Continue here ==> Rewrite the permutations to instead detect the collisions and change them. 
function summation_collisions_to_index_order(statespace::StateSpace, where_defined::Vector{BitVector}, curr_indexes::Vector{SubSpaceIndex})::Tuple{Vector{Int}, Vector{Tuple{Int, Int}}}
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
    i = 0
    for (outer, w, c)  in zip(where_ensembles, where_defined, ensemble_indexes) # iterate over ensemble subspaces
        i += 1
        w_order = sortperm(w, rev=true)
        op_inds[c] = op_inds[c][w_order]
        inner_perms = permutation_moves(w_order)
        for inner_perm in inner_perms  
            curr_perm_params = map_by_subspace(SummationIndex2SubSpaceIndex(outer, ensemble, inner_perm[1], subspace_info), SummationIndex2SubSpaceIndex(outer, ensemble, inner_perm[2], subspace_info), param_info)
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

function decollision_QSum(q::QSum)::QSum 
    where_defined::Vector{BitVector} = [zeros(Bool, s) for s in subspace_info.how_many_sum_by_ensemble]
    return decollision_QSum(q, where_defined)
end

function decollision_QSum(q::QSum, where_defined::Vector{BitVector})::QSum
    qspace = q.statespace
    info   = qspace.subspace_info
    # clone where_defined to mutate
    new_where_defined = copy.(where_defined)

    # mark all indexes in all blocks as defined
    for index in iter_all_indexes(q)
        ensemble = Index2Ensemble(index, info)
        if new_where_defined[ensemble][index.inner]
            index_str = Index2String(index, info)
            error("Summation index $index_str already defined, cannot sum over defined indexes!")
        end
        new_where_defined[ensemble][index.inner] = true
    end
    op_ind, var_tuples_new = where_defined_to_index_order(qspace, new_where_defined)

    # remap indexes in each block to the new order
    new_indexes= SubSpaceIndex[]
    for index in q.eq_indexes
        new_expanded = findfirst(==(index.expanded), op_ind)
        diff = new_expanded - index.expanded
        new_inner = index.inner + diff
        if new_inner < 1
            error("New index no longer in the same ensemble subspace.")
        end
        push!(new_indexes, SubSpaceIndex(index.outer, new_inner, new_expanded))
    end
    new_blocks = Vector{Vector{SubSpaceIndex}}()
    for blk in q.blocks
        new_blk = SubSpaceIndex[]
        for index in blk
            new_expanded = findfirst(==(index.expanded), op_ind)
            diff = new_expanded - index.expanded
            new_inner = index.inner + diff
            if new_inner < 1
                error("New index no longer in the same ensemble subspace.")
            end
            push!(new_blk, SubSpaceIndex(index.outer, new_inner, new_expanded))
        end
        push!(new_blocks, new_blk)
    end
    return modify_expr_indexing( q, repartition(q.expr, false, new_where_defined, op_ind, var_tuples_new), new_indexes, new_blocks )
end