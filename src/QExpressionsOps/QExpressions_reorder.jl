export reorder!

function where_defined_to_index_order(statespace::StateSpace, where_defined::Vector{Vector{Bool}})::Tuple{Vector{Int}, Vector{Int}}
    # takes where_defined and the continuum indexes to determine the new order for both operators and variables 
    # for each element of where_defined, we shift all the true elements to the left, and all false elements to the right, we want to get the indexes of the permutation that achieves that 
    n_vars = length(statespace.vars_str)
    n_ops = length(statespace.neutral_op)
    op_inds = collect(1:n_ops)
    var_inds = collect(1:n_vars)
    continuum_indexes = statespace.continuum_indexes
    variable_indexes = statespace.where_by_continuum_var
    for (w, c, vs)  in zip(where_defined, continuum_indexes, variable_indexes)
        w_order = sortperm(w, rev=true)
        for v in vs
            var_inds[v] = var_inds[v][w_order]
        end
        op_inds[c] = op_inds[c][w_order]
    end
    return op_inds, var_inds
end

### FFunction
function reorder(f::FAtom, var_index_order::Vector{Int})::FAtom
    var_exponents = f.var_exponents[var_index_order]
    return FAtom(copy(f.coeff), var_exponents)
end
function reorder(f::FSum, var_index_order::Vector{Int})::FSum
    f.terms = [reorder(ff, var_index_order) for ff in f.terms]
    return f
end
function reorder(f::FRational, var_index_order::Vector{Int})::FRational
    f.numer = reorder(f.numer, var_index_order)
    f.denom = reorder(f.denom, var_index_order)
    return f
end

# QObj
function reorder(q::QTerm, index_order::Vector{Int})::QTerm
    op_indices = q.op_indices[index_order]
    return QTerm(op_indices)
end
function reorder(q::QAtomProduct, add_at_sum::Bool,  where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_index_order::Vector{Int})::QAtomProduct
    q.expr = [reorder(qq, index_order) for qq in q.expr]
    q.coeff_fun = reorder(q.coeff_fun, var_index_order)
    return q
end
function reorder(q::qExpr, add_at_sum::Bool, where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_index_order::Vector{Int})::qExpr
    q.terms = [reorder(qq, add_at_sum, where_defined, index_order, var_index_order) for qq in q.terms]
    return q
end
function reorder(q::T, add_at_sum::Bool, where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_index_order::Vector{Int})::T where T <: QComposite
    q.expr = reorder(q.expr, add_at_sum, where_defined, index_order, var_index_order)
    return q
end
function reorder(q::T, add_at_sum::Bool, where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_index_order::Vector{Int})::T where T <: QMultiComposite
    q.expr = [reorder(qq, add_at_sum, where_defined, index_order, var_index_order) for qq in q.expr]
    return q
end
function reorder(q::QSum, add_at_sum::Bool, where_defined::Vector{Vector{Bool}}, index_order::Vector{Int}, var_index_order::Vector{Int})::QSum
    # define improved index_order and var_index_order
    qspace = q.statespace
    if add_at_sum
        subsystem = q.subsystem_index 
        element_indexes = q.element_indexes
        # check where subsystem is among qspace.where_continuum
        outer_ind = findfirst(x -> x == subsystem, q.statespace.where_continuum)
        if outer_ind == nothing
            error("Subsystem $subsystem not found in qspace.subspace_by_ind!")
        end
        if !all(where_defined[outer_ind][element_indexes] .== false)
            error("Summation indexes already defined, cannot sum over defined indexes!")
        end
        new_where_defined = deepcopy(where_defined)
        new_where_defined[outer_ind][element_indexes] .= true
        op_ind, var_inds = where_defined_to_index_order(q.statespace, new_where_defined)

        # we need to change the other parameters of sum aswell determining what is summed over 
        #qspace.continuum_indexes, qspace.neutral_continuum_op 
        curr_subspace = qspace.subspaces[subsystem]
        new_element_indexes::Vector{Int} = []
        new_indexes::Vector{String} = []
        for (i, element_index) in enumerate(element_indexes)
            prev_index = curr_subspace.op_index_inds[element_index]
            new_ind = findfirst(x -> x == prev_index, op_ind)
            if new_ind == nothing
                error("Not all element_indexes found in new indexing!")
            end
            # find in subspace 
            new_subind = findfirst(x -> x == new_ind, curr_subspace.op_index_inds)
            if new_subind == nothing
                error("Not all element_indexes found in new indexing!")
            end
            push!(new_element_indexes, new_subind)
            push!(new_indexes, curr_subspace.keys[new_subind])
        end
        q.expr = reorder(q.expr, add_at_sum, new_where_defined, op_ind, var_inds)
        q.indexes = new_indexes
        q.element_indexes = new_element_indexes
    else
        q.expr = reorder(q.expr, add_at_sum, where_defined, index_order, var_index_order)
    end
    return q
end

"""
    reorder!(q::diff_QEq) -> diff_QEq

Reorders the indexes of continuum-subspaces to the left, so that present indexes are i,j,k and not i,k,m.
This allows simplify to further simplify expressions, by removing 
"""
function reorder!(q::diff_QEq)::diff_QEq
    # check index order on left side 
    q = copy(q)#deepcopy(q)
    where_defined_lhs = which_continuum_acting(q.left_hand_side)
    op_inds, var_inds = where_defined_to_index_order(q.statespace, where_defined_lhs)
    # check if op_inds is not sorted (i.e. not equal to 1:length(op_inds))
    if op_inds != 1:length(op_inds) 
        # first we need to sort without changing at sums 
        left_hand_side = reorder(q.left_hand_side, false, where_defined_lhs, op_inds, var_inds)
        # then expr 
        expr = reorder(q.expr, false, where_defined_lhs, op_inds, var_inds)
        where_defined_lhs = which_continuum_acting(left_hand_side)
        op_inds = collect(1:length(op_inds))
        var_inds = collect(1:length(var_inds))
        q = diff_QEq(q.statespace, left_hand_side, expr, q.braket)
    end
    expr = reorder(q.expr, true, where_defined_lhs, op_inds, var_inds)
    return diff_QEq(q.statespace, q.left_hand_side, expr, q.braket)
end
