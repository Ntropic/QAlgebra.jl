export decollision_QSum


# Repartitioning for collision free QSum indexing: in a  QSum:
# 0th initialize empty where_defined
# 1st create subsitutions in case of collision of current indexes with where_defined
# 2nd update where_defined  -> copy where defined, op_tuples, and var_tuples everywhere
# 3rd apply repartition to terms.

struct QSumDecollisionInds 
    init::Bool
    where_acting::Vector{BitVector}  
    op_tuples::Vector{Tuple{Int, Int}}
    var_tuples::Vector{Tuple{Int, Int}}
    var_inds::Vector{Int}
    function QSumDecollisionInds(init::Bool, where_acting::Vector{BitVector}, op_tuples::Vector{Tuple{Int, Int}}, var_tuples::Vector{Tuple{Int, Int}}, var_inds::Vector{Int})
        new(init, where_acting, op_tuples, var_tuples, var_inds) 
    end
    function QSumDecollisionInds(q::QSum)
        qspace = q.qspace
        subspace_info = qspace.subspace_info
        where_acting::Vector{BitVector} = [falses(n) for n in subspace_info.how_many_sum_by_ensemble]
        for subspace_ind in iter_all_indexes(q)
            ensemble_ind, summation_ind = Index2Ensemble_and_Summation(subspace_ind, subspace_info)
            @assert ensemble_ind != 0 "Found an invalid Summation index: $(Index2String(subspace_ind)). "
            where_acting[ensemble_ind][summation_ind] = true 
        end
        op_tuples = Vector{Tuple{Int, Int}}()
        var_tuples = Vector{Tuple{Int, Int}}()
        var_inds = Vector{Int}() #collect(1:length(qspace.params))
        new(false, where_acting, op_tuples, var_tuples, var_inds) 
    end
end 

# Detect summation index collisions, and find the first ensemble index that is free, to switch to, returns collision::Bool, the (potentially) updated where_acting, and a potentially new summation index, 
function collision_find_first_free(where_acting::Vector{BitVector}, subspace_ind::SubSpaceIndex, subspace_info::SubSpaceInfo)::Tuple{Bool, Vector{BitVector}, SubSpaceIndex}
    ensemble, summation = Index2Ensemble_and_Summation(subspace_ind, subspace_info)
    if where_acting[ensemble][summation] # found a collision 
        # find first free 
        free_sum_ind = findfirstfreeafterbefore(where_acting[ensemble], summation)
        @assert !isnothing(free_sum_ind) "No free summation index left in ensemble. "
        where_acting[ensemble][free_sum_ind] = true # QSums need to copy this! to not everright each other 
        return true, where_acting, SummationIndex2SubSpaceIndex(subspace_ind.outer, ensemble, free_sum_ind, subspace_info)
    else 
        where_acting[ensemble][summation] = true 
        return false, where_acting, subspace_ind
    end
end

@inline function _apply_permutation!(vec::Vector{T}, perm::Vector{Int}) where {T}
    length(vec) <= 1 && return
    if all(@inbounds perm[i] == i for i in eachindex(perm))
        return
    end
    tmp = vec[perm]
    copyto!(vec, tmp)
end

function _sort_block!(block::ConstrainedIndexBlock)
    idxs = block.indexes
    length(idxs) <= 1 && return block
    perm = sortperm(idxs; by=expanded)
    _apply_permutation!(idxs, perm)
    _apply_permutation!(block.constraints, perm)
    _assert_no_duplicate_indexes(idxs)
    return block
end

function _relocate_index!(block::ConstrainedIndexBlock, pos::Int, old_index::SubSpaceIndex, new_index::SubSpaceIndex)
    block.indexes[pos] = new_index
    old_inner = old_index.inner
    new_inner = new_index.inner
    old_inner == new_inner && return block
    self_row = block.constraints[pos]
    old_self = self_row[old_inner]
    _swap_constraint_columns!(block.constraints, old_inner, new_inner)
    self_row[new_inner] = old_self
    return block
end

@inline function _assert_no_duplicate_indexes(indexes::Vector{SubSpaceIndex})
    length(indexes) <= 1 && return
    @inbounds for i in 2:length(indexes)
        prev = indexes[i - 1]
        curr = indexes[i]
        prev.expanded == curr.expanded && error("QSum decollision failed: duplicate summation index detected after relabelling.")
    end
end

function update_QSumDecollisionInds(q::QSum, d::QSumDecollisionInds)::Tuple{QSumDecollisionInds, Vector{ConstrainedIndexBlock}}
    qspace = q.qspace
    subspace_info = qspace.subspace_info
    param_info = qspace.param_info

    new_op_tuples = Tuple{Int,Int}[]
    inds_tuples = Tuple{SubSpaceIndex,SubSpaceIndex}[]
    new_where = copy.(d.where_acting)

    blocks = Vector{ConstrainedIndexBlock}(undef, length(q.blocks))

    @inbounds for (block_idx, original) in enumerate(q.blocks)
        working_block = nothing
        for (pos, index) in enumerate(original.indexes)
            collision, new_where, new_index = collision_find_first_free(new_where, index, subspace_info)
            if collision
                push!(new_op_tuples, (index.expanded, new_index.expanded))
                push!(inds_tuples, (index, new_index))
                if working_block === nothing
                    working_block = _clone_block(original)
                end
                _relocate_index!(working_block, pos, index, new_index)
            end
        end
        if working_block === nothing
            blocks[block_idx] = original
        else
            _sort_block!(working_block)
            blocks[block_idx] = working_block
        end
    end

    if isempty(new_op_tuples)
        return QSumDecollisionInds(d.init, new_where, d.op_tuples, d.var_tuples, d.var_inds), blocks
    end

    var_inds = collect(1:length(qspace.params))
    @inbounds for (index, new_index) in Base.Iterators.reverse(inds_tuples)
        curr_perm_params = map_by_subspace(index, new_index, param_info)
        var_inds = var_inds[curr_perm_params]
    end
    if d.init
        var_inds = var_inds[d.var_inds]
    end

    var_tuples = [(i, var_inds[i]) for i in eachindex(var_inds) if var_inds[i] != i]

    return QSumDecollisionInds(true, new_where, vcat(new_op_tuples, d.op_tuples), var_tuples, var_inds), blocks
end

function decollision_QSum_product(q1::QSum, q2::QSum)::Vector{QComposite}
    # assumes that each QSum is already internally decollisioned! 
    qspace = q1.qspace
    subspace_info = qspace.subspace_info
    where_acting::Vector{BitVector} = which_summations_acting(q1, subspace_info) 
    decollision = QSumDecollisionInds(false, where_acting, Tuple{Int, Int}[], Tuple{Int, Int}[], Int[])
    inner_terms = decollision_QSum(q2, decollision)

    base_terms  = QComposite[]
    nested_sums = QSum[]
    for term in inner_terms
        if term isa QSum
            push!(nested_sums, term)
        else
            push!(base_terms, term)
        end
    end

    out_terms = QComposite[]
    if !isempty(base_terms)
        new_inner_expr = q1.expr * QExpr(qspace, base_terms)
        push!(out_terms, QSum(qspace, new_inner_expr, q1.blocks))
    end
    for nested in nested_sums
        merged_expr = q1.expr * nested.expr
        merged_blocks = merge_blocks(q1.blocks, nested.blocks)
        push!(out_terms, QSum(qspace, merged_expr, merged_blocks))
    end
    return out_terms
end

"""
    decollision_QSum(q::QSum) -> Vector{QComposite}

Resolve collisions between summation indexes by relabelling clashing
indices and repartitioning coefficients. The returned vector contains the
collision-free terms that replace the original `QSum`.
"""
function decollision_QSum(q::QSum)::Vector{QComposite} 
    decollision = QSumDecollisionInds(q)
    decollision, blocks = update_QSumDecollisionInds(q, decollision)
    new_q = QSum(q.qspace, q.expr, blocks)
    return decollision_QSum(new_q, decollision, Val(:noupdate))
end
function decollision_QSum(q::QSum, decollision::QSumDecollisionInds, ::Val{:noupdate})::Vector{QComposite} #assume it is already updated 
    qspace = q.qspace
    base_terms  = QComposite[]
    nested_sums = QSum[]
    inner = decollision_QSum(q.expr, decollision)
    for term in inner.terms
        if term isa QSum
            push!(nested_sums, term)
        else
            push!(base_terms, term)
        end
    end
    out_terms = QComposite[]
    if !isempty(base_terms)
        push!(out_terms, QSum(qspace, QExpr(qspace, base_terms), q.blocks)) # no more clone blocks
    end
    for nested in nested_sums
        merged_blocks = merge_blocks(q.blocks, nested.blocks)
        push!(out_terms, QSum(qspace, nested.expr, merged_blocks))
    end
    return out_terms
end

function decollision_QSum(q::QSum, decollision::QSumDecollisionInds)::Vector{QComposite}
    decollision, blocks = update_QSumDecollisionInds(q, decollision)
    new_q = QSum(q.qspace, q.expr, blocks)
    return decollision_QSum(new_q, decollision, Val(:noupdate))
end


function decollision_QSum(q::QTerm, decollision::QSumDecollisionInds, qspace::QSpace)::QTerm
    op_indices = copy(q.op_indices)
    @inbounds for (old_ind, new_ind) in decollision.op_tuples
        @assert isnumeric(q, new_ind, qspace) "Cannot decollision QTerm, because new summation index is already in use, albeit undefined!"
        op_indices[old_ind], op_indices[new_ind] = op_indices[new_ind], op_indices[old_ind]
    end
    return QTerm(op_indices)
end
function decollision_QSum(q::QAbstract, decollision::QSumDecollisionInds, qspace::QSpace)::QAbstract
    return add_to_index_map(q, decollision.op_tuples)
end
function decollision_QSum(q::QAtomProduct, decollision::QSumDecollisionInds)::Vector{QComposite}
    if decollision.init
        return QComposite[modify_coeff_expr(q, repartition(q.coeff_fun, decollision.var_tuples), QAtom[decollision_QSum(x, decollision, q.qspace) for x in q.expr])]
    else
        return QComposite[q] 
    end
end
function decollision_QSum(q::QExpr, decollision::QSumDecollisionInds)::QExpr
    new_terms = QComposite[]
    sizehint!(new_terms, length(q.terms))
    for t in q.terms
        append!(new_terms, decollision_QSum(t, decollision))
    end
    return QExpr(q.qspace, new_terms)
end
function decollision_QSum(q::T, decollision::QSumDecollisionInds)::Vector{QComposite} where T <: QComposite
    return QComposite[modify_coeff_expr(q, repartition(q.coeff_fun, decollision.var_tuples), decollision_QSum(q.expr, decollision))]
end
function decollision_QSum(q::T, decollision::QSumDecollisionInds)::Vector{QComposite} where T <: QMultiComposite
    return QComposite[modify_coeff_expr(q, repartition(q.coeff_fun, decollision.var_tuples), [decollision_QSum(qq, decollision) for qq in q.expr])]
end
