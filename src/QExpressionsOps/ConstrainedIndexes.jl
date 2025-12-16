export are_all_neq
import ..QAlgebra: sorted_push_unique!

"""
    NeqConstraint(lhs, rhs)

Inequality constraint between two summation indices. Accepts any `SumIndexInput`
at construction; indices are normalised to `SubSpaceIndex` objects when a
`QSpace` is available.
"""

struct NeqConstraint{T<:SumIndexInput}
    lhs::T
    rhs::T
end
function neq(lhs::T, rhs::T) where {T<:SumIndexInput}  # convenience constructor
    return NeqConstraint(lhs, rhs)
end
function neqconstraint_of_SubSpaceIndex(qspace::QSpace, constraint::NeqConstraint)::NeqConstraint{SubSpaceIndex}
    lhs = _to_subspace_index(qspace, constraint.lhs)
    rhs = _to_subspace_index(qspace, constraint.rhs)
    @assert lhs.outer == rhs.outer "Neq-constraints must be of the same ensemble subspace, got outer subspace indices $(lhs.outer) and $(rhs.outer) for indices constraints $(constraint.lhs) and $(constraint.rhs)." 
    return NeqConstraint{SubSpaceIndex}(lhs, rhs)
end

function sort_indices_and_constraints_into_ensemble_blocks(qspace::QSpace, indices::Vector{SubSpaceIndex}, neq::Bool, constraints::Vector{NeqConstraint{SubSpaceIndex}})::Vector{ConstrainedIndexBlock}
    subspace_info = qspace.subspace_info 
    where_ensembles = subspace_info.where_ensembles
    how_many_by_ensemble = subspace_info.how_many_by_ensemble
    non_sum = subspace_info.how_many_non_sum_by_ensemble
    # initialize empty blocks for the ensembles
    blocks::Vector{ConstrainedIndexBlock} = [ConstrainedIndexBlock(outer, ensemble_size, how_many_non_sum) for (outer, ensemble_size, how_many_non_sum) in zip(where_ensembles, how_many_by_ensemble, non_sum)]
    for idx in indices 
        ensemble_ind, summation_ind = Index2Ensemble_and_Summation(idx, subspace_info)
        # checks 
        if summation_ind < 1
            subspace = qspace.subspaces[idx.outer]
            available = subspace.keys[non_sum[ensemble_ind]+1:end]
            if !isempty(available)
                error("Index $(Index2String(idx, subspace_info)) is not a summation index. Available summation keys: $available.")
            else
                error("No summation indices defined for ensemble $(subspace.key).");
            end
        end
        if neq
            push_index_falses!(blocks[ensemble_ind], idx, summation_ind)
        else
            push_index_trues!(blocks[ensemble_ind], idx, summation_ind)
        end
    end
    if neq && !isempty(constraints)
        error("You shouldn't define non-equality via neq=true and define inequality constraints.")
    else
        for constraint in constraints
            lhs, rhs = constraint.lhs, constraint.rhs
            ensemble_lhs = subspace_info.ensemble_index_by_subspace_index[lhs.outer] 
            @assert ensemble_lhs != 0 "Neq-constraints must be of an ensemble subspace. got subspace index $(lhs.outer) for $(Index2String(lhs, subspace_info))." 
            how_many_non_sum = non_sum[ensemble_lhs]
            if constraint.lhs.inner <= how_many_non_sum && constraint.lhs.inner <= how_many_non_sum + ensemble_lhs      
                subspace = qspace.subspaces[lhs.outer]
                available = subspace.keys[non_sum[ensemble_lhs]+1:end]
                if !isempty(available)
                    error("Index $(Index2String(lhs, subspace_info)) and $(Index2String(rhs, subspace_info)) are not summation indices. Available summation keys: $available.")
                else
                    error("No summation indices defined for ensemble $(subspace.key).");
                end
            else # Apply to block 
                apply_neq!(blocks[ensemble_lhs], constraint)
            end
        end
    end
    return blocks 
end


# =================================================> Index & Constraint Blocks <=============================================================


struct ConstrainedIndexBlock
    outer::Int
    ensemble_size::Int
    how_many_non_sum::Int
    indices::Vector{SubSpaceIndex}
    constraints::Vector{BitVector}    # true -> can be equal , false -> is neq 
    function ConstrainedIndexBlock(outer::Int, ensemble_size::Int, how_many_non_sum::Int, indices::Vector{SubSpaceIndex}=SubSpaceIndex[], constraints::Vector{BitVector}=BitVector[])
        # assume lengths are equal! 
        return new(outer, ensemble_size, how_many_non_sum, indices, constraints)
    end
end
Base.length(block::ConstrainedIndexBlock) = length(block.indices)
Base.getindex(block::ConstrainedIndexBlock, i::Int) = block.indices[i]
copy_empty(block::ConstrainedIndexBlock)::ConstrainedIndexBlock = ConstrainedIndexBlock(block.outer, block.ensemble_size, block.how_many_non_sum)
Base.copy(block::ConstrainedIndexBlock)::ConstrainedIndexBlock = ConstrainedIndexBlock(
    block.outer,
    block.ensemble_size, block.how_many_non_sum, copy(block.indices), map(copy, block.constraints))

@inline _lt_subspaceindex(a::SubSpaceIndex, b::SubSpaceIndex) =
    a.expanded == b.expanded ? (a.outer == b.outer ? a.inner < b.inner : a.outer < b.outer) : a.expanded < b.expanded

function _insert_index_unique!(indices::Vector{SubSpaceIndex}, idx::SubSpaceIndex)::Int
    old_len = length(indices)
    sorted_push_unique!(indices, idx; lt=_lt_subspaceindex)
    if length(indices) == old_len
        throw(ArgumentError("Duplicate index $(idx) found in constrained index block."))
    end
    return searchsortedfirst(indices, idx; lt=_lt_subspaceindex)
end

function push_index_trues!(block::ConstrainedIndexBlock, idx::SubSpaceIndex, summation::Int)::ConstrainedIndexBlock # assume idx belongs into this block! 
    i = _insert_index_unique!(block.indices, idx)
    insert!(block.constraints, i, trues(block.ensemble_size))
    return block
end
function push_index_falses!(block::ConstrainedIndexBlock, idx::SubSpaceIndex, summation::Int)::ConstrainedIndexBlock # assume idx belongs into this block! 
    i = _insert_index_unique!(block.indices, idx)
    insert!(block.constraints, i, falses(block.ensemble_size))
    return block
end
function apply_neq!(block::ConstrainedIndexBlock, constraint::NeqConstraint{SubSpaceIndex})::ConstrainedIndexBlock
    block_index_lhs = findfirst(x -> x == (constraint.lhs), block.indices)
    block_index_rhs = findfirst(x -> x == (constraint.rhs), block.indices)
    # one of them needs to be an index, for each that is an index we remove the others outer from the constraints (i.e. make it false )
    done = false 
    if !isnothing(block_index_lhs)
        block.constraints[block_index_lhs][constraint.rhs.inner] = false
        done = true
    end
    if !isnothing(block_index_rhs)
        block.constraints[block_index_rhs][constraint.lhs.inner] = false 
        done = true 
    end
    @assert done "None of the indices described in Neq-Condition was present in block."
    return block 
end

function merge_blocks(blocks_a::Vector{ConstrainedIndexBlock}, blocks_b::Vector{ConstrainedIndexBlock})
    @assert length(blocks_a) == length(blocks_b) "Blocks need to be acting on the same set of Ensemble subspaces, and have to have the same number of blocks."
    merged = Vector{ConstrainedIndexBlock}(undef, length(blocks_a))
    @inbounds for (i, (a, b)) in enumerate(zip(blocks_a, blocks_b))
        if isempty(a.indices)
            merged[i] = b
            continue
        elseif isempty(b.indices)
            merged[i] = a
            continue
        end

        @assert a.outer == b.outer "Ensemble mismatch-. outer differs."
        @assert a.ensemble_size == b.ensemble_size  "Ensemble size mismatch."
        @assert a.how_many_non_sum == b.how_many_non_sum  "Non-sum count mismatch."

        # start from a’s sorted unique base
        combined_indices     = copy(a.indices)
        combined_constraints = [copy(c) for c in a.constraints]

        # check previous conditions for neq condition transitivity
        for (bi, bc) in zip(b.indices, b.constraints)
            for (ai, ac) in zip(a.indices, a.constraints)
                if xor(bc[ai.inner], ac[bi.inner])
                    error("Inconsistent inequality conditions. Cannot merge blocks. ")
                end
            end
        end   

        # insert b’s indices one by one
        for (idx, constr) in zip(b.indices, b.constraints)
            pos = _insert_index_unique!(combined_indices, idx)
            insert!(combined_constraints, pos, copy(constr))
        end
        merged[i] = ConstrainedIndexBlock(a.outer, a.ensemble_size, a.how_many_non_sum, combined_indices, combined_constraints,)
    end

    return merged
end

# Helpers for print functions !!!
# ================================================>  Some final helpers <==================================================================================

# flatten all block indices into a single vector in block order
function _flatten_indices(blocks::Vector{ConstrainedIndexBlock})::Vector{SubSpaceIndex}
    # concatenate all indices across blocks preserving order
    total = sum(length(block) for block in blocks)
    result = Vector{SubSpaceIndex}(undef, total)
    cursor = 1
    for block in blocks
        for idx in block.indices
            result[cursor] = idx
            cursor += 1
        end
    end
    return result
end
function _build_blocks(qspace::QSpace, indices::Vector{SubSpaceIndex}, constraints::Vector{BitVector})::Vector{ConstrainedIndexBlock}
    info = qspace.subspace_info
    where_ensembles = info.where_ensembles
    ensemble_sizes = info.how_many_by_ensemble
    non_sum = info.how_many_non_sum_by_ensemble
    n_ensembles = length(where_ensembles)
    temp_indices = [SubSpaceIndex[] for _ in 1:n_ensembles]
    temp_constraints = [BitVector[] for _ in 1:n_ensembles]

    @assert length(indices) == length(constraints) "Index and constraint vectors must match in length."

    @inbounds for (idx, row) in zip(indices, constraints)
        ensemble = info.ensemble_index_by_subspace_index[idx.outer]
        ensemble != 0 || error("Index $(Index2String(idx, info)) does not belong to an ensemble subspace.")
        push!(temp_indices[ensemble], idx)
        push!(temp_constraints[ensemble], BitVector(row))
    end

    blocks = Vector{ConstrainedIndexBlock}(undef, n_ensembles)
    @inbounds for ensemble in 1:n_ensembles
        idxs = temp_indices[ensemble]
        rows = temp_constraints[ensemble]
        if !isempty(idxs)
            perm = sortperm(idxs; by=expanded)
            idxs = idxs[perm]
            rows = rows[perm]
        end
        blocks[ensemble] = ConstrainedIndexBlock(where_ensembles[ensemble], ensemble_sizes[ensemble], non_sum[ensemble], idxs, rows)
    end
    return blocks
end

# Find which indices have neq conditions, (and find the conditions) and which can be equal
function eq_counter_and_neq_indices_by_block(block::ConstrainedIndexBlock, where_acting_block::BitVector, subspace_info::SubSpaceInfo)::Tuple{Int, Vector{NeqConstraint{SubSpaceIndex}}}
    non_sum = block.how_many_non_sum
    neq_constraints::Vector{NeqConstraint{SubSpaceIndex}} = Vector{NeqConstraint{SubSpaceIndex}}()
    eq_counter = 0
    for (ind, constraint) in zip(block.indices, block.constraints)
        curr_inner = ind.inner
        for i in vcat(1:non_sum, curr_inner+1:length(where_acting_block))
            if where_acting_block[i] 
                if !constraint[i]
                    push!(neq_constraints, neq(ind, SubSpaceIndex(ind.outer, i, subspace_info)))
                else
                    eq_counter += 1
                end
            end
        end
    end
    return (eq_counter, neq_constraints)
end
function eq_counter_and_neq_indices(blocks::Vector{ConstrainedIndexBlock}, where_acting::Vector{BitVector}, subspace_info::SubSpaceInfo)::Tuple{Int, Vector{NeqConstraint{SubSpaceIndex}}}
    eq_counter = 0
    neq_constraints::Vector{NeqConstraint{SubSpaceIndex}} = Vector{NeqConstraint{SubSpaceIndex}}()
    for (block, where_acting_block) in zip(blocks, where_acting)
        new_count, new_inds = eq_counter_and_neq_indices_by_block(block, where_acting_block, subspace_info)
        eq_counter += new_count
        append!(neq_constraints, new_inds)
    end 
    return (eq_counter, neq_constraints)
end

""" 
    are_all_neq(block::ConstrainedIndexBlock, where_acting_block::BitVector, subspace_info::SubSpaceInfo) -> Bool 
    are_all_neq(blocks::Vector{ConstrainedIndexBlock}, where_acting_blocks::Vector{BitVector}, subspace_info::SubSpaceInfo -> Bool
    are_all_neq(q::AbstractQSum{A}, where_acting_blocks::Vector{BitVector}) -> Bool

Are all conditions neq in ConstrainedIndexBlock or Vector of ConstrainedIndexBlocks or AbstractQSum (such as QSum and QInt).
"""
function are_all_neq(block::ConstrainedIndexBlock, where_acting_block::BitVector, subspace_info::SubSpaceInfo)::Bool
    non_sum = block.how_many_non_sum
    for (ind, constraint) in zip(block.indices, block.constraints)
        curr_inner = ind.inner
        for i in vcat(1:non_sum, curr_inner+1:length(where_acting_block))
            if where_acting_block[i] 
                if constraint[i]
                    return false
                end
            end
        end
    end
    return true
end
function are_all_neq(blocks::Vector{ConstrainedIndexBlock}, where_acting_blocks::Vector{BitVector}, subspace_info::SubSpaceInfo)::Bool
    return all(are_all_neq(block, where_acting_block, subspace_info) for (block, where_acting_block) in zip(blocks, where_acting_blocks))
end


const NeqAction = Tuple{SubSpaceIndex, Int}
const NeqBranch = Tuple{ConstrainedIndexBlock, BitVector, Vector{NeqAction}}
# Push this into QSum_modify
"""
    neq_expand(block, where_defined) -> Vector{NeqBranch}

Enumerate all ways to resolve allowed equalities for `block` against the provided
`where_defined` mask. Each returned tuple contains:
  * a cloned block with updated constraint rows,
  * a copy of `where_defined` describing the remaining available ensemble slots,
  * the list of equality actions `(index, column)` applied while descending this branch.
"""
function neq_expand(block::ConstrainedIndexBlock, where_defined::BitVector)::Vector{NeqBranch}
    length(where_defined) == block.ensemble_size || error("where_defined length must equal block ensemble size.")
    results = NeqBranch[]
    _neq_expand!(results, copy(block), copy(where_defined), 1, NeqAction[])
    return results
end

function _neq_expand!(results::Vector{NeqBranch}, block::ConstrainedIndexBlock, where_defined::BitVector, row_idx::Int, actions::Vector{NeqAction})::Nothing
    if row_idx > length(block.indices)
        push!(results, (block, copy(where_defined), copy(actions)))
        return nothing
    end

    row = block.constraints[row_idx]
    idx = block.indices[row_idx]
    candidates = Int[]
    max_col = idx.inner - 1
    max_col < 1 || @inbounds for col in 1:max_col
        if row[col] && where_defined[col]
            push!(candidates, col)
        end
    end

    @inbounds for col in candidates
        eq_block = copy(block)
        eq_where = copy(where_defined)
        eq_idx = eq_block.indices[row_idx]
        deleteat!(eq_block.indices, row_idx)
        deleteat!(eq_block.constraints, row_idx)
        for row in eq_block.constraints
            row[eq_idx.inner] = true
        end
        eq_where[col] = false
        eq_where[idx.inner] = false
        new_actions = copy(actions)
        push!(new_actions, (eq_idx, col))
        _neq_expand!(results, eq_block, eq_where, row_idx, new_actions)
    end

    next_block = copy(block)
    if !isempty(candidates)
        row_next = next_block.constraints[row_idx]
        for col in candidates
            row_next[col] = false
        end
    end
    _neq_expand!(results, next_block, where_defined, row_idx + 1, actions)
    return nothing
end

# Helper to copy blocks <=============================== Remove dependency of these functions
@inline function _swap_constraint_columns!(rows::Vector{BitVector}, a::Int, b::Int)
    a == b && return
    @inbounds for row in rows
        row[a], row[b] = row[b], row[a]
    end
end
