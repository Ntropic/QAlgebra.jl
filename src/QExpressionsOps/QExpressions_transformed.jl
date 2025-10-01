export QAtomOrdered, QNeutral, OrderbyOperator

# ===============> Sorting op_indices by subspaces, and returning the ensemble permutations 
"""
    decompose_sorted_blocks(op_indices::Vector{Is}, qspace::QSpace)
        -> (Vector{Vector{Is}}, Vector{Vector{Vector{Int}}})

Decompose a full operator index vector into subspace blocks.

- `blocks`:
  - One `Vector{Is}` per subspace.
  - For ensemble subspaces: contains only non-neutral operators, sorted by operator type.
  - For non-ensemble subspaces: contains the single operator as a 1-element vector.

- `ensemble_indexes`:
  - One entry per ensemble subspace (in order).
  - Each entry is a vector of vectors of positions, grouped by operator identity.
  - Example: for `[X, I, Y, X]` with neutral `I`, result is `[[1,4],[3]]`.
"""
function decompose_sorted_blocks(op_indices::Vector{Is}, qspace::QSpace)::Tuple{Vector{Vector{Is}}, Vector{Vector{Vector{Int}}}}
    
    nsub::Int = length(qspace.subspaces)
    blocks::Vector{Vector{Is}} = Vector{Vector{Is}}(undef, nsub)
    ensemble_indexes::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}()

    index::Int = 1
    for (sidx, subspace) in enumerate(qspace.subspaces)
        if subspace.is_ensemble_ss
            ensemble_size::Int = subspace.ensemble_size
            neutral_element::Is = subspace.op_set.neutral_element

            curr_nontrivial::Vector{Is} = Is[]
            curr_positions::Vector{Int} = Int[]
            for j in 1:ensemble_size
                op::Is = op_indices[index + j - 1]
                if op != neutral_element
                    push!(curr_nontrivial, op)
                    push!(curr_positions, j)
                end
            end

            pairs::Vector{Tuple{Is,Int}} = collect(zip(curr_nontrivial, curr_positions))
            sort!(pairs, by = first)

            blocks[sidx] = [p[1] for p in pairs]

            grouped_pos::Vector{Vector{Int}} = Vector{Vector{Int}}()
            last_op::Union{Nothing,Is} = nothing
            for (op, pos) in pairs
                if last_op === nothing || op != last_op
                    push!(grouped_pos, [pos])
                    last_op = op
                else
                    push!(grouped_pos[end], pos)
                end
            end
            push!(ensemble_indexes, grouped_pos)

            index += ensemble_size
        else
            blocks[sidx] = [op_indices[index]]
            index += 1
        end
    end

    return blocks, ensemble_indexes
end
"""
    recompose_op_indices(blocks::Vector{Vector{Is}}, 
                         ensemble_indexes::Vector{Vector{Vector{Int}}}, 
                         qspace::QSpace) -> Vector{Is}

Rebuild the original `op_indices` vector from its block decomposition.

- Fills ensemble subspaces with their neutral element first,
  then restores non-trivial operators at the recorded positions.
- Non-ensemble subspaces are copied directly.
"""
function recompose_op_indices(blocks::Vector{Vector{Is}}, ensemble_indexes::Vector{Vector{Vector{Int}}}, qspace::QSpace)::Vector{Is}
    op_indices::Vector{Is} = Is[]
    ens_counter::Int = 1

    for (sidx, subspace) in enumerate(qspace.subspaces)
        if subspace.is_ensemble_ss
            ensemble_size::Int = subspace.ensemble_size
            neutral::Is = subspace.op_set.neutral_element

            curr_ops::Vector{Is} = fill(neutral, ensemble_size)

            flat_ops::Vector{Is} = blocks[sidx]
            grouped_pos::Vector{Vector{Int}} = ensemble_indexes[ens_counter]

            pos_counter::Int = 1
            for group in grouped_pos
                for pos in group
                    curr_ops[pos] = flat_ops[pos_counter]
                    pos_counter += 1
                end
            end

            append!(op_indices, curr_ops)
            ens_counter += 1
        else
            append!(op_indices, blocks[sidx])
        end
    end

    return op_indices
end



""" 
QAtomOrdered is a QTerm, but where the operators within ensembles are sorted by operator index.  
"""
struct QAtomOrdered <: QComposite
    qspace::QSpace
    coeff_fun::CFunction
    op_indices::Vector{Vector{Is}}
    ensemble_indexes::Vector{Vector{Int}}
    time_index::Int
end
struct QNeutral <: QComposite 
    qspace::QSpace
    coeff_fun::CFunction 
    time_index::Int 
end

"""
    OrderbyOperator(q::QObj; lt=isless)

Replace simple QAtomProducts (consisting only of QTerms) into QAtomOrdered, to sort 
"""
OrderbyOperator(q::QAtomOrdered) = q
OrderbyOperator(q::QAbstract) = error("Cannot order QAbstract. Must be substituted before.")
OrderbyOperator(q::QTerm) = 

function OrderbyOperator(q::T; lt=isless)::T where T <: QComposite 
    ordered_expr = OrderbyOperator.(q.expr, lt=isless)
    return modify_expr(q, ordered_expr)
end
function OrderbyOperator(q::T; lt=isless)::T where T <: QMultiComposite
    ordered_factors = [OrderbyOperator(term; lt=lt) for term in q.expr]
    return modify_expr(q, ordered_factors)
end

function OrderbyOperator(q::QExpr; lt=isless)
    ordered_terms = Vector{QComposite}(undef, length(q.terms))
    @inbounds for i in eachindex(q.terms)
        ordered_terms[i] = OrderbyOperator(q.terms[i]; lt=lt)
    end
    return QExpr(q.qspace, ordered_terms, Val(:nosimp))
end

# Core OrderbyOperator here!
function OrderbyOperator(q::QAtomProduct; lt=isless)
    n = length(q.expr)
    if n == 0
        return QNeutral(q.qspace, q.coeff_fun, q.time_index)
    elseif n == 1
        blocks, ensemble_indexes = decompose_sorted_blocks(q.expr[1].op_indices, )
        return QAtomOrdered(q.qspace, q.coeff_fun, blocks, ensemble_indexes, q.expr[1].time_index)
    end
    # shouldn'T be needed, but whatever
    ordered_atoms = Vector{QAtomOrdered}(undef, n)
    @inbounds for i in 1:n
        blocks, ensemble_indexes = decompose_sorted_blocks(q.expr[1].op_indices, )
        ordered_atoms[i] =  QAtomOrdered(q.qspace, q.coeff_fun, blocks, ensemble_indexes, q.expr[1].time_index)
    end
    return QCompositeProduct(q.qspace, q.coeff_fun, ordered_atoms)
end