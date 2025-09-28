#### Flatten 
export neq, flatsums, complexsums

""" 
    flatsums(q::QObj)::Bool 

Checks if QSums are flattened in QObj. Returns true if there are no nested QSums. 
"""
flatsums(q::QSum, has_sum::Bool=false)::Bool = has_sum ? false : flatsums(q.expr, has_sum) 
(flatsums(q::T, has_sum::Bool=false)::Bool) where {T<:QComposite} = flatsums(q.expr, has_sum) 
flatsums(q::QAtomProduct, has_sum::Bool=false) = true 
(flatsums(q::T, has_sum::Bool=false)::Bool) where {T<:QMultiComposite} = all(flatsums.(q.expr, has_sum))  
flatsums(q::QExpr, has_sum::Bool=false)::Bool = all(flatsums.(q.terms, has_sum))  
flatsums(q::diffQEq)::Bool = flatsums(q.expr)


""" 
    complexsums(q::QObj)::Bool 

Returns ture if any QSum is nested or in a QComposite. 
"""
complexsums(q::QSum, in_complex::Bool=false)::Bool = in_complex ? true : complexsums(q.expr, true)  
(complexsums(q::T, in_complex::Bool=false)::Bool) where {T<:QComposite} = complexsums(q.expr, true) 
complexsums(q::QAtomProduct, in_complex::Bool=false) = false 
(complexsums(q::T, in_complex::Bool=false)::Bool) where {T<:QMultiComposite} = any(complexsums.(q.expr, true))  
complexsums(q::QExpr, in_complex::Bool=false)::Bool = any(complexsums.(q.terms, in_complex))  
complexsums(q::diffQEq)::Bool = complexsums(q.expr)



#### first output is (changed), then vectors of terms and then of coefficients 
# change from index1 to index2
function term_equal_indexes(expr, args...) # Base method to error
    throw(MethodError(term_equal_indexes, (typeof(expr), args...)))
end
# multiplies from the left 
function term_equal_indexes(term::QTerm, index1::SubSpaceIndex, index2::SubSpaceIndex, subspace::SubSpace)::Tuple{Bool, Vector{QTerm}, Vector{ComplexRational}}
    ind1, ind2 = index1.expanded, index2.expanded
    op1 = term.op_indices[ind1]
    op2 = term.op_indices[ind2]
    neutral = subspace.op_set.neutral_element
    if op1 == neutral && op2 == neutral
        return false, QTerm[term], ComplexRational[ComplexRational(1,0,1)]
    end
    results = subspace.op_set.op_product(op1, op2)
    new_terms = QTerm[]
    new_coeffs = ComplexRational[]
    for (coeff, op) in results
        op_indices = copy(term.op_indices)
        op_indices[ind2] = op
        op_indices[ind1] = neutral
        push!(new_terms, QTerm(op_indices, term.time_index))
        push!(new_coeffs, coeff)
    end
    return true, new_terms, new_coeffs  
end 

function term_equal_indexes(abstract::QAbstract, index1::SubSpaceIndex, index2::SubSpaceIndex, subspace::SubSpace)::Tuple{Bool, Vector{QAbstract}, Vector{ComplexRational}}
    expanded_ss_acting = abstract.operator_type.expanded_ss_acting
    if expanded_ss_acting[index2.expanded]
        return true, QAbstract[add_to_index_map(abstract, (index1, index2))], ComplexRational[ComplexRational(1,0,1)]
    end
    # append this rule to the index map 
    return false, QAbstract[abstract], ComplexRational[ComplexRational(1,0,1)]
end

function term_equal_indexes(q::QAtomProduct, index1::SubSpaceIndex, index2::SubSpaceIndex, subspace::SubSpace, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, Vector{QAtomProduct}}
    changed_any = false
    term_variants = Vector{Vector{QAtom}}()
    coeff_variants = Vector{Vector{ComplexRational}}()
    for atom in q.expr
        changed, variants, coeffs = term_equal_indexes(atom, index1, index2, subspace)
        push!(term_variants, variants)
        push!(coeff_variants, coeffs)
        changed_any |= changed  # Check if any term was changed
    end
    changed, new_coeff_fun = CFunctions.term_equal_indexes(q.coeff_fun, coeff_ind_order)
    changed_any |= changed
    if !changed_any
        return false, [q]
    end
    # Generate all combinations (cartesian product) of updated terms
    combinations = Iterators.product(term_variants...)
    coeff_combinations = Iterators.product(coeff_variants...)
    simplified_products = QAtomProduct[]
    for (combo, coeff_combo) in zip(combinations, coeff_combinations)
        new_expr = collect(combo)
        factor = prod(coeff_combo)
        new_prod = QAtomProduct(q.qspace, new_coeff_fun*factor, new_expr, q.separate_expectation_values, q.braket)
        push!(simplified_products, new_prod)
    end
    return true, simplified_products
end

function term_equal_indexes(qexpr::QExpr, index1::SubSpaceIndex, index2::SubSpaceIndex, subspace::SubSpace, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, Vector{QExpr}}
    changed_any = false
    elements = QComposite[]
    for t in qexpr.terms
        changed, variants = term_equal_indexes(t, index1, index2, subspace, coeff_ind_order)
        append!(elements, variants)
        changed_any |= changed  # Check if any term was changed
    end
    if !changed_any
        return false, [qexpr]
    end
    return true, [QExpr(qexpr.qspace, elements)]
end
#T <: QComposite case
function term_equal_indexes(q::T, index1::SubSpaceIndex, index2::SubSpaceIndex, subspace::SubSpace, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, Vector{T}} where T<:QComposite
    changed, variants = term_equal_indexes(q.expr, index1, index2, subspace, coeff_ind_order)
    if !changed
        return false, [q]
    end
    results = Vector{T}()
    for v in variants
        for new_term in modify_expr(q, v)
            new_term isa T || error("modify_expr returned $(typeof(new_term)), expected $(T).")
            push!(results, new_term)
        end
    end
    return true, results
end
#T <: QMultiComposite case
function term_equal_indexes(q::T, index1::SubSpaceIndex, index2::SubSpaceIndex, subspace::SubSpace, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, Vector{T}} where T<:QMultiComposite
    changed, variants = term_equal_indexes(q.expr, index1, index2, subspace, coeff_ind_order)
    if !changed
        return false, [q]
    end
    results = Vector{T}()
    for v in variants
        for new_term in modify_expr(q, v)
            new_term isa T || error("modify_expr returned $(typeof(new_term)), expected $(T).")
            push!(results, new_term)
        end
    end
    return true, results
end

"""
    changed_indices(mapvec::Vector{Int}) -> Vector{Tuple{Int,Int}}

Given a mapping vector `mapvec` where `mapvec[i]` is the new index for old index `i`,
return a list of `(i, mapvec[i])` pairs for all `i` where the mapping changes
(i.e. `mapvec[i] != i`).
"""
function changed_indices(mapvec::Vector{Int})::Vector{Tuple{Int, Int}}
    changes = Tuple{Int,Int}[]
    for i in eachindex(mapvec)
        if mapvec[i] != i
            push!(changes, (i, mapvec[i]))
        end
    end
    return changes
end

"""
    neq(qeq::QExpr) -> QExpr

Transform sums into neq sums, where all indexes are different from each other, and returns a flattened QExpr with neq sums. 
Considers all cases of the sums, simplifying the cases in which indexes are the same, which then reduces the order of the sum (i.e. a sum_{j} x_i y_j => sum_{j} x_i y_j + im*z_i, where we used x_i*y_i=im*z_i).
"""
function neq(q::QObj, do_abstract::Bool=false)::QObj
    return q
end
function neq(q::QAtomProduct, do_abstract::Bool=false)::QAtomProduct
    return q 
end
function neq(q::T, do_abstract::Bool=false)::T where {T<:QComposite}
    return only(modify_expr(q, neq(q.expr, do_abstract)))
end
function neq(q::T, do_abstract::Bool=false)::T where {T<:QMultiComposite}
    return only(modify_expr(q, neq.(q.expr, do_abstract)))
end

struct NeqState
    expr::QExpr
    blocks::Vector{ConstrainedIndexBlock}
    where_defined::Vector{BitVector}
end

@inline function _is_all_distinct(q::QSum)::Bool
    for block in q.blocks
        for (idx, row) in zip(block.indexes, block.constraints)
            @inbounds begin
                row[idx.inner] || return false
                for (col, flag) in enumerate(row)
                    col == idx.inner && continue
                    flag && return false
                end
            end
        end
    end
    return true
end

function _enforce_all_distinct(blocks::Vector{ConstrainedIndexBlock})::Vector{ConstrainedIndexBlock}
    new_blocks = clone_blocks(blocks)
    for block in new_blocks
        for (i, idx) in enumerate(block.indexes)
            row = falses(block.ensemble_size)
            row[idx.inner] = true
            block.constraints[i] = row
        end
    end
    return new_blocks
end

function _apply_equalities(expr::QExpr, actions::Vector{NeqAction}, qspace::QSpace)::Vector{QExpr}
    variants = QExpr[expr]
    info = qspace.subspace_info
    for (from_idx, column) in actions
        new_variants = QExpr[]
        target_idx = SubSpaceIndex(from_idx.outer, column, info)
        subspace = qspace.subspaces[from_idx.outer]
        coeff_inds = changed_indices(map_by_subspace(from_idx, target_idx, qspace.param_info))
        for variant in variants
            _, exprs = term_equal_indexes(variant, from_idx, target_idx, subspace, coeff_inds)
            append!(new_variants, exprs)
        end
        variants = new_variants
    end
    return variants
end

function _process_block(state::NeqState, ensemble::Int, qspace::QSpace)::Vector{NeqState}
    block = state.blocks[ensemble]
    branches = neq_expand(block, state.where_defined[ensemble])
    new_states = NeqState[]
    for (branch_block, branch_where, actions) in branches
        expr_variants = _apply_equalities(state.expr, actions, qspace)
        for expr_variant in expr_variants
            new_blocks = copy(state.blocks)
            new_blocks[ensemble] = branch_block
            new_where = copy(state.where_defined)
            new_where[ensemble] = branch_where
            push!(new_states, NeqState(expr_variant, new_blocks, new_where))
        end
    end
    return new_states
end

function _expand_qsum_states(q::QSum, where_defined::Vector{BitVector})::Vector{NeqState}
    initial_blocks = clone_blocks(q.blocks)
    initial_where = copy.(where_defined)
    states = NeqState[NeqState(q.expr, initial_blocks, initial_where)]
    for ensemble in eachindex(q.blocks)
        next_states = NeqState[]
        for state in states
            append!(next_states, _process_block(state, ensemble, q.qspace))
        end
        states = next_states
    end
    return states
end

function _collect_results(q::QSum, states::Vector{NeqState})::Vector{QComposite}
    results = QComposite[]
    for state in states
        if all(isempty(block.indexes) for block in state.blocks)
            append!(results, state.expr.terms)
        else
            distinct_blocks = _enforce_all_distinct(state.blocks)
            push!(results, QSum(q.qspace, state.expr, distinct_blocks))
        end
    end
    return results
end

function neq_qsum(s::QSum, do_abstract::Bool=false)
    _is_all_distinct(s) && return QExpr(s.expr.qspace, [s])
    where_defined = which_ensemble_acting(s, do_abstract=do_abstract)
    states = _expand_qsum_states(s, where_defined)
    terms = _collect_results(s, states)
    return QExpr(s.qspace, terms)
end

function neq_qsum(s::QSum, where_defined::Vector{BitVector})
    _is_all_distinct(s) && return QExpr(s.expr.qspace, [s])
    combined_where = vecvec_or(which_ensemble_acting(s, do_abstract=true), where_defined)
    states = _expand_qsum_states(s, combined_where)
    terms = _collect_results(s, states)
    return QExpr(s.qspace, terms)
end
function neq(qeq::QExpr, do_abstract::Bool=false)::QExpr
    if length(qeq) == 0
        return qeq
    end
    if isa(qeq.terms[1], QSum) 
        out = neq_qsum(qeq.terms[1], do_abstract)
    else
        out = QExpr(qeq.qspace, neq(qeq.terms[1], do_abstract))
    end
    for t in qeq.terms[2:end]
        if isa(t, QSum)
            # expand this sum into distinct + diag parts
            out += neq_qsum(t, do_abstract)
        else
            out += neq(t, do_abstract)
        end
    end
    return out
end

#### where defined variants 
function neq(qeq::QExpr, where_defined::Vector{BitVector})::QExpr
    if length(qeq) == 0
        return qeq
    end
    if isa(qeq.terms[1], QSum) 
        out = neq_qsum(qeq.terms[1], where_defined)
    else
        out = QExpr(qeq.qspace, neq(qeq.terms[1], where_defined))
    end
    for t in qeq.terms[2:end]
        if isa(t, QSum)
            # expand this sum into distinct + diag parts
            out += neq_qsum(t, where_defined)
        else
            out += neq(t, where_defined)
        end
    end
    return out
end
function neq(q::QObj, where_defined::Vector{BitVector})::QObj
    return q
end
function neq(q::QAtomProduct, where_defined::Vector{BitVector})::QAtomProduct
    return q 
end
function neq(q::T, where_defined::Vector{BitVector})::T where {T<:QComposite}
    return only(modify_expr(q, neq(q.expr, where_defined)))
end
function neq(q::T, where_defined::Vector{BitVector})::T where {T<:QMultiComposite}
    return only(modify_expr(q, neq.(q.expr, where_defined)))
end

function neq(q::diffQEq)
    if !contains_abstract(q.left_hand_side)
        where_acting = which_ensemble_acting(q.left_hand_side)
        new_rhs = neq(q.expr, where_acting)
        return diffQEq(q.qspace, q.left_hand_side, new_rhs, Val(:nosimp))
    else
        error("Cannot neq a differential Equation with a QAbstract on the left hand side.")
    end
end
