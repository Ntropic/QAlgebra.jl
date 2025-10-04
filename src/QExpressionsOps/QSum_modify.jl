#### Flatten 
export neq, flatsums, complexsums, Sum2Int

""" 
    flatsums(q::QObj)::Bool 

Checks if AbstractQSums are flattened in QObj. Returns true if there are no nested aggregators. 
"""
flatsums(q::AbstractQSum, has_sum::Bool=false)::Bool = has_sum ? false : flatsums(q.expr, has_sum) 
(flatsums(q::T, has_sum::Bool=false)::Bool) where {T<:QComposite} = flatsums(q.expr, has_sum) 
flatsums(q::QAtomProduct, has_sum::Bool=false) = true 
(flatsums(q::T, has_sum::Bool=false)::Bool) where {T<:QMultiComposite} = all(flatsums.(q.expr, has_sum))  
flatsums(q::QExpr, has_sum::Bool=false)::Bool = all(flatsums.(q.terms, has_sum))  
flatsums(q::diffQEq)::Bool = flatsums(q.expr)


""" 
    complexsums(q::QObj)::Bool 

Returns true if any AbstractQSum is nested or appears inside a QComposite. 
"""
complexsums(q::AbstractQSum, in_complex::Bool=false)::Bool = in_complex ? true : complexsums(q.expr, true)  
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

@inline function _is_all_distinct(q::AbstractQSum)::Bool
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
    new_blocks = copy.(blocks)
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

function _expand_qsum_states(q::AbstractQSum, where_defined::Vector{BitVector})::Vector{NeqState}
    initial_blocks = copy.(q.blocks)
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

function _collect_results(q::AbstractQSum, states::Vector{NeqState})::Vector{QComposite}
    results = QComposite[]
    for state in states
        if all(isempty(block.indexes) for block in state.blocks)
            append!(results, state.expr.terms)
        else
            distinct_blocks = _enforce_all_distinct(state.blocks)
            push!(results, QSum_like(q, state.expr, distinct_blocks))
        end
    end
    return results
end

function neq_qsum(s::AbstractQSum, do_abstract::Bool=false)
    _is_all_distinct(s) && return QExpr(s.expr.qspace, [s])
    where_defined = which_ensemble_acting(s, do_abstract=do_abstract)
    states = _expand_qsum_states(s, where_defined)
    terms = _collect_results(s, states)
    return QExpr(s.qspace, terms)
end

function neq_qsum(s::AbstractQSum, where_defined::Vector{BitVector})
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
    if isa(qeq.terms[1], AbstractQSum) 
        out = neq_qsum(qeq.terms[1], do_abstract)
    else
        out = QExpr(qeq.qspace, neq(qeq.terms[1], do_abstract))
    end
    for t in qeq.terms[2:end]
        if isa(t, AbstractQSum)
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
    if isa(qeq.terms[1], AbstractQSum) 
        out = neq_qsum(qeq.terms[1], where_defined)
    else
        out = QExpr(qeq.qspace, neq(qeq.terms[1], where_defined))
    end
    for t in qeq.terms[2:end]
        if isa(t, AbstractQSum)
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


##### Sum2Int <=================================================================================
@inline function which_continuums_ensembles(qspace::QSpace)::BitVector
    where_ensembles = qspace.subspace_info.where_ensembles
    init::BitVector = falses(length(where_ensembles))
    @inbounds for (i, curr_ensmble_index) in enumerate(where_ensembles)  
        ss = qspace.subspaces[curr_ensmble_index]
        init[i] = ss.as_continuum 
    end
    return init 
end

"""
    Sum2Int(q::QObj) -> QObj

Detects which ensembles in `q.qspace` are marked as continuous and rewrites
`QSum` terms accordingly: discrete-only sums stay unchanged, continuous-only
sums become `QInt`, and mixed sums are split into an integral and a remaining
discrete sum. Returns the transformed object.
"""
function Sum2Int(q::T)::T where T <: QObj 
    where_continuums = which_continuums_ensembles(q.qspace)
    return _Sum2Int(q, where_continuums)[1]
end
Sum2Int(q::T) where T <: QAtom = error("Cannot apply Sum2Int to QAtoms. ")

# --------------------------------------------------> walkers <----------------------------------------------------
function _Sum2Int(q::QExpr, where_continuums::BitVector)::Tuple{QExpr, Bool} 
    terms::Vector{QComposite} = q.terms 
    new_terms::Vector{QComposite} = Vector{QComposite}(undef, length(terms))
    any_changed::Bool = false
    @inbounds for (i, term) in enumerate(terms)
        new_term, changed = _Sum2Int(term, where_continuums)
        if changed 
            new_terms[i] = new_term 
            any_changed = true
        else
            new_terms[i] = term 
        end
    end
    if any_changed
        return QExpr(q.qspace, new_terms), any_changed 
    else
        return q, any_changed
    end
end
function _Sum2Int(q::T, where_continuums::BitVector)::Tuple{T, Bool} where T<:QComposite 
    new_expr, changed = _Sum2Int(q.expr, where_continuums)
    if changed
        return modify_expr(q, new_expr,Val(:nosimp))[1], changed 
    else
        return q, changed
    end
end
function _Sum2Int(q::T, where_continuums::BitVector)::Tuple{T, Bool} where T<:QMultiComposite 
    terms::Vector{QExpr} = q.expr 
    new_terms::Vector{QExpr} = Vector{QExpr}(undef, length(terms))
    any_changed::Bool = false
    @inbounds for (i, term) in enumerate(terms)
        new_term, changed = _Sum2Int(t, where_continuums)
        if changed 
            new_terms[i] = new_term 
            any_changed = true
        else
            new_terms[i] = term
        end
    end
    if any_changed
        return  modify_expr(q, new_terms,Val(:nosimp))[1], any_changed 
    else
        return q, any_changed
    end
end
function _Sum2Int(q::QCompositeProduct, where_continuums::BitVector)::Tuple{QCompositeProduct, Bool} 
    terms::Vector{QComposite} = q.expr 
    new_terms::Vector{QComposite} = Vector{QComposite}(undef, length(terms))
    any_changed::Bool = false
    @inbounds for (i, term) in enumerate(terms)
        new_term, changed = _Sum2Int(t, where_continuums)
        if changed 
            new_terms[i] = new_term 
            any_changed = true
        else
            new_terms[i] = term 
        end
    end
    if any_changed
        return  modify_expr(q, new_terms, Val(:nosimp))[1], any_changed 
    else
        return q, any_changed
    end
end
function _Sum2Int(q::QAtomProduct, ::BitVector)::Tuple{QAtomProduct, Bool}
     return q, false 
end

function _Sum2Int(q::QSum, where_continuums::BitVector)::Tuple{Union{QSum, QInt},Bool}
    blocks = q.blocks
    new_expr, changed = _Sum2Int(q.expr, where_continuums)
    has_int::Bool = false
    has_sum::Bool = false
    @inbounds for (block, is_cont) in zip(blocks, where_continuums)
        if length(block) > 0
            if is_cont
                has_int = true
            else
                has_sum = true
            end
            (has_int & has_sum) && break
        end
    end
    
    # Case 1: only discrete → unchanged
    if has_sum && !has_int
        if changed
            return _QSum_(q.qspace, new_expr, blocks), true
        end
        return q, false
    end

    # Case 2: only continuum → whole thing is a single QInt (copy blocks to avoid aliasing)
    if has_int && !has_sum
        copied_blocks = copy.(blocks)  # deep enough for your block type
        return _QInt_(q.qspace, q.expr, copied_blocks), true
    end

    n = length(blocks)
    int_blocks  = Vector{ConstrainedIndexBlock}(undef, n)
    sum_blocks = Vector{ConstrainedIndexBlock}(undef, n)
    @inbounds for (i, (block, is_cont)) in enumerate(zip(blocks, where_continuums))
        if is_cont
            int_blocks[i]  = copy(block)
            sum_blocks[i] = copy_empty(block)  # same shape, empty content
        else
            int_blocks[i]  = copy_empty(block)
            sum_blocks[i] = copy(block)
        end
    end
    return _QSum_(q.qspace, QExpr(q.qspace, QComposite[_QInt_(q.qspace, q.expr, int_blocks)]), sum_blocks), true
end