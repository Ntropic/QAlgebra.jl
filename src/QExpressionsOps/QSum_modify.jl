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
flatsums(q::diff_QEq)::Bool = flatsums(q.expr)


""" 
    complexsums(q::QObj)::Bool 

Returns ture if any QSum is nested or in a QComposite. 
"""
complexsums(q::QSum, in_complex::Bool=false)::Bool = in_complex ? true : complexsums(q.expr, true)  
(complexsums(q::T, in_complex::Bool=false)::Bool) where {T<:QComposite} = complexsums(q.expr, true) 
complexsums(q::QAtomProduct, in_complex::Bool=false) = false 
(complexsums(q::T, in_complex::Bool=false)::Bool) where {T<:QMultiComposite} = any(complexsums.(q.expr, true))  
complexsums(q::QExpr, in_complex::Bool=false)::Bool = any(complexsums.(q.terms, in_complex))  
complexsums(q::diff_QEq)::Bool = complexsums(q.expr)



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
        new_prod = QAtomProduct(q.qspace, new_coeff_fun*factor, new_expr)
        push!(simplified_products, new_prod)
    end
    return true, simplified_products
end

function term_equal_indexes(qexpr::QExpr, index1::SubSpaceIndex, index2::SubSpaceIndex, subspace::SubSpace, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, Vector{QExpr}}
    changed_any = false
    elements = []
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

# -------------------------------------------------------------------
# handle one QSum
@inline function _sum_indexes_sorted(q::QSum)::Vector{SubSpaceIndex}
    idxs = SubSpaceIndex[]
    append!(idxs, q.eq_indexes)
    for block in q.neq_blocks
        append!(idxs, block)
    end
    sort!(idxs, by=expanded)
    return idxs
end

@inline function _block_membership(q::QSum)::Dict{Int,Int}
    membership = Dict{Int,Int}()
    for (blk_idx, block) in enumerate(q.neq_blocks)
        for idx in block
            membership[expanded(idx)] = blk_idx
        end
    end
    return membership
end

@inline function _can_coincide(membership::Dict{Int,Int}, idx1::SubSpaceIndex, idx2::SubSpaceIndex)::Bool
    blk1 = get(membership, expanded(idx1), nothing)
    blk2 = get(membership, expanded(idx2), nothing)
    return blk1 === nothing || blk2 === nothing || blk1 != blk2
end

function _clone_indexes(q::QSum)
    return copy(q.eq_indexes), [copy(block) for block in q.neq_blocks]
end

function _block_lt(a::Vector{SubSpaceIndex}, b::Vector{SubSpaceIndex})::Bool
    ea = expanded(a)
    eb = expanded(b)
    minlen = min(length(ea), length(eb))
    for i in 1:minlen
        if ea[i] != eb[i]
            return ea[i] < eb[i]
        end
    end
    return length(ea) < length(eb)
end

function _canonicalize_indexes!(eq_indexes::Vector{SubSpaceIndex}, neq_blocks::Vector{Vector{SubSpaceIndex}})
    sort!(eq_indexes, by=expanded)
    unique!(eq_indexes)
    i = 1
    while i <= length(neq_blocks)
        block = neq_blocks[i]
        block_sorted = sort(block, by=expanded)
        if length(block_sorted) <= 1
            if length(block_sorted) == 1
                push!(eq_indexes, block_sorted[1])
            end
            deleteat!(neq_blocks, i)
        else
            neq_blocks[i] = block_sorted
            i += 1
        end
    end
    sort!(eq_indexes, by=expanded)
    unique!(eq_indexes)
    sort!(neq_blocks, lt=_block_lt)
    return eq_indexes, neq_blocks
end

function _remove_index_from_qsum(q::QSum, idx::SubSpaceIndex)
    new_eq = SubSpaceIndex[]
    removed = false
    for existing in q.eq_indexes
        if existing == idx
            removed = true
        else
            push!(new_eq, existing)
        end
    end
    new_blocks = Vector{Vector{SubSpaceIndex}}()
    for block in q.neq_blocks
        if any(existing -> existing == idx, block)
            filtered = SubSpaceIndex[]
            for existing in block
                if existing == idx
                    removed = true
                else
                    push!(filtered, existing)
                end
            end
            if length(filtered) >= 2
                push!(new_blocks, filtered)
            elseif length(filtered) == 1
                push!(new_eq, filtered[1])
            end
        else
            push!(new_blocks, copy(block))
        end
    end
    removed || error("neq: index $(idx) not found in QSum during removal.")
    _canonicalize_indexes!(new_eq, new_blocks)
    return new_eq, new_blocks
end

function _enforce_all_distinct(q::QSum)::QSum
    idxs = collect(iter_all_indexes(q))
    isempty(idxs) && return q
    sorted_idxs = sort(idxs, by=expanded)
    return QSum(q.qspace, q.expr, SubSpaceIndex[], Vector{Vector{SubSpaceIndex}}([sorted_idxs]))
end

@inline function _is_all_distinct(q::QSum)::Bool
    isempty(q.eq_indexes) || return false
    length(q.neq_blocks) == 1 || return false
    return length(q.neq_blocks[1]) == length_all_indexes(q)
end

function neq_qsum(s::QSum, do_abstract::Bool=false)
    if _is_all_distinct(s)
        return QExpr(s.expr.qspace, [s])
    end
    ordered = _sum_indexes_sorted(s)
    membership = _block_membership(s)
    where_defined = which_ensemble_acting(s, do_abstract=do_abstract)
    return _neq_qsum_recursive(s, ordered, membership, 1, where_defined)
end

function neq_qsum(s::QSum, where_defined::Vector{BitVector})
    if _is_all_distinct(s)
        return QExpr(s.expr.qspace, [s])
    end
    ordered = _sum_indexes_sorted(s)
    membership = _block_membership(s)
    where_defined = vecvec_or(which_ensemble_acting(s, do_abstract=true), where_defined)
    return _neq_qsum_recursive(s, ordered, membership, 1, where_defined)
end

function _neq_qsum_recursive(s::QSum, ordered::Vector{SubSpaceIndex}, membership::Dict{Int,Int}, index::Int, where_defined::Vector{BitVector})::QExpr
    n = length(ordered)
    qspace = s.expr.qspace
    if n == 0
        return QExpr(qspace, QComposite[s])
    end
    post_expr = if index < n
        _neq_qsum_recursive(s, ordered, membership, index + 1, where_defined)
    else
        QExpr(qspace, QComposite[s])
    end
    curr_index = ordered[index]
    info = qspace.subspace_info
    curr_ensemble_index = Index2Ensemble(curr_index, info)
    check_subindexes = curr_index.inner == 1 ? Int[] : findall(where_defined[curr_ensemble_index][1:curr_index.inner-1])
    curr_subspace = qspace.subspaces[curr_index.outer]
    pieces = QExpr(qspace, [term isa QSum ? _enforce_all_distinct(term) : term for term in post_expr.terms])
    for new_ind_sum in check_subindexes
        new_statespace_sum = curr_subspace.ss_inner_ind[new_ind_sum]
        new_index = SubSpaceIndex(curr_index.outer, new_ind_sum, new_statespace_sum)
        _can_coincide(membership, curr_index, new_index) || continue
        curr_coeff_inds = changed_indices(map_by_subspace(curr_index, new_index, qspace.param_info))
        for expr in post_expr.terms
            if expr isa QSum
                qsum_expr = expr::QSum
                for t in qsum_expr.expr.terms
                    not_neutral, new_terms = term_equal_indexes(t, curr_index, new_index, curr_subspace, curr_coeff_inds)
                    new_inner_expr = QExpr(qspace, new_terms)
                    if not_neutral
                        new_eq, new_blocks = _remove_index_from_qsum(qsum_expr, curr_index)
                        if isempty(new_eq) && isempty(new_blocks)
                            for new_term in new_terms
                                pieces += new_term
                            end
                        else
                            tmp_qsum = QSum(qsum_expr.qspace, new_inner_expr, new_eq, new_blocks)
                            pieces += _enforce_all_distinct(tmp_qsum)
                        end
                    else
                        eq_copy, blocks_copy = _clone_indexes(qsum_expr)
                        tmp_qsum = QSum(qsum_expr.qspace, new_inner_expr, eq_copy, blocks_copy)
                        pieces += _enforce_all_distinct(tmp_qsum)
                    end
                end
            else
                not_neutral, new_terms = term_equal_indexes(expr, curr_index, new_index, curr_subspace, curr_coeff_inds)
                if not_neutral
                    error("Unsupported: Element that isn't part of a Sum should no longer contain sum indexes")
                end
                for new_term in new_terms
                    pieces += new_term
                end
            end
        end
    end
    return pieces
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

function neq(q::diff_QEq)
    if !contains_abstract(q.left_hand_side)
        where_acting = which_ensemble_acting(q.left_hand_side)
        new_rhs = neq(q.expr, where_acting)
        return diff_QEq(q.qspace, q.left_hand_side, new_rhs, q.do_braket)
    else
        error("Cannot neq a differential Equation with a QAbstract on the left hand side.")
    end
end
