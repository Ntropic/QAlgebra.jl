#### Flatten 
import Base.Iterators: flatten
export neq

"""
flatten(qeq::QExpr) -> QExpr

Flattens nested Sums in quantum Equations (QExpr).
Does not support QSums within QComposites within QSums!
"""
function flatten(s::QSum, in_sum::Bool = false, in_sum_comp::Bool = false)
    # first, fully flatten the body
    if in_sum_comp && in_sum
        error("Unsupported: QSum found inside QComposite structure within an outer QSum.")
    end
    inner = flatten(s.expr, true, in_sum_comp)

    # pull out bare terms vs. sums
    base_terms::Vector{QComposite} = []
    nested_sums::Vector{QComposite} = []
    for t in inner.terms
        if t isa QSum
            push!(nested_sums, t)
        else
            push!(base_terms, t)
        end
    end

    out_terms = QComposite[]

    # if there were any base qTerms directly under `s`, keep a sum
    if !isempty(base_terms)
        push!(out_terms, QSum(QExpr(base_terms), s.indexes, s.neq))
    end

    # for each nested sum, merge its indexes onto `s`'s
    for n in nested_sums
        dup = intersect(s.indexes, n.indexes)
        if !isempty(dup)
            error("Unsupported: duplicate summation indexes detected: $(dup)")
        end
        merged_inds = vcat(s.indexes, n.indexes)
        push!(out_terms, QSum(n.expr, merged_inds, s.neq))
    end

    return QExpr(inner.statespace, out_terms)
end 

function flatten(q::QAtomProduct, in_sum::Bool = false, in_sum_comp::Bool = false)
    return [q] 
end
function flatten(q::T, in_sum::Bool = false, in_sum_comp::Bool = false) where T<:QComposite
    return [modify_expr(q, flatten(q.expr))]
end
function flatten(q::QMultiComposite, in_sum::Bool = false, in_sum_comp::Bool = false)
    return [modify_expr(q, flatten.(q.expr))]
end

function flatten(qeq::QExpr, in_sum::Bool = false, in_sum_comp::Bool = false)::QExpr
    new_terms = QComposite[]
    for s in qeq.terms
        append!(new_terms, flatten(s, in_sum, in_sum_comp))
    end
    return QExpr(qeq.statespace, new_terms)
end

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
        push!(new_terms, QTerm(op_indices, Val(:nocopy)))
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
        new_prod = QAtomProduct(q.statespace, new_coeff_fun*factor, new_expr)
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
    return true, [QExpr(qexpr.statespace, elements)]
end
#T <: QComposite case
function term_equal_indexes(q::T, index1::SubSpaceIndex, index2::SubSpaceIndex, subspace::SubSpace, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, Vector{T}} where T<:QComposite
    changed, variants = term_equal_indexes(q.expr, index1, index2, subspace, coeff_ind_order)
    if !changed
        return false, [q]
    end
    results::Vector{T} = []
    for v in variants
        push!(results, modify_expr(q, v))
    end
    return true, results
end
#T <: QMultiComposite case
function term_equal_indexes(q::T, index1::SubSpaceIndex, index2::SubSpaceIndex, subspace::SubSpace, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, Vector{T}} where T<:QMultiComposite
    changed, variants = term_equal_indexes(q.expr, index1, index2, subspace, coeff_ind_order)
    if !changed
        return false, [q]
    end
    results::Vector{T} = []
    for v in variants
        push!(results, modify_expr(q, v))
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
    return modify_expr(q, neq(q.expr, do_abstract))
end
function neq(q::T, do_abstract::Bool=false)::T where {T<:QMultiComposite}
    return modify_expr(q, neq.(q.expr, do_abstract))
end

# -------------------------------------------------------------------
# handle one QSum
function neq_qsum(s::QSum, do_abstract::Bool=false)
    # determine where defined 
    where_defined = which_ensemble_acting(s, do_abstract=do_abstract) 
    return neq_qsum(s, 1, where_defined )
end
function neq_qsum(s::QSum, where_defined::Vector{Vector{Bool}})
    # determine where defined 
    where_defined = vecvec_or(which_ensemble_acting(s, do_abstract=true) , where_defined )
    return neq_qsum(s, 1, where_defined )
end
function neq_qsum(s::QSum, index::Int, where_defined::Vector{Vector{Bool}})::QExpr
    if s.neq
        return QExpr(s.expr.statespace, [s])   # skip
    end
    n = length(s.indexes) # is at least 1
    if n < index
        error("neq: index $index is out of range for this QSum (with n=$n)")
    end

    statespace = s.expr.statespace
    info = statespace.subspace_info
    curr_index::SubSpaceIndex = s.indexes[index]
    curr_ensemble_index = Index2Ensemble(curr_index, info)
    check_subindexes = findall(where_defined[curr_ensemble_index][1:curr_index.inner-1])
    curr_subspace = statespace.subspaces[curr_index.outer] # subspace 
    curr_statespace_ind = curr_subspace.ss_inner_ind[curr_index.inner] # should be the same as curr_index.expanded 

    # consider only one possible equality, then recursively process untill all possibilities have been checked
    if index < n # recursively execute neq_qsum for higher possible indexes
        post_expr = neq_qsum(s, index + 1, where_defined)
    else
        post_expr = QExpr(statespace, QComposite[s])
    end
    pieces = QExpr(statespace, [isa(t, QSum) ? QSum(t.expr, t.indexes, true) : t for t in post_expr.terms])

    for (new_ind_sum, new_statespace_sum) in zip(check_subindexes, curr_subspace.ss_inner_ind[check_subindexes])
        new_index = SubSpaceIndex(curr_index.outer, new_ind_sum, new_statespace_sum)   # one way to do it that doesn't require accessing subspace_info
        
        curr_coeff_inds = changed_indices(map_by_subspace(curr_index, new_index, statespace.param_info))
        new_statespace_ind = curr_subspace.ss_inner_ind[new_ind_sum]

        # check for each term in the subspace if curr_statespace_ind and new_statespace_ind are the neutral_element  
        for expr in post_expr.terms
           if isa(expr, QSum)
                for t in expr.expr.terms
                    not_neutral, new_terms = term_equal_indexes(t, curr_index, new_index, curr_subspace, curr_coeff_inds)
                    # add new terms as QSum(s) with corrected indexing 
                    if not_neutral # something got changed?
                        new_indexes = vcat(expr.indexes[1:index-1], expr.indexes[index+1:end])
                        if length(new_indexes) == 0
                            for new_term in new_terms
                                pieces += new_term
                            end
                        else
                            pieces += QSum(QExpr(statespace, new_terms, Val(:simp)), new_indexes, true)
                        end
                    else # no change to sum structure
                        pieces += QSum(QExpr(statespace, new_terms, Val(:simp)), expr.indexes, true)
                    end
                end
            else ## Old - no longer sufficient: if isa(expr, QTerm)
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
    # flatten first 
    qeq = flatten(qeq)
    if length(qeq) == 0
        return qeq
    end
    if isa(qeq.terms[1], QSum) 
        out = neq_qsum(qeq.terms[1], do_abstract)
    else
        out = QExpr(qeq.statespace, neq(qeq.terms[1], do_abstract))
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
function neq(qeq::QExpr, where_defined::Vector{Vector{Bool}})::QExpr
    # flatten first 
    qeq = flatten(qeq)
    if length(qeq) == 0
        return qeq
    end
    if isa(qeq.terms[1], QSum) 
        out = neq_qsum(qeq.terms[1], 1, where_defined)
    else
        out = QExpr(qeq.statespace, neq(qeq.terms[1], where_defined))
    end
    for t in qeq.terms[2:end]
        if isa(t, QSum)
            # expand this sum into distinct + diag parts
            out += neq_qsum(t, 1, where_defined)
        else
            out += neq(t, where_defined)
        end
    end
    return out
end
function neq(q::QObj, where_defined::Vector{Vector{Bool}})::QObj
    return q
end
function neq(q::QAtomProduct, where_defined::Vector{Vector{Bool}})::QAtomProduct
    return q 
end
function neq(q::T, where_defined::Vector{Vector{Bool}})::T where {T<:QComposite}
    return modify_expr(q, neq(q.expr, where_defined))
end
function neq(q::T, where_defined::Vector{Vector{Bool}})::T where {T<:QMultiComposite}
    return modify_expr(q, neq.(q.expr, where_defined))
end

function neq(q::diff_QEq)
    if !contains_abstract(q.left_hand_side)
        where_acting = which_ensemble_acting(q.left_hand_side)
        new_rhs = neq(q.expr, where_acting)
        return diff_QEq(q.statespace, q.left_hand_side, new_rhs, q.do_braket)
    else
        error("Cannot neq a differential Equation with a QAbstract on the left hand side.")
    end
end