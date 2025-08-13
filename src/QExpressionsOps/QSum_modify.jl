#### Flatten 
"""
flatten(qeq::QExpr) -> QExpr

Flattens nested Sums in quantum Equations (QExpr).
Does not support QSums within QComposites within QSums!
"""
function flatten(s::QSum, in_sum::Bool = false, in_sum_comp::Bool = false)
    s = copy(s)
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
        push!(out_terms, QSum(inner.statespace, QExpr(base_terms), s.indexes, s.subsystem_index, s.element_indexes, s.neq, Val(:simp)))
    end

    # for each nested sum, merge its indexes onto `s`'s
    for n in nested_sums
        dup = intersect(s.indexes, n.indexes)
        if !isempty(dup)
            error("Unsupported: duplicate summation indexes detected: $(dup)")
        end
        merged_idxs = vcat(s.indexes, n.indexes)
        if n.subsystem_index != s.subsystem_index
            error("Unsupported: nested sum with different subsystem index: $(n.subsystem_index)")
        end
        merged_einds = vcat(s.element_indexes, n.element_indexes)
        push!(out_terms, QSum(s.statespace, n.expr, merged_idxs, s.subsystem_index, merged_einds, s.neq, Val(:simp)))
    end

    return QExpr(inner.statespace, out_terms)
end 

function flatten(q::QAtomProduct, in_sum::Bool = false, in_sum_comp::Bool = false)
    return [copy(q)] 
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

# change from index1 to index2
function term_equal_indexes(expr, args...) # Base method to error
    throw(MethodError(term_equal_indexes, (typeof(expr), args...)))
end
# multiplies from the left 
function term_equal_indexes(term::QTerm, index1::Int, index2::Int, subspace::SubSpace)::Tuple{Bool, Vector{QTerm}, Vector{ComplexRational}}
    op1 = term.op_indices[index1]
    op2 = term.op_indices[index2]
    neutral = subspace.op_set.neutral_element
    if op1 === neutral && op2 === neutral
        return false, QTerm[term], ComplexRational[ComplexRational(1,0,1)]
    end
    results = subspace.op_set.op_product(op1, op2)
    new_terms = QTerm[]
    new_coeffs = ComplexRational[]
    for (coeff, op) in results
        new_term = copy(term)
        new_term.op_indices[index2] = op
        new_term.op_indices[index1] = neutral
        push!(new_terms, new_term)
        push!(new_coeffs, coeff)
    end
    return true, new_terms, new_coeffs
end 

function term_equal_indexes(abstract::QAbstract, index1::Int, index2::Int, subspace::SubSpace)::Tuple{Bool, Vector{QAbstract}, Vector{ComplexRational}}
    op_type = abstract.operator_type
    non_trivial_op_indices = op_type.non_trivial_op_indices
    if non_trivial_op_indices[index2]
        new_abstract = copy(abstract)
        push!(new_abstract.index_map, (index1, index2))
        return true, QAbstract[new_abstract], ComplexRational[ComplexRational(1,0,1)]
    end
    # append this rule to the index map 
    return false, QAbstract[abstract], ComplexRational[ComplexRational(1,0,1)]
end

function term_equal_indexes(q::QAtomProduct, index1::Int, index2::Int, subspace::SubSpace, coeff_inds1::Vector{Int}, coeff_inds2::Vector{Int})::Tuple{Bool, Vector{QAtomProduct}}
    changed_any = false
    term_variants = Vector{Vector{QAtom}}()
    coeff_variants = Vector{Vector{ComplexRational}}()
    for atom in q.expr
        changed, variants, coeffs = term_equal_indexes(atom, index1, index2, subspace)
        push!(term_variants, variants)
        push!(coeff_variants, coeffs)
        changed_any |= changed  # Check if any term was changed
    end
    changed, new_coeff_fun = CFunctions.term_equal_indexes(q.coeff_fun, coeff_inds1, coeff_inds2)
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

function term_equal_indexes(QExpr::QExpr, index1::Int, index2::Int, subspace::SubSpace, coeff_inds1::Vector{Int}, coeff_inds2::Vector{Int})::Tuple{Bool, Vector{QExpr}}
    changed_any = false
    elements = []
    for t in QExpr.terms
        changed, variants = term_equal_indexes(t, index1, index2, subspace, coeff_inds1, coeff_inds2)
        append!(elements, variants)
        changed_any |= changed  # Check if any term was changed
    end
    if !changed_any
        return false, [QExpr]
    end
    return true, [QExpr(QExpr.statespace, elements)]
end
#T <: QComposite case
function term_equal_indexes(q::T, index1::Int, index2::Int, subspace::SubSpace, coeff_inds1::Vector{Int}, coeff_inds2::Vector{Int})::Tuple{Bool, Vector{T}} where T<:QComposite
    changed, variants = term_equal_indexes(q.expr, index1, index2, subspace, coeff_inds1, coeff_inds2)
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
function term_equal_indexes(q::T, index1::Int, index2::Int, subspace::SubSpace, coeff_inds1::Vector{Int}, coeff_inds2::Vector{Int})::Tuple{Bool, Vector{T}} where T<:QMultiComposite
    changed, variants = term_equal_indexes(q.expr, index1, index2, subspace, coeff_inds1, coeff_inds2)
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
    neq(qeq::QExpr) -> QExpr

Transform sums into neq sums, where all indexes are different from each other, and returns a flattened QExpr with neq sums. 
Considers all cases of the sums, simplifying the cases in which indexes are the same, which then reduces the order of the sum (i.e. a sum_{j} x_i y_j => sum_{j} x_i y_j + im*z_i, where we used x_i*y_i=im*z_i).
"""
function neq(q::QObj)::QObj
    return q
end
function neq(q::QAtomProduct)::QAtomProduct
    return q 
end
function neq(q::T)::T where {T<:QComposite}
    return modify_expr(q, neq(q.expr))
end
function neq(q::T)::T where {T<:QMultiComposite}
    return modify_expr(q, neq.(q.expr))
end

function neq(qeq::QExpr)::QExpr
    # flatten first 
    qeq = flatten(qeq)
    if length(qeq) == 0
        return qeq
    end
    if isa(qeq.terms[1], QSum) 
        out = neq_qsum(qeq.terms[1])
    else
        out = QExpr(qeq.statespace, neq(qeq.terms[1]))
    end
    for t in qeq.terms[2:end]
        if isa(t, QSum)
            # expand this sum into distinct + diag parts
            out += neq_qsum(t)
        else
            out += neq(t)
        end
    end
    return out
end

# -------------------------------------------------------------------
# handle one QSum
function neq_qsum(s::QSum, index::Int=1)::QExpr
    if s.neq
        return QExpr(s.expr.statespace, [s])   # skip
    end
    n = length(s.element_indexes) # is at least 1
    if n < index
        error("neq: index $index is out of range for this QSum (with n=$n)")
    end
    # 1) the “all distinct” piece # add to the lower part and remove it here 

    # 2) we consider for each sum index combination all possible 
    ss = s.expr.statespace
    sub = ss.subspaces[s.subsystem_index]
    n_sub = length(sub.statespace_inds)
    # consider only one possible equality, then recursively process untill all possibilities have been checked
    curr_element = s.element_indexes[index]
    if index < n # recursively execute neq_qsum for higher possible indexes
        post_expr = neq_qsum(s, index + 1)
    else
        post_expr = QExpr(ss, QComposite[s])
    end
    pieces = copy(post_expr)

    # 3) now we assume index is equal to each of the parameters in subspace, with smaller index than the curr index of the sum 
    curr_ind_sum::Int = s.element_indexes[index]
    curr_statespace_ind::Int = sub.statespace_inds[curr_ind_sum]

    coeffs_of_subspace = ss.where_by_continuum[s.subsystem_index]
    curr_coeff_inds = [coeffs_of_subspace[i][curr_ind_sum] for i in 1:length(coeffs_of_subspace)]
    #println("\nnew run: ($index) : ", pieces) 
    for (new_ind_sum, new_statespace_sum) in zip(1:curr_ind_sum-1, sub.statespace_inds[1:curr_ind_sum-1])
        new_coeff_inds = [coeffs_of_subspace[i][new_ind_sum] for i in 1:length(coeffs_of_subspace)]
        new_statespace_ind = sub.statespace_inds[new_ind_sum]
        # check for each term in the subspace if curr_statespace_ind and new_statespace_ind are the neutral_element  
        for expr in post_expr.terms
           if isa(expr, QSum)
                #println("    ($index) - expr: ", expr)
                for t in expr.expr.terms
                    not_neutral, new_terms = term_equal_indexes(t, curr_statespace_ind, new_statespace_ind, sub, curr_coeff_inds, new_coeff_inds)
                    # add new terms as QSum(s) with corrected indexing 
                    if not_neutral
                        # remove expr.indexes[index] and similarly expr.element_indexes[index]
                        new_indexes = vcat(expr.indexes[1:index-1], expr.indexes[index+1:end])
                        new_element_indexes = vcat(expr.element_indexes[1:index-1], expr.element_indexes[index+1:end])
                        if length(new_indexes) == 0
                            for new_term in new_terms
                                pieces += new_term
                            end
                        else
                            pieces += QSum(s.statespace, QExpr(ss, new_terms, Val(:simp)), new_indexes, expr.subsystem_index, new_element_indexes, true)
                        end
                    else # no change to sum structure
                        pieces += QSum(s.statespace, QExpr(ss, new_terms, Val(:simp)), copy(expr.indexes), expr.subsystem_index, copy(expr.element_indexes), true)
                    end
                end
            else ## Old - no longer sufficient: if isa(expr, QTerm)
                not_neutral, new_terms = term_equal_indexes(expr, curr_statespace_ind, new_statespace_ind, sub, curr_coeff_inds, new_coeff_inds)
                if not_neutral
                    error("Unsupported: Element that isn't part of a Sum should no longer contain sum indexes")
                end
                for new_term in new_terms
                    pieces += new_term
                end
            end
        end
        #println("  Result for ($index => $curr_ind_sum, $new_ind_sum | $curr_statespace_ind, $new_statespace_ind):  " , pieces)
    end
    return pieces
end
