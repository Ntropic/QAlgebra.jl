#### Flatten 
"""
flatten(qeq::qExpr) -> qExpr

Flattens nested Sums in quantum Equations (qExpr).
Does not support qSums within qComposites within qSums!
"""
<<<<<<< HEAD
function flatten(qeq::qExpr)::qExpr
    new_terms = qComposite[]
    for t in qeq.terms
        println(typeof(t))
        if t isa qSum
            flat_eq = flatten_qSum(t)
            append!(new_terms, flat_eq.terms)
        elseif t isa qAtomProduct 
            push!(new_terms, copy(t))
        elseif t isa qCompositeProduct
            push!(new_terms, copy(t))
        elseif t isa qMultiComposite
            new_t = copy(t) 
            terms = t.expr 
            flattened_terms = [flatten(s) for s in terms]
            new_t.expr = flattened_terms 
            push!(new_terms, new_t)
        elseif t isa qComposite
            new_t = copy(t) 
            new_t.expr = flatten(t.expr)
            push!(new_terms, new_t)
        else
            push!(new_terms, t)
        end
    end
    return qExpr(qeq.statespace, new_terms)
end

"""
flatten_qSum(s::qSum) -> qExpr

Take one `qSum` `s`.  First do `inner = flatten(s.expr)` so that all
deeper-nested sums are already one-level.  Split `inner.terms` into
• `base_terms` (just the `qTerm`’s)  
• `nested_sums` (any `qSum`’s).

Emit up to one “parent” sum over the `base_terms` (if non-empty), then
for each nested sum `n` emit a new `qSum` whose index-list is
`vcat(s.indexes, n.indexes)`.  Any duplicate index names will error.
"""
# No docstring here
function flatten_qSum(s::qSum)::qExpr
=======
function flatten(s::qSum, in_sum::Bool = false, in_sum_comp::Bool = false)
>>>>>>> b06581e (Updated neq, and abstract substitutions)
    # first, fully flatten the body
    if in_sum_comp && in_sum
        error("Unsupported: qSum found inside qComposite structure within an outer qSum.")
    end
    inner = flatten(s.expr, true, in_sum_comp)

    # pull out bare terms vs. sums
    base_terms::Vector{qComposite} = []
    nested_sums::Vector{qComposite} = []
    for t in inner.terms
        if t isa qSum
            push!(nested_sums, t)
        else
            push!(base_terms, t)
        end
    end

    out_terms = qComposite[]

    # if there were any base qTerms directly under `s`, keep a sum
    if !isempty(base_terms)
        push!(out_terms, qSum(inner.statespace, qExpr(base_terms), s.indexes, s.subsystem_index, s.element_indexes, s.neq))
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
        push!(out_terms, qSum(s.statespace, n.expr, merged_idxs, s.subsystem_index, merged_einds, s.neq))
    end
<<<<<<< HEAD
    return qExpr(inner.statespace, out_terms)
end 
=======
    return out_terms
end 
function flatten(q::qAtomProduct, in_sum::Bool = false, in_sum_comp::Bool = false)
    return [q] 
end
function flatten(q::qComposite, in_sum::Bool = false, in_sum_comp::Bool = false)
    new_q = copy(q)
    new_q.expr = flatten(q.expr)
    return [new_q]
end
function flatten(q::qMultiComposite, in_sum::Bool = false, in_sum_comp::Bool = false)
    new_q = copy(q)
    new_q.expr = [flatten(x) for x in q.expr]
    return [new_q]
end

function flatten(qeq::qExpr, in_sum::Bool = false, in_sum_comp::Bool = false)::qExpr
    new_terms = qComposite[]
    for s in qeq.terms
        append!(new_terms, flatten(s, in_sum, in_sum_comp))
    end
    return qExpr(qeq.statespace, new_terms)
end
>>>>>>> b06581e (Updated neq, and abstract substitutions)

# change from index1 to index2
function term_equal_indexes(expr, args...) # Base method to error
    throw(MethodError(term_equal_indexes, (typeof(expr), args...)))
end
# multiplies from the left 
function term_equal_indexes(term::qTerm, index1::Int, index2::Int, subspace::SubSpace)::Tuple{Bool, Vector{qTerm}}
    op1 = term.op_indices[index1]
    op2 = term.op_indices[index2]
    neutral = subspace.op_set.neutral_element
    if op1 === neutral && op2 === neutral
        return false, [term]
    end
    results = subspace.op_set.op_product(op1, op2)
    new_terms = qTerm[]
    for (coeff, op) in results
        new_term = deepcopy(term)
        new_term.op_indices[index2] = op
        new_term.op_indices[index1] = neutral
        push!(new_terms, new_term)
    end
    return true, new_terms
end 

function term_equal_indexes(abstract::qAbstract, index1::Int, index2::Int, subspace::SubSpace)::Tuple{Bool, Vector{qAbstract}}
    op_type = abstract.operator_type
    non_trivial_op_indices = op_type.non_trivial_op_indices
    if non_trivial_op_indices[index2]
        new_abstract = copy(abstract)
        push!(new_abstract.index_map, (index1, index2))
        return true, [new_abstract]
    end
    # append this rule to the index map 
    return false, [abstract]
end

function term_equal_indexes(prod::qAtomProduct, index1::Int, index2::Int, subspace::SubSpace, coeff_inds1::Vector{Int}, coeff_inds2::Vector{Int})::Tuple{Bool, Vector{qAtomProduct}}
    changed_any = false
    term_variants = Vector{Vector{qAtom}}()
    for atom in prod.expr
        changed, variants = term_equal_indexes(atom, index1, index2, subspace)
        push!(term_variants, variants)
        changed_any |= changed  # Check if any term was changed
    end
    changed, new_coeff_fun = FFunctions.term_equal_indexes(prod.coeff_fun, coeff_inds1, coeff_inds2)
    changed_any |= changed
    if !changed_any
        return false, [prod]
    end
    # Generate all combinations (cartesian product) of updated terms
    combinations = Iterators.product(term_variants...)
    simplified_products = qAtomProduct[]
    for combo in combinations
        new_expr = collect(combo)
        new_prod = qAtomProduct(prod.statespace, new_coeff_fun, new_expr)
        push!(simplified_products, new_prod)
    end
    return true, simplified_products
end

function term_equal_indexes(qExpr::qExpr, index1::Int, index2::Int, subspace::SubSpace, coeff_inds1::Vector{Int}, coeff_inds2::Vector{Int})::Tuple{Bool, Vector{qExpr}}
    changed_any = false
    elements = []
    for t in qExpr.terms
        changed, variants = term_equal_indexes(t, index1, index2, subspace, coeff_inds1, coeff_inds2)
        append!(elements, variants)
        changed_any |= changed  # Check if any term was changed
    end
    if !changed_any
        return false, [qExpr]
    end
    return true, [qExpr(qExpr.statespace, elements)]
end
#T <: qComposite case
function term_equal_indexes(q::T, index1::Int, index2::Int, subspace::SubSpace, coeff_inds1::Vector{Int}, coeff_inds2::Vector{Int})::Tuple{Bool, Vector{T}} where T<:qComposite
    changed, variants = term_equal_indexes(q.expr, index1, index2, subspace, coeff_inds1, coeff_inds2)
    if !changed
        return false, [q]
    end
    results::Vector{T} = []
    for v in variants
        q_new = copy(q)
        q_new.expr = v
        push!(results, q_new)
    end
    return true, results
end
#T <: qMultiComposite case
function term_equal_indexes(q::T, index1::Int, index2::Int, subspace::SubSpace, coeff_inds1::Vector{Int}, coeff_inds2::Vector{Int})::Tuple{Bool, Vector{T}} where T<:qMultiComposite
    changed, variants = term_equal_indexes(q.expr, index1, index2, subspace, coeff_inds1, coeff_inds2)
    if !changed
        return false, [q]
    end
    results::Vector{T} = []
    for v in variants
        q_new = copy(q)
        q_new.expr = v
        push!(results, q_new)
    end
    return true, results
end

"""
    neq(qeq::qExpr) -> qExpr

Transform sums into neq sums, where all indexes are different from each other, and returns a flattened qExpr with neq sums. 
Considers all cases of the sums, simplifying the cases in which indexes are the same, which then reduces the order of the sum (i.e. a sum_{j} x_i y_j => sum_{j} x_i y_j + im*z_i, where we used x_i*y_i=im*z_i).
"""
function neq(q::qObj)::qObj
    return q
end
function neq(q::qAtomProduct)::qAtomProduct
    return q 
end
function neq(q::T)::T where {T<:qComposite}
    q_copy = copy(q)
    q_copy.expr = neq(q.expr)
    return q_copy
end
function neq(q::T)::T where {T<:qMultiComposite}
    q_copy = copy(q)
    q_copy.expr = [neq(t) for t in q.expr]
    return q_copy
end

function neq(qeq::qExpr)::qExpr
    # flatten first 
    qeq = flatten(qeq)
    if length(qeq) == 0
        return qeq
    end
    if isa(qeq.terms[1], qSum) 
        out = neq_qsum(qeq.terms[1])
    else
        out = qExpr(qeq.statespace, neq(qeq.terms[1]))
    end
    for t in qeq.terms[2:end]
        if isa(t, qSum)
            # expand this sum into distinct + diag parts
            out += neq_qsum(t)
        else
            out += neq(t)
        end
    end
    for t in out.terms
        if t isa qSum
            t.neq = true
        end
    end
    return out
end

# -------------------------------------------------------------------
# handle one qSum
function neq_qsum(s::qSum, index::Int=1)::qExpr
    if s.neq
        return qExpr(s.expr.statespace, [s])   # skip
    end
    n = length(s.element_indexes) # is at least 1
    if n < index
        error("neq: index $index is out of range for this qSum (with n=$n)")
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
        post_expr = qExpr(ss, qComposite[s])
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
           if isa(expr, qSum)
                #println("    ($index) - expr: ", expr)
                for t in expr.expr.terms
                    not_neutral, new_terms = term_equal_indexes(t, curr_statespace_ind, new_statespace_ind, sub, curr_coeff_inds, new_coeff_inds)
                    # add new terms as qSum(s) with corrected indexing 
                    if not_neutral
                        # remove expr.indexes[index] and similarly expr.element_indexes[index]
                        new_indexes = vcat(expr.indexes[1:index-1], expr.indexes[index+1:end])
                        new_element_indexes = vcat(expr.element_indexes[1:index-1], expr.element_indexes[index+1:end])
                        if length(new_indexes) == 0
                            for new_term in new_terms
                                pieces += new_term
                            end
                        else
                            pieces += qSum(s.statespace, qExpr(ss, new_terms), new_indexes, expr.subsystem_index, new_element_indexes, expr.neq)
                        end
                    else # no change to sum structure
                        pieces += qSum(s.statespace, qExpr(ss, new_terms), copy(expr.indexes), expr.subsystem_index, copy(expr.element_indexes), expr.neq)
                    end
                end
            else ## Old - no longer sufficient: if isa(expr, qTerm)
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
