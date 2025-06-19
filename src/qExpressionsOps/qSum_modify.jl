#### Flatten 
"""
flatten(qeq::qExpr) -> qExpr

Flattens nested Sums in quantum Equations (qExpr).
"""
function flatten(qeq::qExpr)::qExpr
    new_terms = qComposite[]
    for t in qeq.terms
        if t isa qSum
            flat_eq = flatten_qSum(t)
            append!(new_terms, flat_eq.terms)
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
    # first, fully flatten the body
    inner = flatten(s.expr)

    # pull out bare terms vs. sums
    base_terms = [t for t in inner.terms if t isa qTerm]
    nested_sums = [t for t in inner.terms if t isa qSum]

    out_terms = qComposite[]

    # if there were any base qTerms directly under `s`, keep a sum
    if !isempty(base_terms)
        push!(out_terms,
            qSum(inner.statespace, qExpr(base_terms),
                s.indexes, s.subsystem_index, s.element_indexes, s.neq))
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

    return qExpr(inner.statespace, out_terms)
end

# change from index1 to index2
function term_equal_indexes(expr, args...)
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

function term_equal_indexes(prod::qAtomProduct, index1::Int, index2::Int, subspace::SubSpace)::Tuple{Bool, Vector{qAtomProduct}}
    changed_any = false
    term_variants = Vector{Vector{qAtom}}()
    for atom in prod.expr
        changed, variants = term_equal_indexes(atom, index1, index2, subspace)
        push!(term_variants, variants)
        changed_any |= changed  # Check if any term was changed
    end

    if !changed_any
        return false, [prod]
    end

    # Generate all combinations (cartesian product) of updated terms
    combinations = Iterators.product(term_variants...)

    simplified_products = qAtomProduct[]

    for combo in combinations
        new_expr = collect(combo)
        new_prod = qAtomProduct(prod.statespace, prod.coeff_fun, new_expr)
        push!(simplified_products, simplify(new_prod))
    end

    return true, simplified_products
end


"""
    neq(qeq::qExpr) -> qExpr

Transform sums into neq sums, where all indexes are different from each other, and returns a flattened qExpr with neq sums. 
Considers all cases of the sums, simplifying the cases in which indexes are the same, which then reduces the order of the sum (i.e. a sum_{j} x_i y_j => sum_{j} x_i y_j + im*z_i, where we used x_i*y_i=im*z_i).
"""
function neq(qeq::qExpr)::qExpr
    # flatten first 
    qeq = flatten(qeq)
    if length(qeq) == 0
        return qeq
    end
    if isa(qeq.terms[1], qSum) 
        out = neq_sum(qeq.terms[1])
    else
        out = qExpr(qeq.statespace, copy(qeq.terms[1]))
    end
    for t in qeq.terms[2:end]
        if isa(t, qSum)
            # expand this sum into distinct + diag parts
            out += neq_qsum(t)
        else
            out += copy(t)
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
            if isa(expr, qTerm)
                not_neutral, new_terms = term_equal_indexes(expr, curr_statespace_ind, new_statespace_ind, sub, curr_coeff_inds, new_coeff_inds)
                if not_neutral
                    error("Unsupported: Element that isn't part of a Sum should no longer contain sum indexes")
                end
                for new_term in new_terms
                    pieces += new_term
                end
            elseif isa(expr, qSum)
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
            end
        end
        #println("  Result for ($index => $curr_ind_sum, $new_ind_sum | $curr_statespace_ind, $new_statespace_ind):  " , pieces)
    end
    return simplify(pieces)
end
