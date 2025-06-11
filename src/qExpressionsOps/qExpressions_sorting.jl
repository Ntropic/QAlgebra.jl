# ===========================================================================================================================================================
# --------> Sorting  <---------------------------------------------------------------------------------------------------------------------------------------
# ===========================================================================================================================================================

function qAtom_sort_key(term::qTerm)
    # Here we convert var_exponents (a Vector{Int}) to a tuple so that it compares lexicographically.
    return tuple(0, term.op_indices...)
end
function qAtom_sort_key(term::qAbstract)
    return tuple(1, term.key_index, term.subindex, term.exponent, Int(term.dag))
end

# Use your custom_sort_key for coefficients.
function qobj_sort_key(term::qAtomProduct)
    # Here we convert var_exponents (a Vector{Int}) to a tuple so that it compares lexicographically.
    return (tuple(0, Int[], 0), sort_key(term.coeff_fun), qAtom_sort_key.(term.expr)...)
end

function qobj_sort_key(term::qSum)
    curr_space = term.expr.statespace
    return (tuple(term.subsystem_index, term.element_indexes, length(term.expr)), term.neq)
end
# Sort the terms in a qExpr using the key above.
function sort(qeq::qExpr; kwargs...)
    # first the internal sort 
    terms = [sort(t) for t in qeq.terms]
    # now the outer sort 
    sorted_terms = _sort(terms, kwargs...)
    return qExpr(qeq.statespace, sorted_terms)
end
function sort(qprod::qAtomProduct; kwargs...)  # cannot sort qprod 
    return copy(qprod)
end
function sort(qcomp::qComposite; kwargs...)
    new_qcomp = copy(qcomp)
    new_qcomp.expr = _sort(qcomp.expr; kwargs...)
    return new_qcomp
end

function _sort(qterms::AbstractVector{<:qComposite}; kwargs...)
    qterms2 = [sort(t) for t in qterms]  # recursive
    sorted_terms = Base.sort(qterms2; by= x-> qobj_sort_key(x), kwargs...)
    return sorted_terms
end