# ===========================================================================================================================================================
# --------> Sorting  <---------------------------------------------------------------------------------------------------------------------------------------
# ===========================================================================================================================================================

function qAtom_sort_key(term::QTerm)
    # Here we convert var_exponents (a Vector{Int}) to a tuple so that it compares lexicographically.
    return tuple(0, term.op_indices...)
end
function qAtom_sort_key(term::QAbstract)
    return tuple(1, term.key_index, term.sub_index, term.exponent, Int(term.dag))
end

# Use your custom_sort_key for coefficients.
function qobj_sort_key(term::QAtomProduct)
    # Here we convert var_exponents (a Vector{Int}) to a tuple so that it compares lexicographically.
    return (tuple(0, 0, Int[], 0), sort_key(term.coeff_fun), qAtom_sort_key.(term.expr)...)
end
function qobj_sort_key(term::QSum)
    return (tuple(1, term.subsystem_index, length(term.element_indexes), term.element_indexes, length(term.expr)), term.neq)
end
function qobj_sort_key(term::QCompositeProduct)
    return (tuple(2, 0, Int[], 0), qAtom_sort_key.(term.expr)...) # first component specifies the type of object we are dealing with 
end
qobj_type_index(::Type{QExp}) = 3
qobj_type_index(::Type{QLog}) = 4 
qobj_type_index(::Type{QCommutator}) = 5 
function qobj_sort_key(term::T)  where T<:QComposite
    return (tuple(qobj_type_index(typeof(term)), 0, Int[], 0), qobj_sort_key.(term.expr)...) # first component specifies the type of object we are dealing with 
end
function qobj_sort_key(term::QPower)
    return (tuple(6, term.n, Int[], 0), qobj_sort_key.(term.expr)...) # first component specifies the type of object we are dealing with, second the subtype
end
function qobj_sort_key(term::QRoot)
    return (tuple(7, term.n, Int[], 0), qobj_sort_key.(term.expr)...) # first component specifies the type of object we are dealing with, second the subtype
end

# Sort the terms in a qExpr using the key above.
function sort(qeq::qExpr; kwargs...)
    # first the internal sort 
    terms = [sort(t) for t in qeq.terms]
    # now the outer sort 
    sorted_terms = _sort(terms, kwargs...)
    return qExpr(qeq.statespace, sorted_terms)
end
function sort(qprod::QAtomProduct; kwargs...)  # don't sort qprod 
    return copy(qprod)
end
function sort(qcomp::T; kwargs...) where T<:QMultiComposite # only sort recursively but not expr itself (outer most layer)
    new_qcomp = copy(qcomp)
    new_exprs::Vector{qExpr} = []
    for term in qcomp.expr 
        push!(new_exprs, sort(term))
    end 
    new_qcomp.expr = new_exprs
    return new_qcomp
end
function sort(qcomp::T; kwargs...) where T<:QComposite
    new_qcomp = copy(qcomp)
    new_qcomp.expr = sort(qcomp.expr; kwargs...)
    return new_qcomp
end

function _sort(qterms::AbstractVector{<:QComposite}; kwargs...)
    qterms2 = [sort(t) for t in qterms]  # recursive
    sorted_terms = Base.sort(qterms2; by= x-> qobj_sort_key(x), kwargs...)
    return sorted_terms
end