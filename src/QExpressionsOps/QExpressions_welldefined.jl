export contains_abstract, which_continuum_acting, are_indexes_defined

""" 
    contains_abstract(q::QObj) -> Bool

Checks if the quantum object contains an abstract operator among its leaves.
"""
function contains_abstract(term::QExpr)::Bool
    return any([contains_abstract(t) for t in term.terms])
end
function contains_abstract(term::T)::Bool where T<:QComposite
    return contains_abstract(term.expr)
end
function contains_abstract(term::T)::Bool where T<:QMultiComposite
    return any([contains_abstract(t) for t in term.expr])
end
function contains_abstract(term::QAtomProduct)::Bool
    return any([isa(t, QAbstract) for t in term.expr])
end
function contains_abstract(term::diff_QEq)::Bool 
    return contains_abstract(term.expr) && contains_abstract(term.left_hand_side)
end

function vecvec_or(A::Vector{Vector{Bool}}, B::Vector{Vector{Bool}})
    # Assume they are equally shaped. 
    return broadcast.(|, A, B)
end
cnimp(a::Bool, b::Bool) = b && !a
function converse_nonimplication(A::Vector{Vector{Bool}}, B::Vector{Vector{Bool}})
    # Assume they are equally shaped. 
    return broadcast.(cnimp, A, B) 
end

# CFunctions
""" 
    which_continuum_acting(q::QObj)::Vector{Vector{Bool}}

Returns a vector of vectors of booleans. Each inner vetor specifies which of its subsystem indexes are acted upon by the QObj. 
This includes actions from CFunctions. Th function should only be applied after substituting all QAbstract terms. 
Their present can be checked via `contains_abstract(q)`.
"""
function which_continuum_acting(f::CAtom, where_continuums_f::Vector{Vector{Vector{Int}}})::Vector{Vector{Bool}}
    where_non_trivial::Vector{Vector{Bool}} = []
    for where_f in where_continuums_f
        push!(where_non_trivial, reduce(.|, [f.var_exponents[w] .!= 0 for w in where_f]))
    end
    return where_non_trivial
end
function which_continuum_acting(f::CSum, where_continuums_f::Vector{Vector{Vector{Int}}})::Vector{Vector{Bool}}
    # or of the individual terms 
    return reduce(vecvec_or, [which_continuum_acting(t, where_continuums_f) for t in f.terms])
end
function which_continuum_acting(f::CRational, where_continuums_f::Vector{Vector{Vector{Int}}})::Vector{Vector{Bool}}
    return vecvec_or(which_continuum_acting(f.numer, where_continuums_f), which_continuum_acting(f.denom, where_continuums_f))
end

# QObjs
function which_continuum_acting(q::QAtom, continuum_indexes::Vector{Vector{Int}}, neutral_continuums_op::Vector{Vector{Is}})::Vector{Vector{Bool}}
    my_continuums::Vector{Vector{Bool}} = []
    for (inds, neutral) in zip(continuum_indexes, neutral_continuums_op)
        push!(my_continuums, q.op_indices[inds] .!= neutral)
    end
    return my_continuums
end
function which_continuum_acting(q::QAbstract, continuum_indexes::Vector{Vector{Int}}, neutral_continuums_op::Vector{Vector{Is}})
    error("Which continuum acting should be applied to abstractless expressions!")
end
function which_continuum_acting(q::QAtomProduct)::Vector{Vector{Bool}}
    # xor between vectors of vector of bool 
    qspace = q.statespace
    return vecvec_or(reduce(vecvec_or, [which_continuum_acting(t, qspace.continuum_indexes, qspace.neutral_continuum_op) for t in q.expr]), which_continuum_acting(q.coeff_fun, qspace.where_by_continuum_var))
end

""" 
    are_indexes_defined(q::diff_QEq)::Bool

Checks if all indexes n the differential equation are properly specified, either by the left-hand-side or 
by QSums on the right-hand-side.
"""
function are_indexes_defined(q::QAtomProduct, where_defined::Vector{Vector{Bool}})::Bool
    # check that no true on which_continuum_acting, that isn't also a true on where_defined => converse nonimplication cnimp(a::Bool, b::Bool) = b && !a
    qspace = q.statespace
    if any(any, converse_nonimplication(where_defined, which_continuum_acting(q)))
        x = converse_nonimplication(where_defined, which_continuum_acting(q))
        error("Cannot use an undefined continuums-index on the right hand side of a differential equation! $x")
    end
    return true
end
function are_indexes_defined(q::QExpr, where_defined::Vector{Vector{Bool}})::Bool
    # check element wise if all indexes are defined
    return all([are_indexes_defined(t, where_defined) for t in q.terms])
end
function are_indexes_defined(q::T, where_defined::Vector{Vector{Bool}})::Bool where T <:QComposite
    return are_indexes_defined(q.expr, where_defined)
end
function are_indexes_defined(q::T, where_defined::Vector{Vector{Bool}})::Bool where T <:QMultiComposite
    return all([are_indexes_defined(t, where_defined) for t in q.exprs])
end
function are_indexes_defined(q::QSum, where_defined::Vector{Vector{Bool}})::Bool
    # add the QSum summation indexes 
    subsystem = q.subsystem_index 
    element_indexes = q.element_indexes
    # check where subsystem is among qspace.where_continuum
    outer_ind = findfirst(x -> x == subsystem, q.statespace.where_continuum)
    if outer_ind == nothing
        error("Subsystem $subsystem not found in qspace.subspace_by_ind!")
    end
    if !all(where_defined[outer_ind][element_indexes] .== false)
        error("Summation indexes already defined, cannot sum over defined indexes!")
    end
    new_where_defined = deepcopy(where_defined)
    new_where_defined[outer_ind][element_indexes] .= true
    # check the summation QExpr 
    return are_indexes_defined(q.expr, new_where_defined)
end
function are_indexes_defined(q::diff_QEq)::Bool
    qspace = q.statespace
    # first two arguments for operators 
    # final argument for paramete/variables
    defined = which_continuum_acting(q.left_hand_side)
    # check recursively if there are undefined elements in the right hand side
    return are_indexes_defined(q.expr, defined)
end
