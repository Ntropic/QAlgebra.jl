export are_indexes_defined


cnimp(a::Bool, b::Bool) = b && !a
function converse_nonimplication(A::Vector{Vector{Bool}}, B::Vector{Vector{Bool}})
    # Assume they are equally shaped. 
    return broadcast.(cnimp, A, B) 
end

# CFunctions

import ..CFunctions: which_ensemble_acting
""" 
    which_ensemble_acting(q::QObj)::Vector{Vector{Bool}}

Returns a vector of vectors of booleans. Each inner vetor specifies which of its subsystem indexes are acted upon by the QObj. 
This includes actions from CFunctions. Th function should only be applied after substituting all QAbstract terms. 
Their present can be checked via `contains_abstract(q)`.
"""
function which_ensemble_acting(q::QAtom, subspace_info::SubSpaceInfo, neutral_ensembles_op::Vector{Vector{Is}}; do_abstract::Bool=false)::Vector{Vector{Bool}}
    my_ensembles::Vector{Vector{Bool}} = []
    for (inds, neutral) in zip(subspace_info.ensemble_indexes, neutral_ensembles_op)
        push!(my_ensembles, q.op_indices[inds] .!= neutral)
    end
    return my_ensembles
end
function which_ensemble_acting(q::QAbstract, subspace_info::SubSpaceInfo, neutral_ensembles_op::Vector{Vector{Is}}; do_abstract::Bool=false)
    if do_abstract ## ==> Assume instead that it is among the defined operator types | This is hacky, and probably not the best solution long term!
        return [zeros(Bool, n) for n in subspace_info.how_many_by_ensemble]
    else
        error("Which ensemble acting should be applied to abstractless expressions! ")
    end
end
function which_ensemble_acting(q::QAtomProduct; do_abstract::Bool=false)::Vector{Vector{Bool}}
    # xor between vectors of vector of bool 
    qspace = q.statespace
    return vecvec_or(reduce(vecvec_or, [which_ensemble_acting(t, qspace.subspace_info, qspace.I_ensemble_op; do_abstract=do_abstract) for t in q.expr]), 
                    which_ensemble_acting(q.coeff_fun))
end

function which_ensemble_acting(q::QExpr; do_abstract::Bool=false)::Vector{Vector{Bool}}
    return reduce(vecvec_or, [which_ensemble_acting(t, do_abstract=do_abstract) for t in q.terms])
end
function which_ensemble_acting(q::QSum; do_abstract::Bool=false)::Vector{Vector{Bool}}
    which_ensembles = which_ensemble_acting(q.expr, do_abstract=do_abstract)
    @inbounds @simd for index in q.indexes 
        which_ensembles[index.outer][Index2Ensemble(index)] = true
    end
    return which_ensembles
end

function which_ensemble_acting(q::QComposite; do_abstract::Bool=false)::Vector{Vector{Bool}}
    qspace = q.statespace
    return vecvec_or(which_ensemble_acting(q.expr, do_abstract=do_abstract), 
                    which_ensemble_acting(q.coeff_fun))
end
function which_ensemble_acting(q::QMultiComposite; do_abstract::Bool=false)::Vector{Vector{Bool}}
    qspace = q.statespace
    return vecvec_or(reduce(vecvec_or, [which_ensemble_acting(x, do_abstract=do_abstract) for x in q.expr]), 
                    which_ensemble_acting(q.coeff_fun))
end



""" 
    are_indexes_defined(q::diff_QEq)::Bool

Checks if all indexes n the differential equation are properly specified, either by the left-hand-side or 
by QSums on the right-hand-side.
"""
function are_indexes_defined(q::QAtomProduct, where_defined::Vector{Vector{Bool}})::Bool
    # check that no true on which_ensemble_acting, that isn't also a true on where_defined => converse nonimplication cnimp(a::Bool, b::Bool) = b && !a
    qspace = q.statespace
    if any(any, converse_nonimplication(where_defined, which_ensemble_acting(q)))
        x = converse_nonimplication(where_defined, which_ensemble_acting(q))
        error("Cannot use an undefined ensembles-index on the right hand side of a differential equation! $x")
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
    info = q.statespace.subspace_info
    indexes = q.indexes 
    new_where_defined = copy.(where_defined)
    for index in indexes 
        ensemble = Index2Ensemble(index, info)
        if new_where_defined[outer(index)][ensemble]
            index_str = Index2String(index, info)
            error("Summation index $index_str already defined, cannot sum over defined indexes!")
        end
        new_where_defined[outer(index)][ensemble]
    end 
    new_where_defined = copy.(where_defined)
    new_where_defined[outer_ind][element_indexes] .= true
    # check the summation QExpr 
    return are_indexes_defined(q.expr, new_where_defined)
end
function are_indexes_defined(q::QExpr)::Bool
    # check element wise if all indexes are defined
    where_defined::Vector{Vector{Bool}} = [zeros(Bool, n) for n in q.statespace.subspace_info.how_many_by_ensemble]
    return all([are_indexes_defined(t, where_defined) for t in q.terms])
end
function are_indexes_defined(q::diff_QEq)::Bool
    qspace = q.statespace
    # first two arguments for operators 
    # final argument for paramete/variables
    defined = which_ensemble_acting(q.left_hand_side)
    # check recursively if there are undefined elements in the right hand side
    return are_indexes_defined(q.expr, defined)
end
