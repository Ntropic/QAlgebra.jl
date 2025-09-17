export are_indexes_defined, which_summations_acting, which_summations_to_root


cnimp(a::Bool, b::Bool) = b && !a
function converse_nonimplication(A::Vector{BitVector}, B::Vector{BitVector})
    # Assume they are equally shaped. 
    return broadcast.(cnimp, A, B) 
end

# CFunctions

import ..CFunctions: which_ensemble_acting
""" 
    which_ensemble_acting(q::QObj)::Vector{BitVector}

Returns a vector of vectors of booleans. Each inner vetor specifies which of its subsystem indexes are acted upon by the QObj. 
This includes actions from CFunctions. Th function should only be applied after substituting all QAbstract terms. 
Their present can be checked via `contains_abstract(q)`.
"""
function which_ensemble_acting(q::QAtom, subspace_info::SubSpaceInfo, neutral_ensembles_op::Vector{Vector{Is}}; do_abstract::Bool=false)::Vector{BitVector}
    my_ensembles::Vector{BitVector} = []
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
function which_ensemble_acting(q::QAtomProduct; do_abstract::Bool=false)::Vector{BitVector}
    # xor between vectors of vector of bool 
    qspace = q.statespace
    return vecvec_or!(reduce(vecvec_or!, [which_ensemble_acting(t, qspace.subspace_info, qspace.I_ensemble_op; do_abstract=do_abstract) for t in q.expr]), 
                    which_ensemble_acting(q.coeff_fun))
end
function which_ensemble_acting(q::QExpr; do_abstract::Bool=false)::Vector{BitVector}
    return reduce(vecvec_or, [which_ensemble_acting(t, do_abstract=do_abstract) for t in q.terms])
end
function which_ensemble_acting(q::QSum; do_abstract::Bool=false)::Vector{BitVector}
    which_ensembles = which_ensemble_acting(q.expr; do_abstract=do_abstract)
    info = q.statespace.subspace_info
    @inbounds @simd for index in iter_all_indexes(q)
        which_ensembles[Index2Ensemble(index, info)][index.inner] = true
    end
    return which_ensembles
end

function which_ensemble_acting(q::QComposite; do_abstract::Bool=false)::Vector{BitVector}
    return vecvec_or!(which_ensemble_acting(q.expr, do_abstract=do_abstract),  which_ensemble_acting(q.coeff_fun))
end
function which_ensemble_acting(q::QMultiComposite; do_abstract::Bool=false)::Vector{BitVector}
    return vecvec_or!(reduce(vecvec_or!, [which_ensemble_acting(x, do_abstract=do_abstract) for x in q.expr]), which_ensemble_acting(q.coeff_fun))
end

""" 
    which_summations_acting(q::QObj)::Vector{BitVector}

Check, which Summation indexes are present. Returns a Boolean of 
"""
function which_summations_acting(q::QObj, subspace_info::SubSpaceInfo)::Vector{BitVector}
    empty_vec::Vector{BitVector} = [zeros(Bool, s) for s in subspace_info.how_many_sum_by_ensemble]
    which_summations_acting(q, empty_vec, subspace_info.summation_indexes)
end
function which_summations_acting(q::QAtom, where_acting::Vector{BitVector}, ::Vector{Vector{Int}})
    error("Shouldn't call QAtom for which_summation_acting.")
end
function which_summations_acting(q::QAtomProduct, where_acting::Vector{BitVector}, which_indexes::Vector{Vector{Int}})
    return Vector{BitVector}()
end
function which_summations_acting(q::QExpr, where_acting::Vector{BitVector}, which_indexes::Vector{Vector{Int}})::Vector{BitVector}
    for t in q.terms
        vecvec_or!(where_acting, which_summations_acting(t, where_acting, which_indexes))
    end
    return where_acting
end
function which_summations_acting(q::T, where_acting::Vector{BitVector}, which_indexes::Vector{Vector{Int}})::Vector{BitVector} where T <: QComposite 
    return which_summations_acting(q.expr, where_acting, which_indexes)
end
function which_summations_acting(q::T, where_acting::Vector{BitVector}, which_indexes::Vector{Vector{Int}})::Vector{BitVector} where T <: QMultiComposite
    for x in q.expr
        vecvec_or!(where_acting, which_summations_acting(x,  where_acting, which_indexes))
    end
    return where_acting 
end
function which_summations_acting(q::QSum, where_acting::Vector{BitVector}, which_indexes::Vector{Vector{Int}})::Vector{BitVector}
    which_summations_acting(q.expr, where_acting, which_indexes)
    info = q.statespace.subspace_info
    @inbounds for index in iter_all_indexes(q)
        ensemble, summation = Index2Ensemble_and_Summation(index, info)
        @assert !where_acting[ensemble][summation] "Summation index $(Index2String(index, info)) already defined in stack!"
        where_acting[ensemble][summation] = true
    end
    return where_acting
end


"""
    which_summations_to_root(from::Union{QComposite,QExpr}) -> Vector{BitVector}

recursively navigate to the QObj root via `parent(x)`, collecting summation indices from any
`QSum` encountered on that path. Returns one BitVector per ensemble.
"""
function which_summations_to_root(from::T)::Vector{BitVector} where T <: Union{QComposite, QExpr}
    info = from.statespace.subspace_info
    where_acting = [zeros(Bool, s) for s in info.how_many_sum_by_ensemble]
    return which_summations_to_root(from, where_acting, info)
end

# Internal worker: keep walking parents; add indices when the current node is a QSum.
function which_summations_to_root(cur::Union{QComposite,QExpr,QParent}, where_acting::Vector{BitVector}, info::SubSpaceInfo)
    while true
        if cur isa QSum
            @inbounds for index in iter_all_indexes(cur)
                ensemble, summation = Index2Ensemble_and_Summation(index, info)
                @assert !where_acting[ensemble][summation] "Summation index already defined on path for $(Index2String(index, info))."
                where_acting[ensemble][summation] = true
            end
        end
        p = parent(cur)
        p === nothing && break   # found the root 
        cur = p
    end
    return where_acting
end

""" 
    are_indexes_defined(q::diff_QEq)::Bool

Checks if all indexes n the differential equation are properly specified, either by the left-hand-side or 
by QSums on the right-hand-side.
"""
function are_indexes_defined(q::QAtomProduct, where_defined::Vector{BitVector})::Bool
    # check that no true on which_ensemble_acting, that isn't also a true on where_defined => converse nonimplication cnimp(a::Bool, b::Bool) = b && !a
    qspace = q.statespace
    if any(any, converse_nonimplication(where_defined, which_ensemble_acting(q)))
        x = converse_nonimplication(where_defined, which_ensemble_acting(q))
        error("Cannot use an undefined ensembles-index on the right hand side of a differential equation! $x")
    end
    return true
end
function are_indexes_defined(q::QExpr, where_defined::Vector{BitVector})::Bool
    # check element wise if all indexes are defined
    return all([are_indexes_defined(t, where_defined) for t in q.terms])
end
function are_indexes_defined(q::T, where_defined::Vector{BitVector})::Bool where T <:QComposite
    return are_indexes_defined(q.expr, where_defined)
end
function are_indexes_defined(q::T, where_defined::Vector{BitVector})::Bool where T <:QMultiComposite
    return all([are_indexes_defined(t, where_defined) for t in q.exprs])
end
function are_indexes_defined(q::QSum, where_defined::Vector{BitVector})::Bool
    info = q.statespace.subspace_info
    new_where_defined = copy.(where_defined)
    @inbounds for index in iter_all_indexes(q)
        ensemble = Index2Ensemble(index, info)
        if new_where_defined[ensemble][index.inner]
            index_str = Index2String(index, info)
            error("Summation index $index_str already defined, cannot sum over defined indexes!")
        end
        new_where_defined[ensemble][index.inner] = true
    end

    return are_indexes_defined(q.expr, new_where_defined)
end

function are_indexes_defined(q::QExpr)::Bool
    # check element wise if all indexes are defined
    where_defined::Vector{BitVector} = [zeros(Bool, n) for n in q.statespace.subspace_info.how_many_by_ensemble]
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
