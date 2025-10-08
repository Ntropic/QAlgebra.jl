export are_indexes_defined, which_summations_acting, which_summations_to_root


cnimp(a::Bool, b::Bool) = b && !a
function converse_nonimplication(A::Vector{BitVector}, B::Vector{BitVector})
    # Assume they are equally shaped. 
    return broadcast.(cnimp, A, B) 
end

# CFunctions

import ..CFunctions: which_ensemble_acting, which_ensemble_acting!

@inline function _empty_where_defined(qspace::QSpace)::Vector{BitVector}
    return [falses(n) for n in qspace.subspace_info.how_many_by_ensemble]
end
""" 
    which_ensemble_acting(q::QObj; do_abstract::Bool=false)::Vector{BitVector}

Returns a vector of vectors of booleans. Each inner vetor specifies which of its subsystem indexes are acted upon by the QObj. 
This includes actions from CFunctions. 
Will error on QAbstract if not using do_abstract=true.
"""
function which_ensemble_acting(q::QObj; do_abstract::Bool=false)::Vector{BitVector}
    accum = _empty_where_defined(q.qspace)
    return which_ensemble_acting!(q, accum; do_abstract=do_abstract) 
end
which_ensemble_acting(q::T) where T<: QAtom = error("Cannot extract qspace from QAtom, and hence not determine which ensembles are acting, without provided where_defined")

function which_ensemble_acting!(q::QExpr, accum::Vector{BitVector}; do_abstract::Bool=false)::Vector{BitVector}
    for term in q.terms
        accum = which_ensemble_acting!(term, accum, do_abstract=do_abstract)
    end
    return accum
end
function which_ensemble_acting!(q::AbstractQSum, accum::Vector{BitVector}; do_abstract::Bool=false)::Vector{BitVector}
    accum = which_ensemble_acting!(q.expr, accum; do_abstract=do_abstract)
    info = q.qspace.subspace_info
    for index in iter_all_indexes(q)
        accum[Index2Ensemble(index, info)][index.inner] = true
    end
    return accum
end
function which_ensemble_acting!(q::QComposite, accum::Vector{BitVector}; do_abstract::Bool=false)::Vector{BitVector}
    accum = which_ensemble_acting!(q.expr, accum, do_abstract=do_abstract)
    accum = which_ensemble_acting!(q.coeff_fun, accum)
    return accum
end
function which_ensemble_acting!(q::QMultiComposite, accum::Vector{BitVector}; do_abstract::Bool=false)::Vector{BitVector}
    for expr in q.expr
        accum = which_ensemble_acting!(expr, accum, do_abstract=do_abstract)
    end
    accum = which_ensemble_acting!(q.coeff_fun, accum)
    return accum
end

function which_ensemble_acting!(q::QAtomProduct, accum::Vector{BitVector}; do_abstract::Bool=false)::Vector{BitVector}
    # xor between vectors of vector of bool 
    qspace = q.qspace
    for atom in q.expr
        accum = which_ensemble_acting!(atom, accum, qspace.subspace_info, qspace.I_ensemble_op; do_abstract=do_abstract)
    end
    accum = which_ensemble_acting!(q.coeff_fun, accum)
    return accum
end
function which_ensemble_acting!(q::QAtom, accum::Vector{BitVector}, subspace_info::SubSpaceInfo, neutral_ensembles_op::Vector{Vector{Is}}; do_abstract::Bool=false)::Vector{BitVector}
    for (i, (inds, neutral)) in enumerate(zip(subspace_info.ensemble_indexes, neutral_ensembles_op))
        accum[i] .|=  q.op_indices[inds] .!= neutral
    end
    return accum
end
function which_ensemble_acting!(q::QAbstract, accum::Vector{BitVector}, subspace_info::SubSpaceInfo, neutral_ensembles_op::Vector{Vector{Is}}; do_abstract::Bool=false)
    if do_abstract
        return accum
    end
    if !isempty(q.index_map)
        for (idx1, idx2) in q.index_map
            for idx in (idx1, idx2)
                ensemble = Index2Ensemble(idx, subspace_info)
                accum[ensemble][idx.inner] = true
            end
        end
        return accum
    end
    return accum
end






""" 
    which_summations_acting(q::QObj)::Vector{BitVector}

Check, which Summation indexes are present. Returns a Boolean of 
"""
function which_summations_acting(q::QObj, subspace_info::SubSpaceInfo)::Vector{BitVector}
    empty_vec::Vector{BitVector} = [falses( s) for s in subspace_info.how_many_sum_by_ensemble]
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
function which_summations_acting(q::AbstractQSum, where_acting::Vector{BitVector}, which_indexes::Vector{Vector{Int}})::Vector{BitVector}
    which_summations_acting(q.expr, where_acting, which_indexes)
    info = q.qspace.subspace_info
    @inbounds for index in iter_all_indexes(q)
        ensemble, summation = Index2Ensemble_and_Summation(index, info)
        @assert !where_acting[ensemble][summation] "Summation index $(Index2String(index, info)) already defined in stack!"
        where_acting[ensemble][summation] = true
    end
    return where_acting
end

# gather a unique list of SubSpaceIndexes from nested AbstractQSums 
function gather_summation_indexes!(q::T, curr_inds::Vector{SubSpaceIndex}=SubSpaceIndex[])::Vector{SubSpaceIndex} where T<: QComposite
    return gather_summation_indexes!(q.expr, curr_inds)
end
function gather_summation_indexes!(q::T, curr_inds::Vector{SubSpaceIndex}=SubSpaceIndex[])::Vector{SubSpaceIndex} where T<: QMultiComposite
    for t in q.expr 
        curr_inds = gather_summation_indexes!(t, curr_inds)
    end
    return curr_inds 
end
function gather_summation_indexes!(q::QExpr, curr_inds::Vector{SubSpaceIndex}=SubSpaceIndex[])::Vector{SubSpaceIndex}
    for t in q.terms 
        curr_inds = gather_summation_indexes!(t, curr_inds)
    end
    return curr_inds
end
gather_summation_indexes!(q::QAtomProduct, curr_inds::Vector{SubSpaceIndex}=SubSpaceIndex[])::Vector{SubSpaceIndex} = curr_inds
function gather_summation_indexes!(q::AbstractQSum, curr_inds::Vector{SubSpaceIndex}=SubSpaceIndex[])::Vector{SubSpaceIndex}
    for index in iter_all_indexes(q)
        sorted_unique_push!(curr_inds, index)
    end 
    return curr_inds
end


""" 
    are_indexes_defined(q::diffQEq)::Bool

Checks if all indexes n the differential equation are properly specified, either by the left-hand-side or 
by AbstractQSums on the right-hand-side.
"""
function are_indexes_defined(q::QAtomProduct, where_defined::Vector{BitVector})::Bool
    # check that no true on which_ensemble_acting, that isn't also a true on where_defined => converse nonimplication cnimp(a::Bool, b::Bool) = b && !a
    qspace = q.qspace
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
function are_indexes_defined(q::AbstractQSum, where_defined::Vector{BitVector})::Bool
    info = q.qspace.subspace_info
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
    where_defined::Vector{BitVector} = [falses( n) for n in q.qspace.subspace_info.how_many_by_ensemble]
    return all([are_indexes_defined(t, where_defined) for t in q.terms])
end
function are_indexes_defined(q::diffQEq)::Bool
    qspace = q.qspace
    # first two arguments for operators 
    # final argument for paramete/variables
    defined = which_ensemble_acting(q.left_hand_side)
    # check recursively if there are undefined elements in the right hand side
    return are_indexes_defined(q.expr, defined)
end
