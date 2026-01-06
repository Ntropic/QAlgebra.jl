export are_indices_defined, which_summations_acting, which_summations_to_root

const EnsembleBits = Vector{BitSet}

@inline _ensemble_count(qspace::QSpace) = length(qspace.subspace_info.where_ensembles)
@inline _empty_where_defined(qspace::QSpace)::EnsembleBits = [BitSet() for _ in 1:_ensemble_count(qspace)]
@inline _copy_where_defined(mask::EnsembleBits)::EnsembleBits = [copy(bs) for bs in mask]

@inline function _mark!(mask::EnsembleBits, ensemble::Int, slot::Int)
    (ensemble > 0 && slot > 0 && ensemble <= length(mask)) && push!(mask[ensemble], slot)
    return mask
end

function _missing_slots(defined::EnsembleBits, used::EnsembleBits)
    length(defined) == length(used) ||
        throw(ArgumentError("Ensemble dimensions mismatch ($(length(defined)) vs $(length(used)))."))
    missing = Vector{Vector{Int}}(undef, length(defined))
    @inbounds for i in eachindex(defined)
        missing[i] = [slot for slot in used[i] if !(slot in defined[i])]
    end
    return missing
end

function _missing_message(missing::Vector{Vector{Int}}, info::SubSpaceInfo)
    parts = String[]
    labels = info.ensemble_labels
    @inbounds for (i, slots) in enumerate(missing)
        isempty(slots) && continue
        name = i <= length(labels) ? join(labels[i], "/") : string("ensemble ", i)
        push!(parts, "$(name): " * join(slots, ","))
    end
    return isempty(parts) ? "[]" : join(parts, "; ")
end

# ========== Which ensembles are touched ===========================================================

function which_ensemble_acting(q::QObj; do_abstract::Bool=false)::EnsembleBits
    accum = _empty_where_defined(q.qspace)
    return which_ensemble_acting!(q, accum; do_abstract=do_abstract)
end
which_ensemble_acting(f::CFunctions.CFunction) = CFunctions.which_ensemble_acting(f)
which_ensemble_acting!(f::CFunctions.CFunction, accum::Vector{BitVector}) =
    CFunctions.which_ensemble_acting!(f, accum)
which_ensemble_acting(q::QAtom) = error("Cannot determine acting ensembles for a bare QAtom without its parent expression.")

function which_ensemble_acting!(q::QExpr, accum::EnsembleBits; do_abstract::Bool=false)::EnsembleBits
    for term in q.terms
        accum = which_ensemble_acting!(term, accum; do_abstract=do_abstract)
    end
    return accum
end

function which_ensemble_acting!(q::AbstractQSum, accum::EnsembleBits; do_abstract::Bool=false)::EnsembleBits
    accum = which_ensemble_acting!(q.expr, accum; do_abstract=do_abstract)
    info = q.qspace.subspace_info
    for idx in iter_all_indices(q)
        ensemble = info.ensemble_index_by_subspace_index[idx.outer]
        _mark!(accum, ensemble, idx.inner)
    end
    return accum
end

function which_ensemble_acting!(q::QComposite, accum::EnsembleBits; do_abstract::Bool=false)::EnsembleBits
    return which_ensemble_acting!(q.expr, accum; do_abstract=do_abstract)
end

function which_ensemble_acting!(q::QMultiComposite, accum::EnsembleBits; do_abstract::Bool=false)::EnsembleBits
    for expr in q.expr
        accum = which_ensemble_acting!(expr, accum; do_abstract=do_abstract)
    end
    return accum
end

function which_ensemble_acting!(q::QAtomProduct, accum::EnsembleBits; do_abstract::Bool=false)::EnsembleBits
    qspace = q.qspace
    for atom in q.expr
        accum = which_ensemble_acting!(atom, accum, qspace; do_abstract=do_abstract)
    end
    return accum
end

function which_ensemble_acting!(term::QTerm, accum::EnsembleBits, qspace::QSpace; do_abstract::Bool=false)::EnsembleBits
    subspaces = qspace.subspaces
    for particle in term.op_indices
        idx = particle.index
        ensemble = idx.ensemble
        ensemble > 0 || continue
        slot = idx.slot
        slot > 0 || continue
        subspace = subspaces[idx.subspace]
        acts = particle.operator != subspace.op_set.neutral_element
        acts || continue
        _mark!(accum, ensemble, slot)
    end
    return accum
end

function which_ensemble_acting!(abstract::QAbstract, accum::EnsembleBits, qspace::QSpace; do_abstract::Bool=false)
    do_abstract && return accum
    info = qspace.subspace_info
    for (idx1, idx2) in abstract.index_map
        ensemble1 = info.ensemble_index_by_subspace_index[idx1.outer]
        ensemble2 = info.ensemble_index_by_subspace_index[idx2.outer]
        _mark!(accum, ensemble1, idx1.inner)
        _mark!(accum, ensemble2, idx2.inner)
    end
    return accum
end
which_ensemble_acting!(atom::QAtom, accum::EnsembleBits, ::QSpace; do_abstract::Bool=false) = accum

# ========== Summation bookkeeping ================================================================

function which_summations_acting(q::QObj, subspace_info::SubSpaceInfo)::EnsembleBits
    empty_vec::EnsembleBits = [BitSet() for _ in 1:length(subspace_info.where_ensembles)]
    return which_summations_acting(q, empty_vec, subspace_info)
end
which_summations_acting(::QAtom, ::EnsembleBits, ::SubSpaceInfo) =
    error("Shouldn't call QAtom for which_summations_acting.")

which_summations_acting(q::QAtomProduct, where_acting::EnsembleBits, ::SubSpaceInfo) = where_acting

function which_summations_acting(q::QExpr, where_acting::EnsembleBits, info::SubSpaceInfo)::EnsembleBits
    for t in q.terms
        where_acting = which_summations_acting(t, where_acting, info)
    end
    return where_acting
end

function which_summations_acting(q::QComposite, where_acting::EnsembleBits, info::SubSpaceInfo)::EnsembleBits
    return which_summations_acting(q.expr, where_acting, info)
end

function which_summations_acting(q::QMultiComposite, where_acting::EnsembleBits, info::SubSpaceInfo)::EnsembleBits
    for x in q.expr
        where_acting = which_summations_acting(x, where_acting, info)
    end
    return where_acting
end

function which_summations_acting(q::AbstractQSum, where_acting::EnsembleBits, info::SubSpaceInfo)::EnsembleBits
    where_acting = which_summations_acting(q.expr, where_acting, info)
    @inbounds for index in iter_all_indices(q)
        ensemble = info.ensemble_index_by_subspace_index[index.outer]
        ensemble > 0 || continue
        if index.inner in where_acting[ensemble]
            error("Summation index $(index.inner) for ensemble $(ensemble) already defined in stack!")
        end
        _mark!(where_acting, ensemble, index.inner)
    end
    return where_acting
end

# gather a unique list of SubSpaceIndexes from nested AbstractQSums 
function gather_summation_indices!(q::T, curr_inds::Vector{SubSpaceIndex}=SubSpaceIndex[])::Vector{SubSpaceIndex} where T<: QComposite
    return gather_summation_indices!(q.expr, curr_inds)
end
function gather_summation_indices!(q::T, curr_inds::Vector{SubSpaceIndex}=SubSpaceIndex[])::Vector{SubSpaceIndex} where T<: QMultiComposite
    for t in q.expr 
        curr_inds = gather_summation_indices!(t, curr_inds)
    end
    return curr_inds 
end
function gather_summation_indices!(q::QExpr, curr_inds::Vector{SubSpaceIndex}=SubSpaceIndex[])::Vector{SubSpaceIndex}
    for t in q.terms 
        curr_inds = gather_summation_indices!(t, curr_inds)
    end
    return curr_inds
end
gather_summation_indices!(q::QAtomProduct, curr_inds::Vector{SubSpaceIndex}=SubSpaceIndex[])::Vector{SubSpaceIndex} = curr_inds
function gather_summation_indices!(q::AbstractQSum, curr_inds::Vector{SubSpaceIndex}=SubSpaceIndex[])::Vector{SubSpaceIndex}
    for index in iter_all_indices(q)
        if !any(isequal(index), curr_inds)
            push!(curr_inds, index)
        end
    end 
    return curr_inds
end

# ========== Well-definedness check ===============================================================

function are_indices_defined(q::QAtomProduct, where_defined::EnsembleBits)::Bool
    used = which_ensemble_acting(q)
    missing = _missing_slots(where_defined, used)
    if any(!isempty(s) for s in missing)
        msg = _missing_message(missing, q.qspace.subspace_info)
        error("Cannot use an undefined ensemble index on the right hand side of a differential equation! missing=$msg")
    end
    return true
end

function are_indices_defined(q::QExpr, where_defined::EnsembleBits)::Bool
    for t in q.terms
        are_indices_defined(t, where_defined) || return false
    end
    return true
end
are_indices_defined(q::QComposite, where_defined::EnsembleBits)::Bool = are_indices_defined(q.expr, where_defined)
function are_indices_defined(q::QMultiComposite, where_defined::EnsembleBits)::Bool
    for t in q.expr
        are_indices_defined(t, where_defined) || return false
    end
    return true
end

function are_indices_defined(q::AbstractQSum, where_defined::EnsembleBits)::Bool
    info = q.qspace.subspace_info
    new_where_defined = _copy_where_defined(where_defined)
    @inbounds for index in iter_all_indices(q)
        ensemble = info.ensemble_index_by_subspace_index[index.outer]
        ensemble > 0 || continue
        if index.inner in new_where_defined[ensemble]
            error("Summation index $(index.inner) for ensemble $(ensemble) already defined, cannot sum over defined indices!")
        end
        _mark!(new_where_defined, ensemble, index.inner)
    end

    return are_indices_defined(q.expr, new_where_defined)
end

function are_indices_defined(q::QExpr)::Bool
    where_defined::EnsembleBits = _empty_where_defined(q.qspace)
    return are_indices_defined(q, where_defined)
end
function are_indices_defined(q::diffQEq)::Bool
    defined = which_ensemble_acting(q.left_hand_side)
    return are_indices_defined(q.expr, defined)
end

# Legacy export placeholder – behaves like which_summations_acting for now.
which_summations_to_root(q::QObj) = which_summations_acting(q, q.qspace.subspace_info)
