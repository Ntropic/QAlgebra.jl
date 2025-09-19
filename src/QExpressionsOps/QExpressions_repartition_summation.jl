export decollision_QSum
# Repartitioning for collision free QSum indexing: in a  QSum:
# 0th initialize empty where_defined
# 1st create subsitutions in case of collision of current indexes with where_defined
# 2nd update where_defined  -> copy where defined, op_tuples, and var_tuples everywhere
# 3rd apply repartition to terms.

struct QSumDecollisionInds 
    init::Bool
    where_acting::Vector{BitVector}  
    op_tuples::Vector{Tuple{Int, Int}}
    var_tuples::Vector{Tuple{Int, Int}}
    var_inds::Vector{Int}
    function QSumDecollisionInds(init::Bool, where_acting::Vector{BitVector}, op_tuples::Vector{Tuple{Int, Int}}, var_tuples::Vector{Tuple{Int, Int}}, var_inds::Vector{Int})
        new(init, where_acting, op_tuples, var_tuples, var_inds) 
    end
    function QSumDecollisionInds(q::QSum)
        statespace = q.statespace
        subspace_info = statespace.subspace_info
        where_acting::Vector{BitVector} = [falses(n) for n in subspace_info.how_many_sum_by_ensemble]
        for subspace_ind in iter_all_indexes(q)
            ensemble_ind, summation_ind = Index2Ensemble_and_Summation(subspace_ind, subspace_info)
            where_acting[ensemble_ind][summation_ind] = true 
        end
        op_tuples = Vector{Tuple{Int, Int}}()
        var_tuples = Vector{Tuple{Int, Int}}()
        var_inds = Vector{Int}() #collect(1:length(statespace.params))
        new(false, where_acting, op_tuples, var_tuples, var_inds) 
    end
end 

# Detect summation index collisions, and find the first ensemble index that is free, to switch to, returns collision::Bool, the (potentially) updated where_acting, and a potentially new summation index, 
function collision_find_first_free(where_acting::Vector{BitVector}, subspace_ind::SubSpaceIndex, subspace_info::SubSpaceInfo)::Tuple{Bool, Vector{BitVector}, SubSpaceIndex}
    ensemble, summation = Index2Ensemble_and_Summation(subspace_ind, subspace_info)
    if where_acting[ensemble][summation] # found a collision 
        # find first free 
        free_sum_ind = findfirstfreeafterbefore(where_acting[ensemble], summation)
        @assert !isnothing(free_sum_ind) "No free summation index left in ensemble. "
        where_acting[ensemble][free_sum_ind] = true # QSums need to copy this! to not everright each other 
        return true, where_acting, SummationIndex2SubSpaceIndex(subspace_ind.outer, ensemble, free_sum_ind, subspace_info)
    else 
        where_acting[ensemble][summation] = true 
        return false, where_acting, subspace_ind
    end
end

function update_QSumDecollisionInds(q::QSum, d::QSumDecollisionInds)
    statespace = q.statespace
    subspace_info = statespace.subspace_info
    param_info = statespace.param_info

    new_op_tuples = Tuple{Int,Int}[]
    inds_tuples = Tuple{SubSpaceIndex,SubSpaceIndex}[]

    new_where = copy.(d.where_acting)

    @inbounds for (container, i, index) in iter_all_indexes_with_refs(q)
        collision, new_where, new_index = collision_find_first_free(new_where, index, subspace_info)
        if collision
            push!(new_op_tuples, (index.expanded, new_index.expanded))
            push!(inds_tuples, (index, new_index))
            container[i] = new_index
        end
    end

    if isempty(new_op_tuples)
        return QSumDecollisionInds(d.init, new_where, d.op_tuples, d.var_tuples, d.var_inds)
    end

    var_inds = collect(1:length(statespace.params))
    @inbounds for (index, new_index) in Base.Iterators.reverse(inds_tuples)
        curr_perm_params = map_by_subspace(index, new_index, param_info)
        var_inds = var_inds[curr_perm_params]
    end
    if d.init
        var_inds = var_inds[d.var_inds]
    end

    var_tuples = [(i, var_inds[i]) for i in eachindex(var_inds) if var_inds[i] != i]

    return QSumDecollisionInds(true, new_where, vcat(new_op_tuples, d.op_tuples), var_tuples, var_inds)
end

function decollision_QSum_product(q1::QSum, q2::QSum)::Vector{QComposite}
    # assumes that each QSum is already internally decollisioned! 
    statespace = q1.statespace
    subspace_info = statespace.subspace_info
    where_acting::Vector{BitVector} = which_summations_acting(q1, subspace_info) 
    decollision = QSumDecollisionInds(false, where_acting, Vector{Tuple{Int, Int}}(), Vector{Tuple{Int, Int}}(), Vector{Int}())
    # combine the two sums into one big sum
    base_terms  = QComposite[]
    nested_sums = QSum[]
    inner = decollision_QSum(q2, decollision)
    for t in inner
        if t isa QSum
            push!(nested_sums, t)
        else
            push!(base_terms, t)
        end
    end
    out_terms = QComposite[]
    if !isempty(base_terms)
        push!(out_terms, QSum(statespace, q1.expr*QExpr(statespace, base_terms), q1.eq_indexes, q1.neq_blocks))
    end
    for n in nested_sums
        merged_eq  = sort!(vcat(q1.eq_indexes, n.eq_indexes), by=expanded)
        merged_neq = vcat(q1.neq_blocks, n.neq_blocks)
        if length(merged_neq) > 1
            perm = sortperm(merged_neq; by = blk -> expanded(first(blk)))
            merged_neq = merged_neq[perm]
        end
        push!(out_terms, QSum(statespace, q1.expr*n.expr, merged_eq, merged_neq))
    end
    return out_terms
end


function decollision_QSum(q::QSum)::Vector{QComposite} 
    decollision = QSumDecollisionInds(q)
    return decollision_QSum(q, decollision, Val(:noupdate))
end
function decollision_QSum(q::QSum, decollision::QSumDecollisionInds, ::Val{:noupdate})::Vector{QComposite} #assume it is already updated 
    statespace = q.statespace
    base_terms  = QComposite[]
    nested_sums = QSum[]
    inner = decollision_QSum(q.expr, decollision)
    for t in inner.terms
        if t isa QSum
            push!(nested_sums, t)
        else
            push!(base_terms, t)
        end
    end
    out_terms = QComposite[]
    if !isempty(base_terms)
        push!(out_terms, QSum(statespace, QExpr(statespace, base_terms), q.eq_indexes, q.neq_blocks))
    end
    for n in nested_sums
        merged_eq  = sort!(vcat(q.eq_indexes, n.eq_indexes), by=expanded)
        merged_neq = vcat(q.neq_blocks, n.neq_blocks)
        if length(merged_neq) > 1
            perm = sortperm(merged_neq; by = blk -> expanded(first(blk)))
            merged_neq = merged_neq[perm]
        end
        push!(out_terms, QSum(statespace, n.expr, merged_eq, merged_neq))
    end
    return out_terms
end

function decollision_QSum(q::QSum, decollision::QSumDecollisionInds)::Vector{QComposite}
    decollision = update_QSumDecollisionInds(q, decollision)
    return decollision_QSum(q, decollision, Val(:noupdate))
end


function decollision_QSum(q::QTerm, decollision::QSumDecollisionInds, statespace::StateSpace)::QTerm
    op_indices = copy(q.op_indices)
    @inbounds for (old_ind, new_ind) in decollision.op_tuples
        @assert isnumeric(q, new_ind, statespace) "Cannot decollision QTerm, because new summation index is already in use, albeit undefined!"
        op_indices[old_ind], op_indices[new_ind] = op_indices[new_ind], op_indices[old_ind]
    end
    return QTerm(op_indices)
end
function decollision_QSum(q::QAbstract, decollision::QSumDecollisionInds, statespace::StateSpace)::QAbstract
    return add_to_index_map(q, decollision.op_tuples)
end
function decollision_QSum(q::QAtomProduct, decollision::QSumDecollisionInds)::Vector{QComposite}
    if decollision.init
        return QComposite[modify_coeff_expr(q, repartition(q.coeff_fun, decollision.var_tuples), QAtom[decollision_QSum(x, decollision, q.statespace) for x in q.expr])]
    else
        return QComposite[q] 
    end
end
function decollision_QSum(q::QExpr, decollision::QSumDecollisionInds)::QExpr
    new_terms = QComposite[]
    sizehint!(new_terms, length(q.terms))
     for t in q.terms
        append!(new_terms, decollision_QSum(t, decollision))
    end
    return QExpr(q.statespace, new_terms)
end
function decollision_QSum(q::T, decollision::QSumDecollisionInds)::Vector{QComposite} where T <: QComposite
    return QComposite[modify_coeff_expr(q, repartition(q.coeff_fun, decollision.var_tuples), decollision_QSum(q.expr, decollision))]
end
function decollision_QSum(q::T, decollision::QSumDecollisionInds)::Vector{QComposite} where T <: QMultiComposite
    return QComposite[modify_coeff_expr(q, repartition(q.coeff_fun, decollision.var_tuples), [decollision_QSum(qq, decollision) for qq in q.expr])]
end