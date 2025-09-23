import ..CFunctions: repartition
using ..QSpaces: map_by_tindex
export reorder, reorder_full, reorder_time

struct IndexOrder
    op_order::Vector{Int}
    op_inverse::Vector{Int}
    var_moves::Vector{Tuple{Int, Int}}
    w_orders::Vector{Vector{Int}}
end

function IndexOrder(op_order::Vector{Int}, var_moves::Vector{Tuple{Int, Int}}, w_orders::Vector{Vector{Int}})
    return IndexOrder(op_order, invperm(op_order), var_moves, w_orders)
end

struct ReorderOrders
    base::IndexOrder
    full::IndexOrder
end

@inline permutation_moves(p::Vector{Int})::Vector{Tuple{Int, Int}} = [(i, pi) for (i, pi) in enumerate(p) if pi > i]

@inline function _collect_changed_pairs(mapvec::Vector{Int})::Vector{Tuple{Int, Int}}
    pairs = Tuple{Int,Int}[]
    @inbounds for i in eachindex(mapvec)
        val = mapvec[i]
        val == i && continue
        push!(pairs, (i, val))
    end
    return pairs
end

@inline function _time_order_data(qspace::QSpace, time_priority::Vector{Int})::Tuple{Vector{Int}, Vector{Tuple{Int,Int}}}
    n_t = qspace.max_t_ind + 1
    if n_t == 0
        return Int[], Tuple{Int,Int}[]
    end
    @assert length(time_priority) == n_t "Time priority vector must have length $(n_t)."
    perm = sortperm(time_priority, rev=true)
    time_inverse = invperm(perm)
    coeff_pairs = Tuple{Int,Int}[]
    param_info = qspace.param_info
    @inbounds for old_idx in 0:(n_t-1)
        priority = time_priority[old_idx+1]
        priority > 0 || continue
        new_idx = time_inverse[old_idx+1] - 1
        new_idx == old_idx && continue
        mapvec = map_by_tindex(new_idx, old_idx, param_info)
        append!(coeff_pairs, _collect_changed_pairs(mapvec))
    end
    return time_inverse, coeff_pairs
end

function where_defined_to_index_order_full(qspace::QSpace, where_defined::Vector{BitVector})::IndexOrder
    n_ops = length(qspace.I_op)
    n_vars = length(qspace.params)
    param_info = qspace.param_info
    subspace_info = qspace.subspace_info

    op_inds = collect(1:n_ops)
    var_inds = collect(1:n_vars)

    ensemble_indexes = subspace_info.ensemble_indexes
    where_ensembles = subspace_info.where_ensembles
    w_orders = Vector{Vector{Int}}(undef, length(where_defined))

    for (idx, (outer, w, c)) in enumerate(zip(where_ensembles, where_defined, ensemble_indexes))
        w_order = sortperm(w, rev=true)
        w_orders[idx] = w_order
        op_inds[c] = op_inds[c][w_order]
        for (from_idx, to_idx) in permutation_moves(w_order)
            idx_from = SubSpaceIndex(outer, from_idx, c[from_idx])
            idx_to = SubSpaceIndex(outer, to_idx, c[to_idx])
            curr_perm_params = map_by_subspace(idx_from, idx_to, param_info)
            var_inds = var_inds[curr_perm_params]
        end
    end

    var_moves = Tuple{Int,Int}[]
    @inbounds for i in eachindex(var_inds)
        if var_inds[i] != i
            push!(var_moves, (i, var_inds[i]))
        end
    end

    return IndexOrder(op_inds, var_moves, w_orders)
end

function where_defined_to_index_order(qspace::QSpace, where_defined::Vector{BitVector})::IndexOrder
    n_ops = length(qspace.I_op)
    n_vars = length(qspace.params)
    param_info = qspace.param_info
    subspace_info = qspace.subspace_info

    op_inds = collect(1:n_ops)
    var_inds = collect(1:n_vars)

    ensemble_indexes = subspace_info.ensemble_indexes
    where_ensembles = subspace_info.where_ensembles
    non_sum_counts = subspace_info.how_many_non_sum_by_ensemble
    sum_counts = subspace_info.how_many_sum_by_ensemble
    w_orders = Vector{Vector{Int}}(undef, length(where_defined))

    for (ensemble_idx, (outer, w, c)) in enumerate(zip(where_ensembles, where_defined, ensemble_indexes))
        non_count = non_sum_counts[ensemble_idx]
        sum_count = sum_counts[ensemble_idx]
        total = non_count + sum_count
        @assert total == length(w) "Mismatch between ensemble size and where_defined length."

        if non_count > 0
            non_flags = @view w[1:non_count]
            non_perm = sortperm(non_flags, rev=true)
            non_order = [i for i in non_perm]
        else
            non_order = Int[]
        end

        if sum_count > 0
            sum_offset = non_count
            sum_flags = @view w[sum_offset+1:total]
            sum_perm = sortperm(sum_flags, rev=true)
            sum_order = [sum_offset + i for i in sum_perm]
        else
            sum_order = Int[]
        end

        w_order = vcat(non_order, sum_order)
        if isempty(w_order)
            w_order = collect(1:total)
        end
        w_orders[ensemble_idx] = w_order

        op_inds[c] = op_inds[c][w_order]
        for (from_idx, to_idx) in permutation_moves(w_order)
            idx_from = SubSpaceIndex(outer, from_idx, c[from_idx])
            idx_to = SubSpaceIndex(outer, to_idx, c[to_idx])
            curr_perm_params = map_by_subspace(idx_from, idx_to, param_info)
            var_inds = var_inds[curr_perm_params]
        end
    end

    var_moves = Tuple{Int,Int}[]
    @inbounds for i in eachindex(var_inds)
        if var_inds[i] != i
            push!(var_moves, (i, var_inds[i]))
        end
    end

    return IndexOrder(op_inds, var_moves, w_orders)
end

@inline function build_reorder_orders(qspace::QSpace, where_defined::Vector{BitVector})::ReorderOrders
    base_order = where_defined_to_index_order(qspace, where_defined)
    full_order = where_defined_to_index_order_full(qspace, where_defined)
    return ReorderOrders(base_order, full_order)
end

@inline get_index_order(orders::ReorderOrders, ::Val{:base}) = orders.base
@inline get_index_order(orders::ReorderOrders, ::Val{:full}) = orders.full

@inline function remap_subspace_index(index::SubSpaceIndex, info::SubSpaceInfo, order::IndexOrder)
    new_expanded = order.op_inverse[index.expanded]
    new_index = SubSpaceIndex(new_expanded, info)
    @assert new_index.outer == index.outer "Reordering changed ensemble affiliation."
    return new_index
end

function reorder(q::QTerm, order::IndexOrder)::QTerm
    op_indices = q.op_indices[order.op_order]
    return QTerm(op_indices, q.time_index)
end

reorder(q::QObj, ::IndexOrder) = q

function reorder(q::QAtomProduct, mode::Val{M}, where_defined::Vector{BitVector}, orders::ReorderOrders; add_at_sum::Bool) where M
    order = get_index_order(orders, mode)
    new_coeff = repartition(q.coeff_fun, order.var_moves)
    new_atoms = QAtom[reorder(atom, order) for atom in q.expr]
    return modify_coeff_expr(q, new_coeff, new_atoms)
end

function reorder(q::QTerm, mode::Val{M}, where_defined::Vector{BitVector}, orders::ReorderOrders; add_at_sum::Bool) where M
    order = get_index_order(orders, mode)
    return reorder(q, order)
end

function reorder(q::QExpr, mode::Val{M}, where_defined::Vector{BitVector}, orders::ReorderOrders; add_at_sum::Bool) where M
    terms = QComposite[reorder(t, mode, where_defined, orders; add_at_sum=add_at_sum) for t in q.terms]
    return QExpr(q.qspace, terms)
end

function reorder(q::QComposite, mode::Val{M}, where_defined::Vector{BitVector}, orders::ReorderOrders; add_at_sum::Bool) where M
    q isa QSum && return reorder(q::QSum, mode, where_defined, orders; add_at_sum=add_at_sum)
    order = get_index_order(orders, mode)
    new_coeff = repartition(q.coeff_fun, order.var_moves)
    new_expr = reorder(q.expr, mode, where_defined, orders; add_at_sum=add_at_sum)
    return modify_coeff_expr(q, new_coeff, new_expr)
end

function reorder(q::QMultiComposite, mode::Val{M}, where_defined::Vector{BitVector}, orders::ReorderOrders; add_at_sum::Bool) where M
    order = get_index_order(orders, mode)
    new_coeff = repartition(q.coeff_fun, order.var_moves)
    new_expr = [reorder(qq, mode, where_defined, orders; add_at_sum=add_at_sum) for qq in q.expr]
    return modify_coeff_expr(q, new_coeff, new_expr)
end

function reorder(q::QSum, mode::Val{M}, where_defined::Vector{BitVector}, orders::ReorderOrders; add_at_sum::Bool) where M
    order = get_index_order(orders, mode)
    qspace = q.qspace
    info = qspace.subspace_info

    if !add_at_sum
        inner = reorder(q.expr, mode, where_defined, orders; add_at_sum=false)
        return modify_expr(q, inner, Val(:nodecollision))
    end

    new_where_defined = copy.(where_defined)
    for index in iter_all_indexes(q)
        ensemble = Index2Ensemble(index, info)
        if new_where_defined[ensemble][index.inner]
            index_str = Index2String(index, info)
            error("Summation index $index_str already defined, cannot sum over defined indexes!")
        end
        new_where_defined[ensemble][index.inner] = true
    end

    new_orders = build_reorder_orders(qspace, new_where_defined)
    new_order = get_index_order(new_orders, mode)

    inner = reorder(q.expr, mode, new_where_defined, new_orders; add_at_sum=true)

    new_eq = SubSpaceIndex[remap_subspace_index(index, info, new_order) for index in q.eq_indexes]
    sort!(new_eq, by=expanded)

    new_blocks = Vector{Vector{SubSpaceIndex}}()
    for blk in q.neq_blocks
        remapped = SubSpaceIndex[remap_subspace_index(index, info, new_order) for index in blk]
        sort!(remapped, by=expanded)
        push!(new_blocks, remapped)
    end

    return modify_expr_indexing(q, inner, new_eq, new_blocks, Val(:nodecollision))
end

function reorder(q::QObj, mode::Val{M}, where_defined::Vector{BitVector}, orders::ReorderOrders; add_at_sum::Bool) where M
    return q
end

function reorder(q::diff_QEq, mode::Val{M}) where M
    where_defined_lhs = which_ensemble_acting(q.left_hand_side)
    orders = build_reorder_orders(q.qspace, where_defined_lhs)
    order = get_index_order(orders, mode)

    lhs = q.left_hand_side
    if order.op_order != collect(1:length(order.op_order))
        lhs = reorder(lhs, mode, where_defined_lhs, orders; add_at_sum=false)
        where_defined_lhs = which_ensemble_acting(lhs)
        orders = build_reorder_orders(q.qspace, where_defined_lhs)
    end

    expr = reorder(q.expr, mode, where_defined_lhs, orders; add_at_sum=true)
    return diff_QEq(q.qspace, lhs, expr, Val(:nosimp); do_braket=q.do_braket)
end

"""
    reorder(eq::diff_QEq) -> diff_QEq

Reorder the ensemble indexes of `eq` so that already-defined (non-summation) indexes
stay on the left and remaining summation indexes are packed next to them. The
returned equation preserves the original structure, updating both the left-hand
side operator and the right-hand side expression.
"""
function reorder(q::diff_QEq)::diff_QEq
    return reorder(q, Val(:base))
end
"""
    reorder(expr::QExpr) -> QExpr

Reorder a quantum expression using the partial index ordering (non-summation
indexes first, followed by summation indexes within each ensemble). This keeps
summation indexes inside their dedicated slots while normalising operator order.
"""
function reorder(q::QExpr)::QExpr
    qspace = q.qspace
    where_defined = [falses(n) for n in qspace.subspace_info.how_many_by_ensemble]
    orders = build_reorder_orders(qspace, where_defined)
    return reorder(q, Val(:base), where_defined, orders; add_at_sum=true)
end
"""
    reorder_full(eq::diff_QEq) -> diff_QEq

Apply the legacy full reordering, which also shifts summation indexes out of
their dedicated block. This mirrors the original behaviour of `repartition` and
is only exposed for differential equations.
"""
function reorder_full(q::diff_QEq)::diff_QEq
    return reorder(q, Val(:full))
end

# --- time-index reordering --------------------------------------------------------------------------------------------

struct TimeReorderContext
    time_map::Vector{Int}
    coeff_moves::Vector{Tuple{Int,Int}}
end

@inline function _build_time_context(qspace::QSpace, time_priority::Vector{Int})
    isempty(time_priority) && return nothing
    any(>(0), time_priority) || return nothing
    time_inverse, coeff_moves = _time_order_data(qspace, time_priority)
    isempty(time_inverse) && isempty(coeff_moves) && return nothing
    unchanged = all(i -> time_inverse[i] == i, eachindex(time_inverse)) && isempty(coeff_moves)
    unchanged && return nothing
    time_map = [time_inverse[i] - 1 for i in eachindex(time_inverse)]
    return TimeReorderContext(time_map, coeff_moves)
end

@inline function _reorder_time(q::QObj, ::TimeReorderContext)
    return q
end

function _reorder_time(q::QTerm, ctx::TimeReorderContext)
    ti = q.time_index
    ti < 0 && return q
    ti + 1 > length(ctx.time_map) && return q
    new_ti = ctx.time_map[ti + 1]
    new_ti == ti && return q
    return modify_time_index(q, new_ti)
end

function _reorder_time(q::QAbstract, ctx::TimeReorderContext)
    ti = q.time_index
    ti < 0 && return q
    ti + 1 > length(ctx.time_map) && return q
    new_ti = ctx.time_map[ti + 1]
    new_ti == ti && return q
    return modify_time_index(q, new_ti)
end

function _reorder_time(vec::AbstractVector{T}, ctx::TimeReorderContext) where {T<:QObj}
    new_vec = Vector{T}(undef, length(vec))
    @inbounds for (i, item) in enumerate(vec)
        new_vec[i] = _reorder_time(item, ctx)::T
    end
    return new_vec
end

function _reorder_time(q::QComposite, ctx::TimeReorderContext)
    q isa QSum && return _reorder_time(q::QSum, ctx)
    new_expr = _reorder_time(q.expr, ctx)
    new_coeff = isempty(ctx.coeff_moves) ? q.coeff_fun : repartition(q.coeff_fun, ctx.coeff_moves)
    return modify_coeff_expr(q, new_coeff, new_expr)
end

function _reorder_time(q::QMultiComposite, ctx::TimeReorderContext)
    new_expr = _reorder_time(q.expr, ctx)
    new_coeff = isempty(ctx.coeff_moves) ? q.coeff_fun : repartition(q.coeff_fun, ctx.coeff_moves)
    return modify_coeff_expr(q, new_coeff, new_expr)
end

function _reorder_time(q::QSum, ctx::TimeReorderContext)
    inner = _reorder_time(q.expr, ctx)
    inner === q.expr && return q
    return modify_expr(q, inner, Val(:nodecollision))
end

function _reorder_time(q::QExpr, ctx::TimeReorderContext)
    new_terms = _reorder_time(q.terms, ctx)
    return QExpr(q.qspace, new_terms)
end

"""
    reorder_time(expr::QExpr) -> QExpr

Reorder the time indexes in `expr` so that used indexes are packed starting at 0,
mirroring the coefficient remapping performed by time substitutions.
"""
function reorder_time(q::QExpr)::QExpr
    qspace = q.qspace
    time_priority = Int.(contains_which_t_indexes(q))
    ctx = _build_time_context(qspace, time_priority)
    ctx === nothing && return q
    return _reorder_time(q, ctx)
end

"""
    reorder_time(eq::diff_QEq) -> diff_QEq

Reorder the time indexes in a differential equation, prioritising indexes already
present on the left-hand side before packing remaining indexes on the right-hand side.
The structure of the equation is preserved.
"""
function reorder_time(q::diff_QEq)::diff_QEq
    qspace = q.qspace
    lhs_usage = contains_which_t_indexes(q.left_hand_side)
    rhs_usage = contains_which_t_indexes(q.expr)
    time_priority = Vector{Int}(undef, length(lhs_usage)) # elements on both sides are prioritized
    @inbounds for i in eachindex(lhs_usage)
        lhs_val = lhs_usage[i] ? 2 : 0
        rhs_val = rhs_usage[i] ? 1 : 0
        time_priority[i] = lhs_val + rhs_val
    end
    ctx = _build_time_context(qspace, time_priority)
    ctx === nothing && return q
    new_lhs = _reorder_time(q.left_hand_side, ctx)
    new_rhs = _reorder_time(q.expr, ctx)
    return diff_QEq(qspace, new_lhs, new_rhs, Val(:nosimp); do_braket=q.do_braket)
end
