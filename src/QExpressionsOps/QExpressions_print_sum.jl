function _sum_index_subscript(indexes::Vector{SubSpaceIndex}, info::SubSpaceInfo; do_latex::Bool=false)
    isempty(indexes) && return do_latex ? "\\sum" : "∑"
    labels = Index2String.(indexes, Ref(info))
    if do_latex
        body = join(labels, ",")
        return "\\sum_{" * body * "}"
    else
        seq = join(labels, ",")
        sub = length(labels) == 1 ? str2sub(seq) : str2sub("(" * seq * ")")
        return "∑" * sub
    end
end
function _format_superscript(content::AbstractString; do_latex::Bool)
    isempty(content) && return ""
    if do_latex
        return "^{" * content * "}"
    else
        return str2sup( content )
    end
end
function _format_neq_condition(condition::NeqConstraint{SubSpaceIndex}, subspace_info::SubSpaceInfo; do_latex::Bool=false)
    lhs_str = Index2String(condition.lhs, subspace_info)
    rhs_str = Index2String(condition.rhs, subspace_info)
    if do_latex
        return lhs_str * raw" \neq " * rhs_str
    else 
        return lhs_str * " ≠ " * rhs_str
    end
end

function eq_counter_and_neq_indexes_by_block(block::ConstrainedIndexBlock, where_acting_block::BitVector, subspace_info::SubSpaceInfo)::Tuple{Int, Vector{NeqConstraint{SubSpaceIndex}}}
    non_sum = block.how_many_non_sum
    neq_constraints::Vector{NeqConstraint{SubSpaceIndex}} = Vector{NeqConstraint{SubSpaceIndex}}()
    eq_counter = 0
    for (ind, constraint) in zip(block.indexes, block.constraints)
        curr_inner = ind.inner
        for i in vcat(1:non_sum, curr_inner+1:length(where_acting_block))
            if where_acting_block[i] 
                if !constraint[i]
                    push!(neq_constraints, neq(ind, SubSpaceIndex(ind.outer, i, subspace_info)))
                else
                    eq_counter += 1
                end
            end
        end
    end
    return (eq_counter, neq_constraints)
end
function eq_counter_and_neq_indexes(blocks::Vector{ConstrainedIndexBlock}, where_acting::Vector{BitVector}, subspace_info::SubSpaceInfo)::Tuple{Int, Vector{NeqConstraint{SubSpaceIndex}}}
    eq_counter = 0
    neq_constraints::Vector{NeqConstraint{SubSpaceIndex}} = Vector{NeqConstraint{SubSpaceIndex}}()
    for (block, where_acting_block) in zip(blocks, where_acting)
        new_count, new_inds = eq_counter_and_neq_indexes_by_block(block, where_acting_block, subspace_info)
        eq_counter += new_count
        append!(neq_constraints, new_inds)
    end 
    return (eq_counter, neq_constraints)
end

function sum_symbol_str(term::QSum, where_acting::Vector{BitVector}; do_latex::Bool=false)
    subspace_info = term.qspace.subspace_info
    base = _sum_index_subscript(all_indexes(term), subspace_info; do_latex=do_latex)
    eq_counter, neq_constraints = eq_counter_and_neq_indexes(term.blocks, where_acting, subspace_info) 
    sup = ""
    # three cases eq_counter == 0 => all neq 
    if length(neq_constraints) == 0  # all eq 
        if eq_counter > 0
            sup = _format_superscript("=", do_latex=do_latex)
        end
    elseif eq_counter == 0
        sup = _format_superscript(do_latex ? "\\neq" : "≠", do_latex=do_latex)
    else
        neq_condition_strs = _format_neq_condition.(neq_constraints, Ref(subspace_info), do_latex=do_latex)
        sup = _format_superscript(join(neq_condition_strs, ","), do_latex=do_latex)
    end
    return base * sup
end
