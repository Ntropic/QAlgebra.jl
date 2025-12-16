@inline aggregator_symbol(::Type{SumAggregator}, n::Int; do_latex::Bool=false) = do_latex ? "\\sum" : "∑"
@inline function aggregator_symbol(::Type{IntegralAggregator}, n::Int; do_latex::Bool=false)
    count = max(n, 1)
    if do_latex
        count == 1 && return "\\int"
        count == 2 && return "\\iint"
        count == 3 && return "\\iiint"
        return "\\idotsint"
    else
        count == 1 && return "∫"
        count == 2 && return "∬"
        count == 3 && return "∭"
        return repeat("∫", count)
    end
end
@inline function aggregator_symbol(::Type{A}, n::Int; do_latex::Bool=false) where {A<:AbstractQAggregator}
    error("No aggregator symbol defined for " * string(nameof(A)))
end

@inline aggregator_symbol(term::AbstractQSum, indices::Vector{SubSpaceIndex}; do_latex::Bool=false) =
    aggregator_symbol(aggregator_type(term), length(indices); do_latex=do_latex)

function _sum_index_subscript(term::AbstractQSum, indices::Vector{SubSpaceIndex}, info::SubSpaceInfo; do_latex::Bool=false)
    base = aggregator_symbol(term, indices; do_latex=do_latex)
    isempty(indices) && return base
    labels = Index2String.(indices, Ref(info))
    is_integral = aggregator_type(term) === IntegralAggregator
    if do_latex
        body = join(labels, ",")
        if is_integral
            density = raw"\rho_{" * body * "}"
            return base * "_{" * density * "}"
        else
            return base * "_{" * body * "}"
        end
    else
        seq = join(labels, ",")
        formatted = length(labels) == 1 ? seq : "(" * seq * ")"
        if is_integral
            density = "ρ" * formatted
            return base * str2sub(density)
        else
            sub = str2sub(formatted)
            return base * sub
        end
    end
end

function _integral_measure(term::AbstractQSum, indices::Vector{SubSpaceIndex}, info::SubSpaceInfo; do_latex::Bool=false)
    aggregator_type(term) === IntegralAggregator || return ""
    isempty(indices) && return ""
    pieces = String[]
    for idx in indices
        label = Index2String(idx, info)
        if do_latex
            push!(pieces, raw"\,\mathrm{d}" * label)
        else
            push!(pieces, " d" * label)
        end
    end
    return join(pieces, "")
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



function sum_symbol_str(term::AbstractQSum, where_acting::Vector{BitVector}; do_latex::Bool=false)
    subspace_info = term.qspace.subspace_info
    indices = all_indices(term)
    base = _sum_index_subscript(term, indices, subspace_info; do_latex=do_latex)
    eq_counter, neq_constraints = eq_counter_and_neq_indices(term.blocks, where_acting, subspace_info) 
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
    return base * sup, _integral_measure(term, indices, subspace_info; do_latex=do_latex)
end
