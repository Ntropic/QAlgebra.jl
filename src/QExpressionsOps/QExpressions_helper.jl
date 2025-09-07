export QExprLookup 
"""
    flatmap_to(f, xs, ::Type{T}) where T

Apply `f` to each element of `xs` (where each result is a `Vector{T}`),
and return all results concatenated into a single `Vector{T}`.
"""
function flatmap_to(f, xs, ::Type{T})::Vector{T} where T
    out = Vector{T}()
    for x in xs
        append!(out, f(x))
    end
    return out
end

""" 
    QExprLookup 

Keeps Vectors of QExprs accessible via their Symbols either as a function x(:x,:i) or vector x[:x,:i]. But similarly via string indexing.
""" 
struct QExprLookup
    ops_comb::Vector{Vector{Symbol}}
    ops_vec::Vector{QExpr}
    use_dict::Bool
    dict::Union{Dict{Tuple{Vararg{Symbol}}, QExpr}, Nothing}
end
function QExprLookup(ops_comb::Vector{Vector{Symbol}}, ops_vec::Vector{QExpr};
                     use_dict::Bool=false)
    @assert length(ops_comb) == length(ops_vec)

    # if there’s only one operator, just return it directly
    if length(ops_vec) == 1
        return ops_vec[1]
    end

    out_keys = Vector{Vector{Symbol}}()
    out_vals = QExpr[]

    for (k, op) in zip(ops_comb, ops_vec)
        push!(out_keys, k)
        push!(out_vals, op)

        # if key contains :t0, add alias with that element removed
        if :t0 in k
            k_alias = filter(!=(:t0), k)
            push!(out_keys, k_alias)
            push!(out_vals, op)
        end
    end

    if use_dict
        lookup = Dict{Tuple{Vararg{Symbol}}, QExpr}()
        for (c, op) in zip(out_keys, out_vals)
            lookup[Tuple(c)...] = op
        end
        return QExprLookup(out_keys, out_vals, true, lookup)
    else
        return QExprLookup(out_keys, out_vals, false, nothing)
    end
end

# Make it callable
function (ql::QExprLookup)(params...)
    syms = Tuple(Symbol.(params))  # normalize to Tuple
    if ql.use_dict
        return get(ql.dict, syms) do
            error("No operator found for input $(syms)")
        end
    else
        for (c, op) in zip(ql.ops_comb, ql.ops_vec)
            if Tuple(c) == syms
                return op
            end
        end
        error("No operator found for input $(syms)")
    end
end
function Base.getindex(ql::QExprLookup, params...)
    return ql(params...)   # just forward to the call
end
# Provide available keys
import Base: keys
keys(ql::QExprLookup) = ql.ops_comb

