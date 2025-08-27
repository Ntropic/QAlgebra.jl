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
function QExprLookup(ops_comb::Vector{Vector{Symbol}}, ops_vec::Vector{QExpr}; use_dict::Bool=false)
    @assert length(ops_comb) == length(ops_vec)
    if use_dict
        lookup = Dict{Tuple{Vararg{Symbol}}, QExpr}()
        for (c, op) in zip(ops_comb, ops_vec)
            lookup[Tuple(c)...] = op
        end
        return QExprLookup(ops_comb, ops_vec, true, lookup)
    else
        return QExprLookup(ops_comb, ops_vec, false, nothing)
    end
end

# Make it callable
function (ql::QExprLookup)(vars...)
    syms = Tuple(Symbol.(vars))  # normalize to Tuple
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
function Base.getindex(ql::QExprLookup, vars...)
    return ql(vars...)   # just forward to the call
end
# Provide available keys
import Base: keys
keys(ql::QExprLookup) = ql.ops_comb

