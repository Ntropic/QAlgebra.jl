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
@inline function _as_lookup_symbol(param)::Symbol
    param isa Symbol && return param
    param isa AbstractString && return Symbol(param)
    error("Unsupported lookup key $(param)::$(typeof(param)); expected Symbol or String.")
end

@inline function _canonical_lookup_symbol(sym::Symbol)::Symbol
    sym === :t && return :t0
    return sym
end

@inline function _canonical_lookup_symbols(params)::Tuple{Vararg{Symbol}}
    return Tuple(_canonical_lookup_symbol(_as_lookup_symbol(p)) for p in params)
end

@inline function _format_lookup_combo(combo)::String
    canonical = Tuple(_canonical_lookup_symbol(sym) for sym in combo)
    return "(" * join(string.(canonical), ", ") * ")"
end


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
            lookup[Tuple(_canonical_lookup_symbol(sym) for sym in c)] = op
        end
        return QExprLookup(out_keys, out_vals, true, lookup)
    else
        return QExprLookup(out_keys, out_vals, false, nothing)
    end
end
@inline function _available_lookup_combos(ql::QExprLookup)::Vector{String}
    combos = [_format_lookup_combo(c) for c in ql.ops_comb]
    return unique(combos)
end


# Make it callable
function (ql::QExprLookup)(params...)
    syms = _canonical_lookup_symbols(params)
    combos = _available_lookup_combos(ql)
    if ql.use_dict
        return get(ql.dict, syms) do
            error("No entry for $(syms). Pick one of $(join(combos, ", ")).")
        end
    else
        for (c, op) in zip(ql.ops_comb, ql.ops_vec)
            if Tuple(_canonical_lookup_symbol(sym) for sym in c) == syms
                return op
            end
        end
        error("No entry for $(syms). Pick one of $(join(combos, ", ")).")
    end
end
function Base.getindex(ql::QExprLookup, params...)
    return ql(params...)   # just forward to the call
end
# Provide available keys
import Base: keys
keys(ql::QExprLookup) = ql.ops_comb
