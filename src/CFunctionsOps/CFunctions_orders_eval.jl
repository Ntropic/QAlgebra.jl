export max_exponents, build_xpows, evaluate
"""
    max_exponents(a::CFunction)    -> Vector{Int}
    max_exponents(a::AbstractVector{CFunction}) -> Vector{Int}

Computes the elementwise maximum of variable exponents, 
ideal for getting the maximum requied orders for `build_xpows`.
"""
function max_vec(a::Vector{Int}, b::Vector{Int})::Vector{Int}
    return max.(a, b)
end
max_exponents(a::CAtom) = Vector(abs.(a.var_exponents))
max_exponents(a::CAbstract) = error("max_exponents can not be used for CAbstract.")
function max_exponents(s::CFunction)::Vector{Int}
    m = zeros(Int, dims(s))
    for leaf in leaf_iter
        m = max_vec(m, max_exponents(leaf))
    end
    return m
end

"""
    build_xpows(x::Vector{<:Number}, max_exp::Vector{Int}) -> Vector{Vector}

Precomputes powers of each variable for fast evaluation:
- Returns `xpows` such that `xpows[j][k+1] == x[j]^k` for `k = 0:max_exp[j]`.
- Lengths of `x` and `max_exp` must agree.
"""
function build_xpows(x::Vector{T}, max_exp::Vector{Int})::Vector{Vector{T}}  where T <: Number
    @assert length(x) == length(max_exp)
    xpows = [Vector{typeof(x[i])}(undef, length(x)) for i in eachindex(x)]
    for j in eachindex(x)
        xpows[j] = [ x[j]^k for k in 0:max_exp[j] ]
    end
    return xpows
end

# Rational exponent on a numeric base
@inline function _pow_r(y, q::Rational{Int})
    if denominator(q) == 1
        return y ^ Int(q)
    else
        return y ^ float(q)    # works for real/complex bases
    end
end
ctimes(c::ComplexRational, d::T) where T <: Number = (c.a+im*c.b)/c.c * d

const VType{T<:Number} = Union{Vector{T}, Vector{Vector{T}}}
"""
    evaluate(p::CFunction, x::Vector{<:Number})        -> Number
    evaluate(p::CFunction, xpows::Vector{Vector{T}})   -> Number where T<:Number

Evaluate a `CFunction` by multiplying its coefficient with the evaluations of each factor.
- In the first form, each factor is evaluated directly at the point `x`.
- In the second form, each factor is evaluated using precomputed powers `xpows`, see the `build_xpows` function.
"""
function evaluate(a::CAtom, x::Vector{<:Number})
    ctimes(a.coeff , prod(x[i]^a.var_exponents[i] for i in eachindex(a.var_exponents)))
end

function evaluate(a::CAtom, xpows::Vector{Vector{T}}) where T <: Number
    ctimes(a.coeff , prod(
        a.var_exponents[i] >= 0 ?
            xpows[i][a.var_exponents[i] + 1] :
            1 / xpows[i][abs(a.var_exponents[i]) + 1]
        for i in eachindex(a.var_exponents)))
end

# Navigate the tree 
function evaluate(a::CAbstract, x::VType{T}) where T <: Number
    error("CAbstract requires evaluated values for its substitution to be itself evaluated.")
end
function evaluate(s::CSum, x::VType{T}) where T <: Number
    sum(evaluate(t, x) for t in s.expr)
end
function evaluate(r::CRational, x::VType{T}) where T <: Number
    evaluate(r.numer, expr) / evaluate(r.denom, x)
end
function evaluate(e::CExp, x::VType{T}) where T <: Number
    ctimes(e.coeff , exp(evaluate(e.expr, x)))
end
function evaluate(l::CLog, x::VType{T}) where T <: Number
    ctimes(l.coeff , log(evaluate(l.expr, x)))
end
function evaluate(p::CProd, x::VType{T}) where T <: Number
    # coefficient times the product of all term‐evaluations
    ctimes(p.coeff, prod(evaluate(term, x) for term in p.terms))
end
function evaluate(p::CPower, x::VType{T}) where T <: Number
    ctimes(p.coeff, _pow_r(evaluate(p.expr, x), p.exponent))
end
function evaluate(v::CVector, x::VType{T}) where T <: Number
    [ ctimes(v.coeff, evaluate(e, x)) for e in v.expr ]
end
function evaluate(M::CMatrix, x::VType{T}) where T <: Number
    m, n = size(M.expr)
    reshape([ ctimes(M.coeff, evaluate(e, x)) for e in M.expr ], m, n)
end
function evaluate(c::CCustomType, x::VType{T}) where T <: Number
    arg_values = [evaluate(e, x) for e in c.expr]
    val = evaluate(c.ctype_def.fun, arg_values, c.ctype_def.index_map)
    return ctimes(c.coeff, val)
end

# With substitutions 
function evaluate(a::CAtom, x::Vector{<:Number}, abstract_vals::Vector{<:Number}, index_map::Vector{Int}) 
    ctimes(a.coeff , prod(x[i]^a.var_exponents[i] for i in eachindex(a.var_exponents)))
end
function evaluate(a::CAtom, xpows::Vector{Vector{T}}, abstract_vals::Vector{<:Number}, index_map::Vector{Int}) where T <: Number
    ctimes(a.coeff , prod(
        a.var_exponents[i] >= 0 ?
            xpows[i][a.var_exponents[i] + 1] :
            1 / xpows[i][abs(a.var_exponents[i]) + 1]
        for i in eachindex(a.var_exponents)))
end
function evaluate(a::CAbstract, x::VType{T}, abstract_vals::Vector{<: Number}, index_map::Vector{Int})  where T <: Number
    new_val = abstract_vals[index_map[a.index]]
    if a.dag
        new_val = new_val'
    end
    if a.exponent != 1 
        new_val = new_val^a.exponent
    end
    return ctimes(a.coeff, new_val)
end
function evaluate(s::CSum, x::VType{T}, abstract_vals::Vector{<: Number}, index_map::Vector{Int})  where T <: Number
    sum(evaluate(t, x, abstract_vals, index_map) for t in s.expr)
end
function evaluate(r::CRational, x::VType{T}, abstract_vals::Vector{<: Number}, index_map::Vector{Int})  where T <: Number
    evaluate(r.numer, x, abstract_vals, index_map) / evaluate(r.denom, x, abstract_vals, index_map)
end
function evaluate(e::CExp, x::VType{T}, abstract_vals::Vector{<: Number}, index_map::Vector{Int})  where T <: Number
    ctimes(e.coeff , exp(evaluate(e.expr, x, abstract_vals, index_map)))
end
function evaluate(l::CLog, x::VType{T}, abstract_vals::Vector{<: Number}, index_map::Vector{Int})  where T <: Number
    ctimes(l.coeff , log(evaluate(l.expr, x, abstract_vals, index_map)))
end
function evaluate(p::CProd, x::VType{T}, abstract_vals::Vector{<: Number}, index_map::Vector{Int})  where T <: Number
    # coefficient times the product of all term‐evaluations
    ctimes(p.coeff, prod(evaluate(term, x, abstract_vals, index_map) for term in p.expr))
end
function evaluate(p::CPower, x::VType{T}, abstract_vals::Vector{<: Number}, index_map::Vector{Int})  where T <: Number
    ctimes(p.coeff, _pow_r(evaluate(p.expr, x, abstract_vals, index_map), p.exponent))
end
function evaluate(v::CVector, x::VType{T}, abstract_vals::Vector{<: Number}, index_map::Vector{Int})  where T <: Number
    [ ctimes(v.coeff, evaluate(e, x, abstract_vals, index_map)) for e in v.expr ]
end
function evaluate(M::CMatrix, x::VType{T}, abstract_vals::Vector{<: Number}, index_map::Vector{Int})  where T <: Number
    m, n = size(M.expr)
    reshape([ ctimes(M.coeff, evaluate(e, x, abstract_vals, index_map)) for e in M.expr ], m, n)
end
function evaluate(c::CCustomType, x::VType{T}, abstract_vals::Vector{<: Number}, index_map::Vector{Int})  where T <: Number
    arg_values = [evaluate(e, x, abstract_vals, index_map) for e in c.expr]
    val = evaluate(c.ctype_def.fun, arg_values, c.ctype_def.index_map)
    return ctimes(c.coeff, val)
end
