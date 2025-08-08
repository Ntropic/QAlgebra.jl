export max_exponents, build_xpows, evaluate
"""
    max_exponents(a::CFunction)    -> Vector{Int}
    max_exponents(a::AbstractVector{CFunction}) -> Vector{Int}

Computes the elementwise maximum of variable exponents, 
ideal for getting the maximum requied orders for `build_xpows`.
"""
function max_exponents(vs::Vector{Vector{Int}})
    isempty(vs) && return Int[]
    return reduce((a, b) -> max.(a, b), vs)
end
max_exponents(a::CAtom) = abs.(a.var_exponents)
function max_exponents(s::CSum)
    return max_exponents([max_exponents(t) for t in s.terms])
end
function max_exponents(s::CProd)
    return max_exponents([max_exponents(t) for t in s.terms])
end
function max_exponents(r::CRational)
    max_a = max_exponents(r.numer)
    max_b = max_exponents(r.denom)
    return max_exponents([max_a, max_b])
end
max_exponents(e::CExp) = max_exponents(e.x)
max_exponents(l::CLog) = max_exponents(l.x)
max_exponents(v::AbstractVector{<:CFunction}) = max_exponents([max_exponents(t) for t in v])

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

ctimes(c::ComplexRational, d::T) where T <: Number = (c.a+im*c.b)/c.c * d

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
function evaluate(s::CSum, x::Vector{<:Number})
    sum(evaluate(t, x) for t in s.terms)
end
function evaluate(r::CRational, x::Vector{<:Number})
    evaluate(r.numer, x) / evaluate(r.denom, x)
end
function evaluate(e::CExp, x::Vector{<:Number})
    ctimes(e.coeff , exp(evaluate(e.x, x)))
end
function evaluate(l::CLog, x::Vector{<:Number})
    ctimes(l.coeff , log(evaluate(l.x, x)))
end
function evaluate(p::CProd, x::Vector{<:Number})
    # coefficient times the product of all term‐evaluations
    ctimes(p.coeff,
        prod(evaluate(term, x) for term in p.terms)
    )
end

# Evaluate but with xpows
function evaluate(a::CAtom, xpows::Vector{Vector{T}}) where T <: Number
    ctimes(a.coeff , prod(
        a.var_exponents[i] >= 0 ?
            xpows[i][a.var_exponents[i] + 1] :
            1 / xpows[i][abs(a.var_exponents[i]) + 1]
        for i in eachindex(a.var_exponents)))
end
function evaluate(s::CSum, xpows::Vector{Vector{T}}) where T <: Number
    return sum(evaluate(t, xpows) for t in s.terms)
end
function evaluate(r::CRational, xpows::Vector{Vector{T}}) where T <: Number
    return evaluate(r.numer, xpows) / evaluate(r.denom, xpows)
end
function evaluate(l::CExp, xpows::Vector{Vector{T}}) where T <: Number
    ctimes(l.coeff , exp(evaluate(l.x, xpows)))
end
function evaluate(l::CLog, xpows::Vector{Vector{T}}) where T <: Number
    ctimes(l.coeff , log(evaluate(l.x, xpows)))
end
function evaluate(p::CProd, xpows::Vector{Vector{T}}) where T<:Number
    # coefficient times the product of all term‐evaluations via xpows
    ctimes(p.coeff,
        prod(evaluate(term, xpows) for term in p.terms))
end