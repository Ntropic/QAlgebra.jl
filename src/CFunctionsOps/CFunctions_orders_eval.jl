export max_exponents, build_xpows, evaluate
"""
    max_exponents(a::CAtom)    -> Vector{Int}
    max_exponents(s::CSum)      -> Vector{Int}
    max_exponents(r::CRational) -> Vector{Int}

Computes the elementwise maximum of variable exponents:
- **Atom**: absolute value of its own exponents.
- **Sum**: maximum across all terms.
- **Rational**: maximum between numerator and denominator exponents.
"""
max_exponents(a::CAtom) = abs.(a.var_exponents)
function max_exponents(s::CSum)
    if length(s.terms) == 0
        return []
    else
        max_m = max_exponents(s[1])
        for t in s.terms[2:end]
            curr_max = max_exponents(t)
            for i in eachindex(curr_max)
                if curr_max[i] > max_m[i]
                    max_m[i] = curr_max[i] 
                end
            end
        end
        return max_m
    end
end
function max_exponents(r::CRational)
    max_a = max_exponents(r.numer)
    max_b = max_exponents(r.denom)
    return [max(a, b) for (a, b) in zip(max_a, max_b)]
end

"""
    build_xpows(x::Vector{<:Number}, max_exp::Vector{Int}) -> Vector{Vector}

Precomputes powers of each variable for fast evaluation:
- Returns `xpows` such that `xpows[j][k+1] == x[j]^k` for `k = 0:max_exp[j]`.
- Lengths of `x` and `max_exp` must agree.
"""
function build_xpows(x::Vector{<:Number}, max_exp::Vector{Int})::Vector{Vector}
    @assert length(x) == length(M)
    xpows = [Vector{typeof(x[i])}(undef, length(x)) for i in eachindex(x)]
    for j in eachindex(x)
        xpows[j] = [ x[j]^k for k in 0:max_exp[j] ]
    end
    return xpows
end

"""
    evaluate(a::CAtom, xpows::Vector{Vector})       -> Number
    evaluate(s::CSum, xpows::Vector{Vector})        -> Number
    evaluate(r::CRational, xpows::Vector{Vector})   -> Number

Evaluates the expression using precomputed powers `xpows` (as returned by `build_xpows`) for faster repeated evaluation.
"""
function evaluate(a::CAtom, x::AbstractVector{<:Number})
    a.coeff * prod(x[i]^a.var_exponents[i] for i in eachindex(a.var_exponents))
end
function evaluate(s::CSum, x::AbstractVector{<:Number})
    sum(evaluate(t, x) for t in s.terms)
end
function evaluate(r::CRational, x::AbstractVector{<:Number})
    evaluate(r.numer, x) / evaluate(r.denom, x)
end

# Evaluate but with xpows
function evaluate(a::CAtom, xpows::Vector{Vector})
    a.coeff * prod(
        a.var_exponents[i] >= 0 ?
            xpows[i][a.var_exponents[i] + 1] :
            1 / xpows[i][abs(a.var_exponents[i]) + 1]
        for i in eachindex(a.var_exponents)
    )
end
function evaluate(s::CSum, xpows::Vector{Vector})
    return sum(evaluate(t, xpows) for t in s.terms)
end
function evaluate(r::CRational, xpows::Vector{Vector})
    return evaluate(r.numer, xpows) / evaluate(r.denom, xpows)
end