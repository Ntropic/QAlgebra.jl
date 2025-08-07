export expand
# function expand(r::CRational, order

"""
    expand(f::CExp, order::Int) -> CSum

Taylor-expand `coeff*exp(x)` up to `order`:
  
  coeff * ∑_{n=0}^order x^n / n!
"""
function expand(f::CExp, order::Int)
    # Pre-allocate the output vector
    terms = Vector{CFunction}(undef, order+1)
    # zero-exponent atom for x^0
    zero_exp = zeros(Int, dims(f.x))
    const1 = CAtom(ComplexRational(1,0,1), zero_exp)

    @inbounds for n in 0:order
        # x^n
        base = n == 0 ? const1 : f.x^n
        # coefficient = f.coeff * (1/n!)
        # 1/n! = ComplexRational(1,0,factorial(n))
        coef = f.coeff * ComplexRational(1, 0, factorial(n))
        # multiply by scalar
        terms[n+1] = base * coef
    end

    return CSum(terms)
end

"""
    expand(r::CRational) -> CFunction

Distribute a fraction over a sum in the numerator:
  
  (A + B + …) / D  ⇒  A/D  +  B/D  +  …
"""
function expand(r::CRational)
    if r.numer isa CSum
        # pre-allocate
        nt = length(r.numer.terms)
        out = Vector{CFunction}(undef, nt)
        @inbounds for i in 1:nt
            out[i] = r.numer.terms[i] / r.denom
        end
        return CSum(out)
    else
        # nothing to expand
        return r
    end
end

# an override for generic dispatch: recurse into sums/products,
# but push exponential/fraction expansions when encountered
function expand(f::CFunction, order::Int)
    if f isa CExp
        return expand(f, order)
    elseif f isa CRational
        return expand(f)
    elseif f isa CSum
        # pre-allocate
        ts = f.terms
        nt = length(ts)
        out = Vector{CFunction}(undef, nt)
        @inbounds for i in 1:nt
            out[i] = expand(ts[i], order)
        end
        return CSum(out)
    elseif f isa CProd
        # leave coefficient intact, expand each factor
        tms = f.terms
        nt = length(tms)
        out = Vector{CFunction}(undef, nt)
        @inbounds for i in 1:nt
            out[i] = expand(tms[i], order)
        end
        return CProd(f.coeff, out, Val{:nosimp}())
    else
        # CAtom, CLog, etc. — just return as-is (or recurse into log)
        return f isa CLog ? CLog(f.coeff, expand(f.x, order)) : f
    end
end