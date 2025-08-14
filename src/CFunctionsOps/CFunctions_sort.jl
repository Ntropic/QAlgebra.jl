import Base: isless

# -------------------------------
# Small helpers (boolean style)
# -------------------------------

@inline eq_ordered(a, b) = !(isless(a, b) || isless(b, a))  # works for any type with isless

# For speed on primitives we can specialize equality directly:
@inline eq_ordered(a::Int,  b::Int)  = a == b
@inline eq_ordered(a::Bool, b::Bool) = a == b
@inline eq_ordered(a::Real, b::Real) = a == b   # assumes no NaNs in your sort fields

# Lexicographic compare for vectors of Ints (fast path)
function less_vec(a::AbstractVector{<:Integer}, b::AbstractVector{<:Integer})
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        ai = a[i]; bi = b[i]
        if ai != bi
            return ai < bi
        end
    end
    return length(a) < length(b)
end

# Generic lexicographic compare for vectors whose elements have isless
function less_vec(a::AbstractVector, b::AbstractVector)
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        ai = a[i]; bi = b[i]
        if !eq_ordered(ai, bi)
            return isless(ai, bi)
        end
    end
    return length(a) < length(b)
end

# Tuples (elements may be heterogenous but should have isless)
function less_tuple(a::Tuple, b::Tuple)
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        ai = a[i]; bi = b[i]
        if !eq_ordered(ai, bi)
            return isless(ai, bi)
        end
    end
    return length(a) < length(b)
end

# -------------------------------
# Tags for top-level CFunction kinds
# -------------------------------
cf_tag(::CAtom)     = 0
cf_tag(::CSum)      = 1
cf_tag(::CProd)     = 2
cf_tag(::CRational) = 3
cf_tag(::CExp)      = 4
cf_tag(::CLog)      = 5

# -------------------------------
# Same-kind comparisons (boolean)
# -------------------------------

# CAtom: compare var_exponents lexicographically
isless_same(a::CAtom, b::CAtom) = less_vec(a.var_exponents, b.var_exponents)

# CSum: length first, then terms pairwise
function isless_same(a::CSum, b::CSum)
    la = length(a.terms); lb = length(b.terms)
    la != lb && return la < lb
    return less_vec(a.terms, b.terms)   # elements are CFunction → uses isless below
end

# CProd: length first, then terms pairwise
function isless_same(a::CProd, b::CProd)
    la = length(a.terms); lb = length(b.terms)
    la != lb && return la < lb
    return less_vec(a.terms, b.terms)
end

# CRational: denominator first (priority), then numerator
function isless_same(a::CRational, b::CRational)
    if !eq_ordered(a.denom, b.denom)
        return isless(a.denom, b.denom)
    end
    return isless(a.numer, b.numer)
end

# CExp / CLog: just compare inner expr
isless_same(a::CExp, b::CExp) = isless(a.x, b.x)
isless_same(a::CLog, b::CLog) = isless(a.x, b.x)

# -------------------------------
# Cross-type dispatcher (boolean)
# -------------------------------
function isless(a::CFunction, b::CFunction)
    ta = cf_tag(a); tb = cf_tag(b)
    ta != tb && return ta < tb
    return isless_same(a, b)
end
