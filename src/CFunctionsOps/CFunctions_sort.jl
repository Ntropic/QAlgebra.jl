export sort_key, recursive_sort

import Base: sort!, sort


"""
    sort_key(f::CFunction) -> Vector{<:Real}

Returns a vector used to sort symbolic expressions in a canonical order:
- `CAtom`: variable exponents followed by a coefficient key.
- `CSum`: length and keys of terms.
- `CRational`: denominator keys first (priority), then numerator.
"""
function sort_key(a::CAtom)
    return a.var_exponents #vcat(a.var_exponents, customsort_key(a.coeff))
end

function sort_key(s::CSum)
    keys = [sort_key(t) for t in s.terms]
    return vcat(1, length(keys), keys...)
end
function sort_key(p::CProd)
    keys = [sort_key(t) for t in p.terms]
    return vcat(2, length(keys), keys...)
end
function sort_key(r::CRational)
    dkeys = sort_key(r.denom)
    nkeys = sort_key(r.numer)
    return vcat(3, dkeys, nkeys)
end
function sort_key(e::CExp)
    return vcat(4, sort_key(e.x))
end
function sort_key(l::CLog)
    return vcat(5, sort_key(l.x))
end

import Base: isless
function isless(a::T1, b::T2) where {T1 <: CFunction, T2 <: CFunction}
    return sort_key(a) < sort_key(b)
end


function _sort!(x::AbstractVector{<:CFunction}; rev::Bool=false, alg=InsertionSort)
    Base.sort!(x; by=sort_key, rev=rev, alg=alg)
    return x
end
function _sort(x::AbstractVector{<:CFunction};  rev::Bool=false, alg=InsertionSort)
    return sort(x, by=sort_key, rev=rev, alg=alg)
end