export sort_key, recursive_sort

import Base: sort!, sort


"""
    sort_key(f::CFunction) -> Vector{<:Real}

Returns a vector used to recursive_sort symbolic expressions in a canonical order:
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


function recursive_sort!(s::CAtom; kwargs...)
    return s
end
function recursive_sort!(s::CSum; kwargs...)
    sort!(s.terms; kwargs...)
    return s
end
function recursive_sort!(s::CProd; kwargs...) 
    recursive_sort!(s.terms; kwargs...) 
    return s
end
function recursive_sort!(s::CRational; kwargs...)
    recursive_sort!(s.numer; kwargs...)
    recursive_sort!(s.denom; kwargs...)
    return s
end
function recursive_sort!(s::CExp; kwargs...) 
    recursive_sort!(s.x; kwargs...)
    return s
end
function recursive_sort!(s::CLog; kwargs...)
    recursive_sort!(s.x; kwargs...)
    return s
end
function recursive_sort!(x::AbstractVector{<:CFunction}; kwargs...)
    for i in eachindex(x)
        recursive_sort!(x[i]; kwargs...)
    end
    invoke(Base.sort!, Tuple{AbstractVector}, x; by=sort_key, kwargs...)
    return x
end

function recursive_sort(s::CAtom; by=sort_key)::CAtom
    # do nothing 
    return s
end
function recursive_sort(s::CSum; kwargs...)::CSum
    return CSum(sort(s.terms; kwargs...))
end
function recursive_sort(s::CProd; kwargs...)::CProd
    return CProd(recursive_sort(s.terms; kwargs...))
end
function recursive_sort(s::CRational; kwargs...)::CRational
    numer = recursive_sort(s.numer; kwargs...)
    denom = recursive_sort(s.denom; kwargs...)
    return CRational(numer, denom)
end
function recursive_sort(s::CExp; kwargs...)::CExp
    return CExp(recursive_sort(s.x; kwargs...))
end
function recursive_sort(s::CLog; kwargs...)::CLog
    return CLog(recursive_sort(s.x; kwargs...))
end
function recursive_sort(x::AbstractVector{<:CFunction}; kwargs...)::AbstractVector{<:CFunction}
    return invoke(Base.sort, Tuple{AbstractVector}, x; by=sort_key, kwargs...)
end

function sort!(x::AbstractVector{<:CFunction}; rev::Bool=false, alg=InsertionSort)
     invoke(Base.sort!, Tuple{AbstractVector}, x; by=sort_key, rev=rev, alg=alg)
    return x
end
function sort(x::AbstractVector{<:CFunction};  rev::Bool=false, alg=InsertionSort)
     invoke(Base.sort, Tuple{AbstractVector}, x; by=sort_key, rev=rev, alg=alg)
    return x
end