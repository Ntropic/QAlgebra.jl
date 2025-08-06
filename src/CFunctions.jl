module CFunctions

using ..StringUtils
using ComplexRationals
using ..QAlgebra: get_default, FLIP_IF_FIRST_TERM_NEGATIVE, DO_BRACED

export CFunction, CAtom, CSum, CRational, CProd, CExp, CLog
export isnumeric

import Base: copy, exp, log, length
"""
    CFunction

Abstract supertype for symbolic functions representing atoms (`CAtom`), sums (`CSum`), and rationals (`CRational`).
"""
abstract type CFunction end

"""
    CAtom(coeff::Int, var_exponents::Vector{Int})
    CAtom(coeff::Rational, var_exponents::Vector{Int})
    CAtom(coeff::ComplexRational, var_exponents::Vector{Int})

A single term with a complex‐rational coefficient and integer exponents for each variable.
- The `Int` and `Rational` constructors wrap the coefficient into a `ComplexRational`.
- `var_exponents[j]` is the exponent of variable _j_.
"""
mutable struct CAtom <: CFunction
    coeff::ComplexRational
    var_exponents::Vector{Int}
    function CAtom(coeff::Int, var_exponents::Vector{Int})
        c = ComplexRational(coeff, 0, 1)
        return new(c, copy(var_exponents))
    end
    function CAtom(coeff::Rational, var_exponents::Vector{Int})
        c = ComplexRational(numerator(coeff), 0, denominator(coeff))
        return new(c, copy(var_exponents))
    end
    function CAtom(coeff::Complex, var_exponents::Vector{Int})
        c = crationalize(coeff)
        return new(c, copy(var_exponents))
    end
    function CAtom(coeff::ComplexRational, var_exponents::Vector{Int})
        return new(copy(coeff), copy(var_exponents))
    end
    function CAtom(coeff::Number, var_exponents::Vector{Int})
        c = crationalize(coeff+0im)
        return new(c, copy(var_exponents))
    end
end

function copy(x::CAtom)
    return CAtom(copy(x.coeff), copy(x.var_exponents))
end
coeff(a::CAtom) = [a.coeff]
var_exponents(a::CAtom) = [a.var_exponents]
dims(q::CAtom) = length(q.var_exponents)
length(a::CAtom) = 1


"""
    CSum(terms::AbstractVector{<:CFunction})
    CSum(ts::CFunction...)

Constructs a sum of `CFunction` terms.
- Flattens any nested `CSum` automatically.
- Variadic form `CSum(a, b, c)` is provided for convenience.
"""
mutable struct CSum <: CFunction
    terms::Vector{CFunction}
    function CSum(ts::AbstractVector{<:CFunction})
        flat = CFunction[]
        if length(ts) == 1
            return copy(ts[1])
        end
        for t in ts
            if t isa CSum
                append!(flat, copy.(t.terms))
            else
                push!(flat, copy(t))
            end
        end
        return new(flat)
    end
end

CSum(ts::CFunction...) = CSum(collect(ts))
copy(x::CSum) = CSum(copy.(x.terms))
coeff(x::CSum) = error("Sums don't have a coeff, you likely have a sum in a sum, this shouldn't happen. Please inform the developers. ")
var_exponents(x::CSum) = error("Sums don't have var_exponents, you likely have a sum in a sum, this shouldn't happen. Please inform the developers. ")
dims(q::CSum) = dims(q.terms[1])
length(q::CSum) = length(q.terms)


mutable struct CProd <: CFunction
    coeff::ComplexRational
    terms::Vector{CFunction}
    function CProd(coeff::ComplexRational, terms::AbstractVector{<:CFunction})
        terms_copy = copy.(terms)
        recursive_sort!(terms_copy)
        if length(terms_copy) == 1
            return terms_copy[1] / coeff
        end
        return new(copy(coeff), terms_copy)
    end
end
function CProd(terms::Vector{CFunction})
    CProd(ComplexRational(1, 0, 1), terms)
end
function CProd(ts::CFunction...)
    CProd(ComplexRational(1, 0, 1), collect(ts))
end
copy(x::CProd) = CProd(copy(x.coeff), copy.(x.terms))
coeff(x::CProd) = [x.coeff]
var_exponents(x::CProd) = vcat(var_exponents.(x.terms)...)
dims(q::CProd) = dims(q.terms[1])
length(q::CProd) = max(length.(q.terms)...)


"""
    CRational(numer::CSum, denom::CSum)

Represents a rational function with numerator `numer` and denominator `denom`, both sums of `CFunction` terms.
"""
mutable struct CRational <: CFunction
    numer::CFunction
    denom::CFunction
    function CRational(numer::CFunction, denom::CFunction)
        new(copy(numer), copy(denom))
    end
end
copy(x::CRational) = CRational(copy(x.numer), copy(x.denom))
coeff(x::CRational) = coeff.(x.numer) #/coeff(x.denom)
var_exponents(x::CRational) = var_exponents(x.numer)   #vcat(var_exponents.(x.numer), var_exponents.(var_exponents.(x.denom)))
dims(q::CRational) = dims(q.numer) 
length(q::CRational) = max(length(q.numer), length(q.denom))


mutable struct CExp <: CFunction
    coeff::ComplexRational
    x::CFunction
    function CExp(coeff::ComplexRational, x::CFunction)
        new(copy(coeff), copy(x))
    end
    function CExp(x::CFunction)
        new(ComplexRational(1,0,1), copy(x))
    end
end
function exp(x::CFunction)
    return CExp(x) 
end
copy(x::CExp) = CExp(copy(x.coeff), copy(x.x))
coeff(x::CExp) = [x.coeff]
var_exponents(x::CExp) = [zeros(Int, dims(x))]
dims(q::CExp) = dims(q.x)
length(q::CExp) = 1


mutable struct CLog <: CFunction
    coeff::ComplexRational
    x::CFunction
    function CLog(coeff::ComplexRational, x::CFunction)
        new(copy(coeff), copy(x))
    end
    function CLog(x::CFunction)
        new(ComplexRational(1,0,1), copy(x))
    end
end
function log(x::CFunction)
    return CLog(x) 
end
copy(x::CLog) = CLog(copy(x.coeff), copy(x.x))
coeff(x::CLog) = [x.coeff] 
var_exponents(x::CLog) = [zeros(Int, dims(x))]
dims(q::CLog) = dims(q.x)
length(q::CLog) = 1







#### Some basic functions ##############################################################################################


_terms(f::T) where T <: CFunction = [f]
_terms(f::CSum) =  f.terms 

import Base: iszero, isempty, isone
"""
    iszero(a::CAtom)     -> Bool
    iszero(s::CSum)      -> Bool
    iszero(r::CRational) -> Bool

Returns `true` if the expression is identically zero:
- **Atom**: zero coefficient.
- **Sum**: all terms zero or empty.
- **Rational**: zero numerator.
"""
iszero(a::CAtom)        = iszero(a.coeff)
iszero(s::CSum)         = isempty(s.terms) || all(iszero, s.terms)
iszero(p::CProd)        = iszero(p.coeff) || all(iszero, p.terms)
iszero(r::CRational)    = iszero(r.numer)
iszero(x::CExp)         = false
iszero(x::CLog)         = false

isempty(s::CSum)        = isempty(s.terms)

"""
    isnumeric(a::CAtom)    -> Bool
    isnumeric(s::CSum)     -> Bool
    isnumeric(r::CRational) -> Bool

Returns `true` if the expression contains no variables (i.e., all exponents are zero in atoms, and both numerator and denominator are numeric sums).
"""
isnumeric(a::CAtom)    = all(e->e==0, a.var_exponents)
isnumeric(s::CSum)     = all(isnumeric, s.terms)
isnumeric(p::CProd)    = all(isnumeric, p.terms)
isnumeric(r::CRational)= isnumeric(r.numer) && isnumeric(r.denom)
isnumeric(x::CExp)     = isnumeric(x.x)
isnumeric(x::CLog)     = isnumeric(x.x)

isone(c::CFunction) = false 
isone(a::CAtom)     = isnumeric(a) && isone(a.coeff)

allnegative(a::CAtom) = is_negative(a.coeff)
allnegative(s::CSum)  = !isempty(s.terms) && all(allnegative, s.terms)
allnegative(p::CProd) = is_negative(p.coeff)
allnegative(r::CRational) = allnegative(r.numer)
allnegative(x::CExp)  = false
allnegative(x::CLog)  = false


min_exponents(a::CAtom)::Vector{Int}    = a.var_exponents
function min_exponents(s::CSum)::Vector{Int}
    if length(s.terms) == 0
        return []
    end
    min_vals = min_exponents(s.terms[1])
    for term in s.terms[2:end]
        curr_min_vals = min_exponents(term)
        min_vals = min.(min_vals, curr_min_vals)
    end
    return min_vals
end
#function min_exponents(p::CProd)::Vector{Int}
#
#end
function min_exponents(r::CRational)::Vector{Int}
    min_vals = min_exponents(r.numer)
    min_vals2 = min_exponents(r.denom)
    return min.(min_vals, min_vals2)
end
function min_exponents(x::CExp)::Vector{Int}
    return zeros(Int, length(min_exponents(x.x)))
end
function min_exponents(x::CLog)::Vector{Int}
    return zeros(Int, length(min_exponents(x.x)))
end

import Base: length, getindex, iterate, deleteat!, reverse
length(p::CFunction)::Int = 1


getindex(p::CSum, i::Int) = p.terms[i]
iterate(p::CSum, state=1) = state > length(p.terms) ? nothing : (p.terms[state], state + 1)
deleteat!(p::CSum, i::Int) = CSum(deleteat!(p.terms, i))
reverse(q::CSum) = CSum(reverse(q.terms))


include("CFunctionsOps/CFunctions_algebra.jl")
include("CFunctionsOps/CFunctions_sort.jl")
include("CFunctionsOps/CFunctions_simplify.jl")
include("CFunctionsOps/CFunctions_orders_eval.jl")
include("CFunctionsOps/CFunctions_helper.jl")
include("CFunctionsOps/CFunctions_print.jl")

end # module CFunctions