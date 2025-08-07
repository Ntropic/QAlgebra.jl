
import Base: +, -, *, /, ^, ==

# addition always builds a flat sum
+(a::T) where {T <: CFunction}  =  a 
+(a::Ta, b::Tb) where {Ta <: CFunction, Tb <: CFunction}  = CSum( vcat(_terms(a), _terms(b)) )
+(a::T, b::Number) where {T <: CFunction}  =  a+CAtom(b, zeros(Int, dims(a)))
+(b::Number, a::T) where {T <: CFunction}  =  a+b

# unary minus & subtraction
import Base: -, +

-(a::CAtom) = CAtom(-a.coeff, a.var_exponents)
-(s::CSum) = CSum([ -t for t in s.terms ])
-(r::CRational) = CRational(-r.numer, r.denom)
-(a::CProd) = CProd(-a.coeff, a.terms)
-(a::CExp) = CExp(-a.coeff, a.x)
-(a::CLog) = CLog(-a.coeff, a.x)
-(a::CFunction, b::CFunction) = a + (-b)
-(a::CFunction, b::Number)  = a + CAtom(-b, zeros(Int, dims(a)))
-(b::Number, a::CFunction)  = CAtom(b, zeros(Int, dims(a))) - a

# distribute * over sums
*(a::CSum, b::CSum) = CSum([ x*y for x in a.terms for y in b.terms ])
*(s::CSum, a::T) where {T<:CFunction}   = CSum([ x*a for x in s.terms ])
*(a::T, s::CSum) where {T<:CFunction}   = s*a

# atom‐level ×
*(a::CAtom, b::CAtom)     = CAtom(crationalize(a.coeff*b.coeff), a.var_exponents .+ b.var_exponents)
*(a::CAtom, r::CRational) = CRational(CSum(a)*r.numer, r.denom)
*(r::CRational, a::CAtom) = CRational(r.numer*CSum(a), r.denom)
*(a::CRational, b::CRational) = CRational(a.numer*b.numer, a.denom*b.denom)

*(a::CSum, b::CRational) = CRational(a*b.numer, b.denom)
*(b::CRational, a::CSum) = CRational(a*b.numer, b.denom)
function *(a::Ta, b::Tb) where {Ta <: CFunction, Tb <: CFunction}  
    ca = coeff(a)
    cb = coeff(b)
    if length(ca) != 1 || length(cb) != 1
        throw(ArgumentError("Cannot multiply two CFunctions with more than one coefficient (this function shouldn't be called by CSum)."))
    end
    return CProd(ca[1]*cb[1], sort!([a/ca[1], b/cb[1]]))
end
function *(a::CProd, b::T) where T <: CFunction 
    ca = coeff(a)
    cb = coeff(b)
    return CProd(ca[1]*cb[1], sort!(vcat(a.terms, b/cb[1])))
end
*(a::T, b::CProd) where T <: CFunction = b*a 
*(a::CSum, b::CProd) = CProd(b.coeff, vcat(a, b.terms))
*(b::CProd, a::CSum) = CProd(b.coeff, vcat(a, b.terms))
*(a::CProd, b::CProd) = CProd(a.coeff*b.coeff, sort!(vcat(a.terms, b.terms)))
*(a::CExp, b::CExp) = CExp(a.coeff*b.coeff, a.x+b.x)

# number 
function *(a::CAtom, b::Number)  
    return CAtom(a.coeff*b, a.var_exponents)
end
*(a::CSum, b::Number)  = CSum([ x*b for x in a.terms ])
*(a::CRational, b::Number) = CRational(a.numer*b, a.denom)
function *(a::CProd, b::Number)
    return CProd(a.coeff*b, a.terms)
end
function *(a::CLog, b::Number)
    return CLog(b*a.coeff, a.x)
end
function *(a::CExp, b::Number)
    return CExp(b*a.coeff, a.x)
end
*(b::Number, a::T) where T <: CFunction = a * b

multiply_one(a::CRational, b::Int) = (a.numer * b) / (a.denom * b)

# division
/(A::CSum, B::CSum)       = CRational(A, B)
/(A::CAtom, B::CSum)      = CRational(A, B)
/(A::CSum, b::CAtom)      = CSum([ x/b for x in A.terms ])
/(a::CAtom, b::CAtom)     = CAtom(a.coeff/b.coeff, a.var_exponents .- b.var_exponents)
/(a::CAtom, r::CRational) = CRational(CSum(a)*r.denom, r.numer)
/(r::CRational, a::CAtom) = CRational(r.numer, r.denom*CSum(a))
/(a::CRational, b::CRational) = CRational(a.numer*b.denom, a.denom*b.numer)
/(a::CRational, b::CSum)  = CRational(a.numer, a.numer*b)
/(a::CSum, b::CRational)  = CRational(a*b.denom, b.numer)
/(a::CAtom, b::Number)  = CAtom(a.coeff/b, a.var_exponents)
/(a::CRational, b::Number)  = CRational(a.numer, a.denom*b)
function /(a::CSum, n::Number) 
    return CSum([ x/n for x in a.terms ])
end
function /(a::CProd, n::Number) 
    return CProd(a.coeff/n, a.terms)
end
function /(a::CLog, n::Number) 
    return CLog(a.coeff/n, a.x)
end
function /(a::CExp, n::Number) 
    return CExp(a.coeff/n, a.x)
end
/(n::Number, a::T)  where T <: CFunction  = CAtom(n, zeros(Int, dims(a))) / a

/(a::T1, b::T2) where {T1 <: CFunction, T2 <: CFunction} = CRational(a, b)
function /(a::CProd, b::CAtom) 
    ind = findfirst(x -> isa(x, Union{CAtom, CSum, CRational}), a.terms)
    if ind === nothing
        return CRational(a, b) 
    end
    cb = b.coeff
    term_ind = a.terms[ind] / b * cb
    return CProd(a.coeff / cb, vcat(a.terms[1:ind-1], term_ind, a.terms[ind+1:end]))
end 


# exponentiation 
^(A::CAtom, n::Int) = CAtom(A.coeff^n, A.var_exponents .* n)
^(A::CSum, n::Int) = CSum([ x^n for x in A.terms ])
^(a::CRational, n::Int) = CRational(a.numer^n, a.denom^n)
^(a::CProd, n::Int) = CProd(a.coeff^n, [ x^n for x in A.terms ])
^(a::CExp, n::Int) = CExp(a.coeff^n, a.x * n)
^(a::CLog, n::Int) = error("Cannot take the power of a logarithm. Waiting to implement CPower type for this.")

import Base: inv 
inv(a::CAtom) = CAtom(a.coeff^(-1), a.var_exponents .*(-1))
inv(a::CSum) = CSum(inv.(a.terms))
inv(a::CRational) = CRational(inv(a.numer), inv(a.denom)) 

import Base: adjoint, conj 
function adjoint(f::CAtom)::CAtom
    return CAtom(conj(f.coeff), copy(f.var_exponents))
end
adjoint(f::CSum) = CSum([adjoint(t) for t in f.terms])
adjoint(f::CRational) = CRational(adjoint(f.numer), adjoint(f.denom))
conj(f::CFunction) = adjoint(f)

function ==(a::CAtom, b::CAtom)
    return (a.coeff == b.coeff && a.var_exponents == b.var_exponents)
end
function ==(a::CSum, b::CSum)
    if length(a) != length(b)
        return false
    end
    for i in 1:length(a) 
        if a[i] != b[i]
            return false
        end
    end
    return true
end
function ==(a::CRational, b::CRational)
    return (a.numer == b.numer && a.denom == b.denom)
end