############ CFunctions_algebra.jl (updated for new structs) ################
import Base: +, -, *, /, ^, sqrt, ==, inv, adjoint, conj, transpose
using ..CFunctions: CFunction, CAtomic, CComposite, CMultiComposite,
                    CAbstract, CCustomType,
                    CAtom, CSum, CRational, CProd, CExp, CLog, CPower,
                    CVector, CMatrix, _CSum, coeff, length, modify_expr, modify_exprs
using ComplexRationals: ComplexRational, crationalize
using ..CFunctions: CR_ZERO, CR_ONE

# -------- helpers -----------------------------------------------------------
@inline pinfo(f::CFunction) = f.param_info
@inline dims(f::CFunction)  = f.param_info.dims

@inline function _ensure_same_param_info(a::CFunction, b::CFunction)
    if pinfo(a) !== pinfo(b)
        error("ParameterInfo mismatch between operands.")
    end
end

@inline function num_atom(f::CFunction, n::Number)
    CAtom(pinfo(f), crationalize(n + 0im), zeros(Int, dims(f)))
end
@inline zero_atom(f::CFunction) = CAtom(pinfo(f), CR_ZERO,  zeros(Int, dims(f)))
@inline one_atom(f::CFunction)  = CAtom(pinfo(f), CR_ONE,   zeros(Int, dims(f)))

# Clean accessors for “terms”
_terms(f::CFunction) = [f]
_terms(s::CSum)      = s.expr

# convenience
const _ScalarLike = Union{CAtom, CAbstract, CSum, CRational, CCustomType}

# -------- addition ----------------------------------------------------------
+(a::T) where {T<:CFunction} = a

function +(a::Ta, b::Tb) where {Ta<:CFunction, Tb<:CFunction}
    _ensure_same_param_info(a, b)
    _CSum(pinfo(a), vcat(_terms(a), _terms(b)))
end

+(a::T, b::Number) where {T<:CFunction} = a + num_atom(a, b)
+(b::Number, a::T) where {T<:CFunction} = a + b

function +(u::CVector, v::CVector)
    u.row == v.row || error("Vector orientations must match for addition.")
    length(u.expr) == length(v.expr) || error("Vector lengths must match.")
    _ensure_same_param_info(u.expr[1], v.expr[1])
    CVector(pinfo(u), [ u.coeff*ui + v.coeff*vi for (ui,vi) in zip(u.expr, v.expr) ]; row=u.row)
end

function +(A::CMatrix, B::CMatrix)
    size(A.expr) == size(B.expr) || error("Matrix sizes must match for addition.")
    m, n = size(A.expr)
    _ensure_same_param_info(A.expr[1,1], B.expr[1,1])
    CMatrix(pinfo(A), [ A.coeff*A.expr[i,j] + B.coeff*B.expr[i,j] for i in 1:m, j in 1:n ])
end

# -------- unary minus & subtraction ----------------------------------------
-(a::CAtom)       = CAtom(pinfo(a), -a.coeff, a.var_exponents)
-(a::CAbstract)   = CAbstract(a.param_info, -a.coeff, a.index, a.exponent, a.dag)
-(f::CCustomType) = CCustomType(f.param_info, -f.coeff, f.expr, f.ctype_def)

-(s::CSum)        = _CSum(pinfo(s), [ -t for t in s.expr ], Val(:nosimp))
-(r::CRational)   = CRational(pinfo(r), -r.numer, r.denom, Val(:nosimp))
-(a::CProd)       = CProd(pinfo(a), -a.coeff, a.expr, Val(:nosimp))
-(a::CExp)        = CExp(pinfo(a), -a.coeff, a.expr, Val(:nosimp))
-(a::CLog)        = CLog(pinfo(a), -a.coeff, a.expr, Val(:nosimp))
-(p::CPower)      = CPower(p.param_info, -p.coeff, p.expr, p.exponent, Val(:nosimp))
-(v::CVector)     = CVector(pinfo(v), -v.coeff, v.expr; row=v.row)
-(M::CMatrix)     = CMatrix(pinfo(M), -M.coeff, M.expr)

-(a::CFunction, b::CFunction) = a + (-b)
-(a::CFunction, b::Number)    = a + num_atom(a, -b)
-(b::Number, a::CFunction)    = num_atom(a, b) - a

function -(u::CVector, v::CVector)
    u.row == v.row || error("Vector orientations must match for subtraction.")
    length(u.expr) == length(v.expr) || error("Vector lengths must match.")
    _ensure_same_param_info(u.expr[1], v.expr[1])
    CVector(pinfo(u), [ u.coeff*ui - v.coeff*vi for (ui,vi) in zip(u.expr, v.expr) ]; row=u.row)
end

function -(A::CMatrix, B::CMatrix)
    size(A.expr) == size(B.expr) || error("Matrix sizes must match for subtraction.")
    m, n = size(A.expr)
    _ensure_same_param_info(A.expr[1,1], B.expr[1,1])
    CMatrix(pinfo(A), [ A.coeff*A.expr[i,j] - B.coeff*B.expr[i,j] for i in 1:m, j in 1:n ])
end

# -------- multiplication (distribute over sums) -----------------------------
*(a::CSum, b::CSum)      = _CSum(pinfo(a), [ x*y for x in a.expr for y in b.expr ])
*(s::CSum, a::CFunction) = _CSum(pinfo(s), [ x*a for x in s.expr ])
*(a::CFunction, s::CSum) = s*a

# atom-level ×
function *(a::CAtom, b::CAtom)
    _ensure_same_param_info(a, b)
    CAtom(pinfo(a), crationalize(a.coeff*b.coeff), a.var_exponents .+ b.var_exponents)
end

# abstract × abstract -> merge if same index & dag
function *(a::CAbstract, b::CAbstract)
    _ensure_same_param_info(a, b)
    if a.index == b.index && a.dag == b.dag
        return CAbstract(a.param_info, a.coeff*b.coeff, a.index, a.exponent + b.exponent, a.dag)
    else
        # generic fallback: becomes a product factor
        return CProd(pinfo(a), a.coeff*b.coeff, sort!([modify_coeff(a, CR_ONE), modify_coeff(b, CR_ONE)]))
    end
end

# CCustomType scaling by number goes through generic rules too, but keep these fast paths:
*(f::CCustomType, k::Number) = CCustomType(f.param_info, f.coeff*k, f.expr, f.ctype_def)
*(k::Number, f::CCustomType) = f * k

# Scale CPower by a Number
*(p::CPower, k::Number) = CPower(p.param_info, p.coeff*k, p.expr, p.exponent, Val(:nosimp))
*(k::Number, p::CPower) = p * k

# Rational interactions
*(a::CFunction, r::CRational) = CRational(pinfo(a), a*r.numer, r.denom)
*(r::CRational, a::CFunction) = CRational(pinfo(a), r.numer*a, r.denom)
*(a::CRational, b::CRational) = CRational(pinfo(a), a.numer*b.numer, a.denom*b.denom)

# generic multiply when both are non-sums
function *(a::Ta, b::Tb) where {Ta<:CFunction, Tb<:CFunction}
    _ensure_same_param_info(a, b)
    ca = coeff(a); cb = coeff(b)
    if length(ca) != 1 || length(cb) != 1
        throw(ArgumentError("Internal: multiply() called on multi-coeff CFunctions."))
    end
    if iszero(cb[1]) || iszero(ca[1])
        return zero_atom(a)
    else
        # Factor out coefficients into CProd
        return CProd(pinfo(a), ca[1]*cb[1], sort!([a/ca[1], b/cb[1]]))
    end
end

function *(a::CProd, b::CFunction)
    _ensure_same_param_info(a, b)
    cb = coeff(b)
    if iszero(cb[1])
        return zero_atom(a)
    else
        return CProd(pinfo(a), a.coeff*cb[1], sort!(vcat(a.expr, b/cb[1])))
    end
end
*(a::CFunction, b::CProd) = b*a

*(a::CSum, b::CProd) = _CSum(pinfo(a), [ x*b for x in a.expr ])
*(b::CProd, a::CSum) = a*b

*(a::CProd, b::CProd) = CProd(pinfo(a), a.coeff*b.coeff, sort!(vcat(a.expr, b.expr)))

# Exp/Log
*(a::CExp, b::CExp) = CExp(pinfo(a), a.coeff*b.coeff, a.expr + b.expr)

# CRational × CProd: push multiplication into first “scalar-like” term if present
function *(a::CRational, b::CProd)
    ind = findfirst(x -> x isa _ScalarLike, b.expr)
    if ind === nothing
        return CRational(pinfo(a), a.numer*b, a.denom)
    end
    term_ind = b.expr[ind] * a
    CProd(pinfo(b), b.coeff, vcat(b.expr[1:ind-1], term_ind, b.expr[ind+1:end]))
end
*(a::CProd, b::CRational) = b * a

# Number scaling for various types
*(a::CAbstract, b::Number) = CAbstract(a.param_info, a.coeff*b, a.index, a.exponent, a.dag)
*(a::CAtom,     b::Number) = CAtom(pinfo(a), a.coeff*b, a.var_exponents)
*(s::CSum,      b::Number) = _CSum(pinfo(s), [ x*b for x in s.expr ])
*(a::CRational, b::Number) = CRational(pinfo(a), a.numer*b, a.denom, Val(:nosimp))
*(a::CProd,     b::Number) = CProd(pinfo(a), a.coeff*b, a.expr, Val(:nosimp))
*(a::CLog,      b::Number) = CLog(pinfo(a),  b*a.coeff, a.expr, Val(:nosimp))
*(a::CExp,      b::Number) = CExp(pinfo(a),  b*a.coeff, a.expr, Val(:nosimp))
*(b::Number, a::CFunction) = a * b

# Vector/Matrix scaling and mixing
*(k::Number, v::CVector)    = CVector(pinfo(v), v.coeff*k, v.expr; row=v.row)
*(v::CVector, k::Number)    = k * v
*(s::CFunction, v::CVector) = CVector(pinfo(v), v.coeff, [s * x for x in v.expr]; row=v.row)
*(v::CVector, s::CFunction) = CVector(pinfo(v), v.coeff, [x * s for x in v.expr]; row=v.row)

*(k::Number, A::CMatrix)    = CMatrix(pinfo(A), A.coeff*k, A.expr)
*(A::CMatrix, k::Number)    = k * A
*(s::CFunction, A::CMatrix) = CMatrix(pinfo(A), A.coeff, [s * x for x in A.expr])
*(A::CMatrix, s::CFunction) = CMatrix(pinfo(A), A.coeff, [x * s for x in A.expr])

function *(u::CVector, v::CVector)
    mu, nu = size(u); mv, nv = size(v)
    if mu == 1 && nv == 1 && nu == mv
        s = _CSum(pinfo(u), [ u.expr[i] * v.expr[i] for i in 1:nu ])
        return s * (u.coeff * v.coeff)
    elseif nu == 1 && mv == 1
        n = mu; m = nv
        return CMatrix(pinfo(u), u.coeff * v.coeff,
                       [ u.expr[i] * v.expr[j] for i in 1:n, j in 1:m ])
    else
        error("Vector * Vector mismatch: size(u)=$(size(u)), size(v)=$(size(v)).")
    end
end

function *(A::CMatrix, v::CVector)
    m, n = size(A.expr)
    size(v) == (n, 1) || error("A*v mismatch: size(A)=($m,$n), size(v)=$(size(v)).")
    CVector(pinfo(A), A.coeff * v.coeff,
            [ _CSum(pinfo(A), [ A.expr[i,j] * v.expr[j] for j in 1:n ]) for i in 1:m ];
            row=false)
end

function *(u::CVector, A::CMatrix)
    m, n = size(A.expr)
    size(u) == (1, m) || error("u*A mismatch: size(u)=$(size(u)), size(A)=($m,$n).")
    CVector(pinfo(A), u.coeff * A.coeff,
            [ _CSum(pinfo(A), [ u.expr[i] * A.expr[i,j] for i in 1:m ]) for j in 1:n ];
            row=true)
end

function *(A::CMatrix, B::CMatrix)
    m, n = size(A.expr); n2, p = size(B.expr)
    n == n2 || error("Inner dimensions must match for matrix multiplication.")
    CMatrix(pinfo(A), A.coeff * B.coeff,
            [ _CSum(pinfo(A), [ A.expr[i,k] * B.expr[k,j] for k in 1:n ])
              for i in 1:m, j in 1:p ])
end

# -------- Hadamard (broadcasted * ) ----------------------------------------
# Helpers
_check_vec_same_shape(u::CVector, v::CVector) =
    (u.row == v.row && length(u.expr) == length(v.expr)) ||
    error("Hadamard requires same vector orientation and length.")

_hadamard(u::CVector, v::CVector) = begin
    _check_vec_same_shape(u, v)
    CVector(pinfo(u), u.coeff * v.coeff, [ ui * vi for (ui, vi) in zip(u.expr, v.expr) ]; row=u.row)
end

_hadamard(A::CMatrix, B::CMatrix) = begin
    size(A.expr) == size(B.expr) || error("Hadamard requires same matrix size.")
    m, n = size(A.expr)
    CMatrix(pinfo(A), A.coeff * B.coeff, [ A.expr[i,j] * B.expr[i,j] for i in 1:m, j in 1:n ])
end

Base.Broadcast.broadcasted(::typeof(*), u::CVector, v::CVector) = _hadamard(u, v)
Base.Broadcast.broadcasted(::typeof(*), v::CVector, k::Number)  = CVector(pinfo(v), v.coeff * k, v.expr; row=v.row)
Base.Broadcast.broadcasted(::typeof(*), k::Number, v::CVector)  = CVector(pinfo(v), v.coeff * k, v.expr; row=v.row)
Base.Broadcast.broadcasted(::typeof(*), v::CVector, s::CFunction) = CVector(pinfo(v), v.coeff, [ x * s for x in v.expr ]; row=v.row)
Base.Broadcast.broadcasted(::typeof(*), s::CFunction, v::CVector) = CVector(pinfo(v), v.coeff, [ s * x for x in v.expr ]; row=v.row)

Base.Broadcast.broadcasted(::typeof(*), A::CMatrix, B::CMatrix) = _hadamard(A, B)
Base.Broadcast.broadcasted(::typeof(*), A::CMatrix, k::Number)  = CMatrix(pinfo(A), A.coeff * k, A.expr)
Base.Broadcast.broadcasted(::typeof(*), k::Number, A::CMatrix)  = CMatrix(pinfo(A), A.coeff * k, A.expr)
Base.Broadcast.broadcasted(::typeof(*), A::CMatrix, s::CFunction) = CMatrix(pinfo(A), A.coeff, [ x * s for x in A.expr ])
Base.Broadcast.broadcasted(::typeof(*), s::CFunction, A::CMatrix) = CMatrix(pinfo(A), A.coeff, [ s * x for x in A.expr ])

# -------- division ----------------------------------------------------------
/(a::CAtom,     b::Number) = CAtom(pinfo(a), a.coeff/b, a.var_exponents)
(/)(a::CAbstract, b::Number) = CAbstract(a.param_info, a.coeff/b, a.index, a.exponent, a.dag)
(/)(f::CCustomType, b::Number) = CCustomType(f.param_info, f.coeff/b, f.expr, f.ctype_def)

function /(A::CSum, b::CAtom)
    _CSum(pinfo(A), [ x/b for x in A.expr ])
end
/(A::CSum, b::CAbstract) = _CSum(pinfo(A), [ x/b for x in A.expr ])

# abstract / abstract -> merge exponents if same index & dag
function /(a::CAbstract, b::CAbstract)
    _ensure_same_param_info(a, b)
    if a.index == b.index && a.dag == b.dag
        return CAbstract(a.param_info, a.coeff/b.coeff, a.index, a.exponent - b.exponent, a.dag)
    else
        return CRational(pinfo(a), a, b)
    end
end

/(a::CAtom, b::CAtom) = begin
    _ensure_same_param_info(a, b)
    CAtom(pinfo(a), a.coeff/b.coeff, a.var_exponents .- b.var_exponents)
end

/(r::CRational, a::CAtom)     = CRational(pinfo(r), r.numer, r.denom*a)
(/)(r::CRational, a::CAbstract)= CRational(pinfo(r), r.numer, r.denom*a)

# generic to rational
/(A::CSum, B::CSum)   = CRational(pinfo(A), A, B)
/(A::CAtom, B::CSum)  = CRational(pinfo(A), A, B)
(/)(a::T1, b::T2) where {T1<:CFunction, T2<:CFunction} = CRational(pinfo(a), a, b)

# rational & sums
/(a::CAtom, r::CRational)      = CRational(pinfo(a), a*r.denom, r.numer)
(/)(a::CRational, b::CRational)= CRational(pinfo(a), a.numer*b.denom, a.denom*b.numer)
(/)(a::CRational, b::CSum)     = CRational(pinfo(a), a.numer, a.denom*b)
(/)(a::CSum, b::CRational)     = CRational(pinfo(a), a*b.denom, b.numer)
(/)(a::CRational, b::Number)   = CRational(pinfo(a), a.numer, a.denom*b)

# Numbers
function /(a::CSum, n::Number)
    @assert !iszero(n) "Cannot divide by zero"
    _CSum(pinfo(a), [ x/n for x in a.expr ])
end
function /(a::CProd, n::Number)
    @assert !iszero(n) "Cannot divide by zero"
    CProd(pinfo(a), a.coeff/n, a.expr, Val(:nosimp))
end
function /(a::CLog, n::Number)
    @assert !iszero(n) "Cannot divide by zero"
    CLog(pinfo(a), a.coeff/n, a.expr, Val(:nosimp))
end
function /(a::CExp, n::Number)
    @assert !iszero(n) "Cannot divide by zero"
    CExp(pinfo(a), a.coeff/n, a.expr, Val(:nosimp))
end
/(n::Number, a::T) where {T<:CFunction} = num_atom(a, n) / a

function /(a::CProd, b::_ScalarLike)
    _ensure_same_param_info(a, b)
    ind = findfirst(x -> x isa _ScalarLike, a.expr)
    if ind === nothing
        return CRational(pinfo(a), a, b, Val(:nosimp))
    end
    cb = coeff(b)[1]
    term_ind = a.expr[ind] / b * cb
    CProd(pinfo(a), a.coeff / cb, vcat(a.expr[1:ind-1], term_ind, a.expr[ind+1:end]))
end

/(v::CVector, k::Number) = CVector(pinfo(v), v.coeff / k, v.expr; row=v.row)
/(v::CVector, s::CFunction) = CVector(pinfo(v), v.coeff, [x / s for x in v.expr]; row=v.row)
/(A::CMatrix, k::Number) = CMatrix(pinfo(A), A.coeff / k, A.expr)
/(A::CMatrix, s::CFunction) = CMatrix(pinfo(A), A.coeff, [x / s for x in A.expr])

# CPower / number
/(p::CPower, k::Number) = CPower(p.param_info, p.coeff/k, p.expr, p.exponent, Val(:nosimp))

# -------- exponentiation ----------------------------------------------------
^(x::CFunction, q::Rational{Int}) = CPower(pinfo(x), x, q)
^(x::CFunction, n::Int)           = CPower(pinfo(x), x, n//1)

^(A::CAtom, n::Int)     = CAtom(pinfo(A), A.coeff^n, A.var_exponents .* n)
^(a::CAbstract, n::Int) = CAbstract(a.param_info, a.coeff^n, a.index, a.exponent*n, a.dag)
^(a::CRational, n::Int) = CRational(pinfo(a), a.numer^n, a.denom^n)
^(a::CProd, n::Int)     = CProd(pinfo(a), a.coeff^n, [ x^n for x in a.expr ])
^(a::CExp, n::Int)      = CExp(pinfo(a), a.coeff^n, a.expr * n)

function ^(p::CPower, n::Int)
    n >= 0 && return CPower(p.param_info, p.coeff^n, p.expr, p.exponent*n, Val(:nosimp))
    # negative: return 1 / p^(-n)
    return CRational(p.param_info, one_atom(p), p^(-n))
end

function ^(s::CSum, n::Int)
    if n < 0
        return CRational(pinfo(s), one_atom(s), s^(-n))
    elseif n == 0
        return one_atom(s)
    elseif n == 1
        return s
    else
        result = s
        for _ in 2:n
            result = expand_prod(result, s)
        end
        return result
    end
end

# Helper: distributive law for sums
expand_prod(a::CSum, b::CSum)      = _CSum(pinfo(a), vcat([x*y for x in a.expr, y in b.expr]...))
expand_prod(a::CSum, b::CFunction) = _CSum(pinfo(a), [ x*b for x in a.expr ])
expand_prod(a::CFunction, b::CSum) = _CSum(pinfo(b), [ a*y for y in b.expr ])

# -------- roots -------------------------------------------------------------
include("ComplexRationals_sqrt.jl")
sqrt(x::CFunction) = CPower(pinfo(x), x, 1//2)
function sqrt(x::CAtom)
    if has_rational_sqrt(x.coeff) && all(e -> e % 2 == 0, x.var_exponents)
        return CAtom(pinfo(x), principal_sqrt(x.coeff), x.var_exponents .÷ 2)
    else
        return CPower(pinfo(x), x, 1//2)
    end
end

# -------- inverses & adjoints ----------------------------------------------
inv(a::CAtom)     = CAtom(pinfo(a), inv(a.coeff), a.var_exponents .* (-1))
inv(a::CAbstract) = CAbstract(a.param_info, inv(a.coeff), a.index, -a.exponent, a.dag)
inv(a::CSum)      = CRational(pinfo(a), one_atom(a), a)
inv(a::CRational) = CRational(pinfo(a), inv(a.numer), inv(a.denom))
inv(p::CPower)    = CPower(p.param_info, inv(p.coeff), p.expr, -p.exponent, Val(:nosimp))

function adjoint(f::CAtom)::CAtom
    CAtom(pinfo(f), conj(f.coeff), copy(f.var_exponents))
end
function adjoint(a::CAbstract)
    CAbstract(a.param_info, conj(a.coeff), a.index, a.exponent, !a.dag)
end
function adjoint(f::CCustomType)
    CCustomType(f.param_info, conj(f.coeff), adjoint.(f.expr), f.ctype_def)
end

adjoint(f::CSum)      = _CSum(pinfo(f), [adjoint(t) for t in f.expr])
adjoint(f::CRational) = CRational(pinfo(f), adjoint(f.numer), adjoint(f.denom))
adjoint(p::CPower)    = CPower(p.param_info, conj(p.coeff), adjoint(p.expr), p.exponent, Val(:nosimp))
conj(f::CFunction)    = adjoint(f)

transpose(v::CVector) = CVector(pinfo(v), v.coeff, v.expr; row = !v.row)
adjoint(v::CVector)   = CVector(pinfo(v), conj(v.coeff), adjoint.(v.expr); row = !v.row)
conj(v::CVector)      = CVector(pinfo(v), conj(v.coeff), adjoint.(v.expr); row = v.row)

transpose(A::CMatrix) = CMatrix(pinfo(A), A.coeff, permutedims(A.expr))
adjoint(A::CMatrix)   = CMatrix(pinfo(A), conj(A.coeff), adjoint.(permutedims(A.expr)))
conj(A::CMatrix)      = CMatrix(pinfo(A), conj(A.coeff), adjoint.(A.expr))

# -------- equality ----------------------------------------------------------
==(a::S, b::T) where {S<:CFunction, T<:CFunction} = false
==(a::S, b::S) where {S<:CFunction} = error("Equality not implemented for type $S")

==(a::CAtom, b::CAtom)           = (a.coeff == b.coeff && a.var_exponents == b.var_exponents)
==(a::CAbstract, b::CAbstract)   = (a.index == b.index && a.dag == b.dag && a.exponent == b.exponent && a.coeff == b.coeff)
==(a::CSum, b::CSum)             = (length(a) == length(b) && a.expr == b.expr)
==(a::CRational, b::CRational)   = (a.numer == b.numer && a.denom == b.denom)
==(a::CProd, b::CProd)           = (length(a.expr) == length(b.expr) && a.coeff == b.coeff && a.expr == b.expr)
==(a::CExp, b::CExp)             = (a.coeff == b.coeff && a.expr == b.expr)
==(a::CLog, b::CLog)             = (a.coeff == b.coeff && a.expr == b.expr)
==(a::CPower, b::CPower)         = (a.exponent == b.exponent && a.coeff == b.coeff && a.expr == b.expr)
==(f::CCustomType, g::CCustomType) = (f.ctype_def === g.ctype_def && f.coeff == g.coeff && f.expr == g.expr)

==(u::CVector, v::CVector) =
    (u.row == v.row) &&
    (length(u.expr) == length(v.expr)) &&
    all( (u.coeff*u.expr[i]) == (v.coeff*v.expr[i]) for i in eachindex(u.expr) )

==(A::CMatrix, B::CMatrix) =
    (size(A.expr) == size(B.expr)) &&
    all( (A.coeff*A.expr[i]) == (B.coeff*B.expr[i]) for i in eachindex(A.expr) )
##############################################################################
