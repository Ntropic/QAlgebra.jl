module FFunctions

using ..StringUtils
using ComplexRationals

export FFunction, FAtom, FSum, FRational, simplify, isnumeric, iszero, max_exponents, build_xpows, evaluate, to_string

"""
    FFunction

Abstract supertype for symbolic functions representing atoms (`FAtom`), sums (`FSum`), and rationals (`FRational`).
"""
abstract type FFunction end

"""
    FAtom(coeff::Int, var_exponents::Vector{Int})
    FAtom(coeff::Rational, var_exponents::Vector{Int})
    FAtom(coeff::ComplexRational, var_exponents::Vector{Int})

A single term with a complex‐rational coefficient and integer exponents for each variable.
- The `Int` and `Rational` constructors wrap the coefficient into a `ComplexRational`.
- `var_exponents[j]` is the exponent of variable _j_.
"""
mutable struct FAtom      <: FFunction
    coeff::ComplexRational
    var_exponents::Vector{Int}
    function FAtom(coeff::Int, var_exponents::Vector{Int})
        c = ComplexRational(coeff, 0, 1)
        return new(c, var_exponents)
    end
    function FAtom(coeff::Rational, var_exponents::Vector{Int})
        c = ComplexRational(numerator(coeff), 0, denominator(coeff))
        return new(c, var_exponents)
    end
    function FAtom(coeff::Complex, var_exponents::Vector{Int})
        c = crationalize(coeff)
        return new(c, var_exponents)
    end
    function FAtom(coeff::ComplexRational, var_exponents::Vector{Int})
        return new(coeff, var_exponents)
    end
end

"""
    FSum(terms::AbstractVector{<:FFunction})
    FSum(ts::FFunction...)

Constructs a sum of `FFunction` terms.
- Flattens any nested `FSum` automatically.
- Variadic form `FSum(a, b, c)` is provided for convenience.
"""
mutable struct FSum       <: FFunction
    terms::Vector{FFunction}
    # (inner) constructor for a Vector{<:FFunction>, with flattening
    function FSum(ts::AbstractVector{<:FFunction})
        flat = FFunction[]
        for t in ts
            if t isa FSum
                append!(flat, (t::FSum).terms)   # flatten nested sums
            else
                push!(flat, t)
            end
        end
        new(flat)
    end
end
# varargs “outer” constructor so you can call FSum(a,b,c) directly
FSum(ts::FFunction...) = FSum(collect(ts))

"""
    FRational(numer::FSum, denom::FSum)

Represents a rational function with numerator `numer` and denominator `denom`, both sums of `FFunction` terms.
"""
mutable struct FRational  <: FFunction
    numer::FSum
    denom::FSum
end

_terms(f::FFunction) = f isa FSum ? (f::FSum).terms : [f]

import Base: isless
function isless(a::FAtom, b::FAtom)
    a.var_exponents < b.var_exponents    # Vector{Int} has lex order
end


import Base: iszero, isempty
"""
    iszero(a::FAtom)     -> Bool
    iszero(s::FSum)      -> Bool
    iszero(r::FRational) -> Bool

Returns `true` if the expression is identically zero:
- **Atom**: zero coefficient.
- **Sum**: all terms zero or empty.
- **Rational**: zero numerator.
"""
iszero(a::FAtom)        = iszero(a.coeff)
iszero(s::FSum)         = isempty(s.terms) || all(iszero, s.terms)
iszero(r::FRational)    = iszero(r.numer)

isempty(s::FSum)        = isempty(s.terms)

"""
    isnumeric(a::FAtom)    -> Bool
    isnumeric(s::FSum)     -> Bool
    isnumeric(r::FRational) -> Bool

Returns `true` if the expression contains no variables (i.e., all exponents are zero in atoms, and both numerator and denominator are numeric sums).
"""
isnumeric(a::FAtom)    = all(e->e==0, a.var_exponents)
isnumeric(s::FSum)     = all(isnumeric, s.terms)
isnumeric(r::FRational)= isnumeric(r.numer) && isnumeric(r.denom)

allnegative(a::FAtom) = is_negative(a.coeff)
allnegative(s::FSum)  = !isempty(s.terms) && all(allnegative, s.terms)
allnegative(r::FRational) = allnegative(r.numer)

min_exponents(a::FAtom)    = a.var_exponents
function min_exponents(s::FSum)
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
function min_exponents(r::FRational)
    min_vals = min_exponents(r.numer)
    min_vals2 = min_exponents(r.denom)
    return min.(min_vals, min_vals2)
end


import Base: length, getindex, iterate, deleteat!, reverse
function length(p::FSum)::Int
    return length(p.terms)
end
getindex(p::FSum, i::Int) = p.terms[i]
iterate(p::FSum, state=1) = state > length(p.terms) ? nothing : (p.terms[state], state + 1)
deleteat!(p::FSum, i::Int) = FSum(deleteat!(p.terms, i))
reverse(q::FSum) = FSum(reverse(q.terms))
dims(q::FAtom) = length(q.var_exponents)
dims(q::FSum) = dims(q.terms[1])
dims(q::FRational) = dims(q.numer) 

import Base: +, -, *, /, ^, ==

# addition always builds a flat sum
+(a::FFunction) =  a 
+(a::FFunction, b::FFunction) = FSum( vcat(_terms(a), _terms(b)) )

# unary minus & subtraction
import Base: -, +

-(a::FAtom) = FAtom(-a.coeff, a.var_exponents)
-(s::FSum) = FSum([ -t for t in s.terms ])
-(r::FRational) = FRational(-r.numer, r.denom)
-(a::FFunction, b::FFunction) = a + (-b)


# distribute * over sums
*(A::FSum, B::FSum)       = FSum([ x*y for x in A.terms,   y in B.terms ])
*(A::FSum, b::FFunction)  = FSum([ x*b for x in A.terms ])
*(a::FFunction, B::FSum)  = FSum([ a*y for y in B.terms ])

# atom‐level ×
*(a::FAtom, b::FAtom)     = FAtom(crationalize(a.coeff*b.coeff), a.var_exponents .+ b.var_exponents)
*(a::FAtom, r::FRational) = FRational(FSum(a)*r.numer, r.denom)
*(r::FRational, a::FAtom) = FRational(r.numer*FSum(a), r.denom)
*(a::FRational, b::FRational) = FRational(a.numer*b.numer, a.denom*b.denom)
*(a::FSum, b::FRational) = FRational(a*b.numer, b.denom)
*(b::FRational, a::FSum) = FRational(a*b.numer, b.denom)
# number 
*(a::FFunction, b::Number)  = a*FAtom(b, zeros(Int, dims(a)))
*(b::Number, a::FFunction)  = FAtom(b, zeros(Int, dims(a))) * a
multiply_one(a::FRational, b::Int) = (a.numer * b) / (a.denom * b)

# division
/(A::FSum, B::FSum) = length(B.terms)==1 ? FSum([ x/B.terms[1] for x in A.terms ]) : FRational(A, B)
/(A::FAtom, B::FSum)     = FRational(FSum(A), B)
/(A::FSum, b::FAtom)     = FSum([ x/b for x in A.terms ])
/(a::FAtom, b::FAtom)     = FAtom(a.coeff/b.coeff, a.var_exponents .- b.var_exponents)
/(a::FAtom, r::FRational) = FRational(FSum(a)*r.denom, r.numer)
/(r::FRational, a::FAtom) = FRational(r.numer, r.denom*FSum(a))
/(a::FRational, b::FRational) = FRational(a.numer*b.denom, a.denom*b.numer)
/(a::FRational, b::FSum) = FRational(a.numer, b.numer*b)
/(a::FSum, b::FRational) = FRational(a*b.denom, b.numer)
/(a::FFunction, n::Number) = a/FAtom(n, zeros(Int, dims(a)))
/(n::Number, a::FFunction) = FAtom(n, zeros(Int, dims(a))) / a

# exponentiation 
^(A::FSum, n::Int) = FSum([ x^n for x in A.terms ])
^(A::FAtom, n::Int) = FAtom(A.coeff^n, A.var_exponents .* n)
^(a::FRational, n::Int) = FRational(a.numer^n, a.denom^n)
import Base: copy 
function copy(x::FAtom)
    return FAtom(copy(x.coeff), copy(x.var_exponents))
end
function copy(x::FSum)
    return FSum(copy.(x.terms))
end
function copy(x::FRational)
    return FRational(copy(x.numer), copy(x.denom))
end

function ==(a::FAtom, b::FAtom)
    return (a.coeff == b.coeff && a.var_exponents == b.var_exponents)
end
function ==(a::FSum, b::FSum)
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
function ==(a::FRational, b::FRational)
    return (a.numer == b.numer && a.denom == b.denom)
end

function unifiable(a::FAtom, b::FAtom)::Bool
    return a.var_exponents == b.var_exponents
end
function unifiable(a::FRational, b::FRational)::Bool
    return a.denom == b.denom
end
function unifiable(a::FSum, b::FSum)::Bool
    true
end
# assume inifiable
function unify(a::FAtom, b::FAtom)::FFunction
    absum = a.coeff+b.coeff
    if iszero(absum)
        var_zeros = zeros(Int, dims(a))
        return FAtom(0, var_zeros)
    end
    return FAtom(absum, a.var_exponents)
end
function unify(a::FRational, b::FRational)::FFunction
    simple_numer = simplify(a.numer+b.numer)
    if iszero(simple_numer)
        var_zeros = zeros(Int, dims) 
        return FAtom(0, var_zeros)
    end
    return FRational(simple_numer, a.denom)
end
function unify(a::FSum, b::FSum)::FFunction
    return simplify(a+b)
end

import Base: sort!, sort

# Replace with the  custom_sort_key ? 
#function complexrational2coeffs(a::ComplexRational)::Vector{Int}
#    return [a.a, a.b, a.c]
#end

function sort_key(a::FAtom)
    return vcat(a.var_exponents, custom_sort_key(a.coeff))
end

function sort_key(s::FSum)
    keys = [sort_key(t) for t in s.terms]
    return vcat(length(keys), keys...)
end
function sort_key(r::FRational)
    dkeys = sort_key(r.denom)
    nkeys = sort_key(r.numer)
    return vcat(length(dkeys), dkeys, nkeys)
end

function sort!(s::FSum; by=sort_key)
    sort!(s.terms, by=by)
end
function sort!(s::FRational, by=sort_key)
    sort!(s.numer.terms, by=by)
    sort!(s.denom.terms, by=by)
end
function sort!(s::FAtom; by=sort_key)
    # do nothing 
    return s
end

function sort(s::FSum; by=sort_key)::FSum
    return FSum(sort(s.terms, by=by))
end
function sort(s::FRational; by=sort_key)::FRational
    numer = sort(s.numer.terms, by=by)
    denom = sort(s.denom.terms, by=by)
    return FRational(numer, denom)
end
function sort(s::FAtom; by=sort_key)::FAtom
    # do nothing 
    return s
end

issimple(f::FFunction) = error("Not implemented for type $(typeof(f)).")
issimple(f::FSum) = all(t -> t isa FAtom, f.terms)
issimple(f::FRational) = issimple(f.numer) && issimple(f.denom)
issimple(f::FAtom) = true

# divisors 
function divisors(a::FFunction)
    error("Not implemented for type $(typeof(a)).")
end
function divisors(a::FAtom)::Vector{Int}
    return [a.coeff.c]
end
function divisors(a::FSum)::Vector{Int}
    return reduce(vcat, [divisors(t) for t in a.terms])
end
function divisors(a::FRational)::Vector{Int}
    return reduce(vcat, [divisors(t) for t in [a.numer, a.denom]])
end

function vec_multiply(x::FAtom, vector::Vector{Int})::FAtom
    return FAtom(x.coeff, x.var_exponents - vector)
end
function vec_multiply(x::FSum, vector::Vector{Int})::FSum
    return FSum([vec_multiply(t, vector) for t in x.terms])
end
function vec_multiply(x::FRational, vector::Vector{Int})::FRational
    return FRational(vec_multiply(x.numer), vec_multiply(x.denom))
end

# Assumes that the divisors are 1 
function coeffs(s::FAtom)::Vector{Int}
    # if subtype complex take real and imaginary parts separately 
    return [s.coeff.a, s.coeff.b]
end
function coeffs(s::FSum)::Vector{Int}
    return vcat([coeffs(t) for t in s.terms]...)
end
function coeffs(s::FRational)::Vector{Int}
    return vcat(coeffs(s.numer), coeffs(s.denom))
end
import Base: gcd
function gcd(s::FFunction)
    return gcd(coeffs(s))   
end
function gcd(s::FFunction, t::FFunction)
    return gcd(vcat(coeffs(s), coeffs(t)))
end
#gcd(a*2+b+a*4)

"""
    simplify(a::FAtom) -> FAtom
    simplify(s::FSum)  -> FSum
    simplify(r::FRational) -> FFunction

Simplifies symbolic expressions:
- **Atom**: returns unchanged.
- **Sum**: recursively simplifies terms, sorts them, and combines like terms (unifies coefficients).
- **Rational**: simplifies numerator and denominator, cancels common factors in coefficients and exponents, and collapses to a sum if the denominator becomes a single term.
"""
function simplify end

function simplify(s::FAtom)
    return s 
end
function simplify(r::FRational)
    n = simplify(r.numer)
    d = simplify(r.denom)
    # remove fractions in coefficients in Rational
    curr_div = vcat(divisors(n), divisors(d))
    factor = lcm(curr_div...)
    n = factor * n
    d = factor * d
    # common denominator
    factor = gcd(n, d)
    n = n / factor
    d = d / factor

    min_n = min_exponents(n)
    min_d = min_exponents(d)
    min_vals = min.(min_n, min_d)
    n = vec_multiply(n, min_vals)
    d = vec_multiply(d, min_vals)
    if allnegative(d) || (get_default(:FIRST_MODE) && allnegative(d[1]))   # prefer negatives on numerator
        n = -n
        d = -d
    end
    # if denom now has exactly one term, collapse back to a sum
    if length(d.terms) == 0
        error("Divinding by zero")
    elseif length(d.terms) == 1
        # use our /‐overload to divide each term in n by d.terms[1]
        return n / d.terms[1]
    else
        return FRational(n, d)
    end
end

function simplify(s::FSum)
    # first simplify the lower levels 
    elements = [simplify(e) for e in s.terms]
    s = FSum(elements)
    sort!(s)
    i = 1
    elements = s.terms
    new_elements = FFunction[]
    curr_element = elements[1]
    for i in 2:length(elements)
        if typeof(curr_element) == typeof(elements[i]) 
            if unifiable(curr_element, elements[i])
                curr_element = unify(curr_element, elements[i])
            else
                if !iszero(curr_element)
                    push!(new_elements, curr_element)
                end
                curr_element = elements[i]
            end 
        else
            if !iszero(curr_element)
                push!(new_elements, curr_element)
            end
            curr_element = elements[i]
        end
    end
    if !iszero(curr_element)
        push!(new_elements, curr_element)
    end
    if length(new_elements) == 0 
        push!(new_elements, FAtom(0, zeros(Int, dims(s))))
    end
    return FSum(new_elements)
end


"""
    max_exponents(a::FAtom)    -> Vector{Int}
    max_exponents(s::FSum)      -> Vector{Int}
    max_exponents(r::FRational) -> Vector{Int}

Computes the elementwise maximum of variable exponents:
- **Atom**: absolute value of its own exponents.
- **Sum**: maximum across all terms.
- **Rational**: maximum between numerator and denominator exponents.
"""
max_exponents(a::FAtom) = abs.(a.var_exponents)
function max_exponents(s::FSum)
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
function max_exponents(r::FRational)
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
    evaluate(a::FAtom, xpows::Vector{Vector})       -> Number
    evaluate(s::FSum, xpows::Vector{Vector})        -> Number
    evaluate(r::FRational, xpows::Vector{Vector})   -> Number

Evaluates the expression using precomputed powers `xpows` (as returned by `build_xpows`) for faster repeated evaluation.
"""
function evaluate(a::FAtom, x::AbstractVector{<:Number})
    a.coeff * prod(x[i]^a.var_exponents[i] for i in eachindex(a.var_exponents))
end
function evaluate(s::FSum, x::AbstractVector{<:Number})
    sum(evaluate(t, x) for t in s.terms)
end
function evaluate(r::FRational, x::AbstractVector{<:Number})
    evaluate(r.numer, x) / evaluate(r.denom, x)
end

# Evaluate but with xpows
function evaluate(a::FAtom, xpows::Vector{Vector})
    a.coeff * prod(
        a.var_exponents[i] >= 0 ?
            xpows[i][a.var_exponents[i] + 1] :
            1 / xpows[i][abs(a.var_exponents[i]) + 1]
        for i in eachindex(a.var_exponents)
    )
end
function evaluate(s::FSum, xpows::Vector{Vector})
    return sum(evaluate(t, xpows) for t in s.terms)
end
function evaluate(r::FRational, xpows::Vector{Vector})
    return evaluate(r.numer, xpows) / evaluate(r.denom, xpows)
end

function sign_string(c::ComplexRational, do_latex::Bool=false)::Tuple{Bool, String}
    if is_negative(c)
        return (true, string(c, do_latex=do_latex)[2:end])
    else
        return (false, string(c, do_latex=do_latex))
    end
end
is_abs_one(c::ComplexRational)::Bool = (abs(c.a) == abs(c.b))
function is_abs_one(c::FFunction)
    if isnumeric(c)
        if isa(c, FAtom)
            return is_abs_one(c.coeff)
        elseif isa(c, FSum)
            if length(c) == 1
                return is_abs_one(c[1])
            else
                return false
            end
        else
            return false
        end
        return true
    else
        return false
    end
end

# --- Generic string constructor ---
function stringer(f::FFunction)::Tuple{String,String}
    error("No fallback method for FFunction type: $(typeof(f))")
end

# --- FAtom ---
function stringer(a::FAtom)::Tuple{Bool,String}  # true = minus, false = plus
    if isnumeric(a)
        return sign_string(a.coeff) 
    else
        c = a.coeff
        vec = string(a.var_exponents)
        sign, c_str = sign_string(c) 
        if is_abs_one(c)
            return (sign, vec)
        else
            return sign,  c_str*"*"*vec
        end
    end
end

# --- FSum ---
function stringer(s::FSum; braced::Bool=false) ::Tuple{Bool,String}
    # braced specifies whether the terms will be grouped, so that an external sign is needed
    s = simplify(s)
    terms = s.terms
    if braced 
        if allnegative(s) || (get_default(:FIRST_MODE) && allnegative(s[1]))
            # if all negative and braced, we can just negate the whole thing
            sig = true
            _, body = stringer(-s)
            return sig, body
        else
            sig = false
            body = join(stringer(t) for t in terms)
        end
    end

    # process each term into (sign, body)
    parts = String[]
    for (i, t) in enumerate(terms)
        sig, body = stringer(t)
        if i == 1
            # first term keeps its sign, but no space if positive
            push!(parts, sig == true ? "-$body" : body)
        else
            push!(parts, sig == true ? "-$body" : "+$body")
        end
    end
    return false, join(parts, "")
end


# --- FRational ---
function stringer(r::FRational)::Tuple{Bool,String}
    r = simplify(r)
    n = r.numer
    d = r.denom

    # check first term in numerator
    n_sig, n_body = stringer(n, braced=true)
    d_sig, d_body = stringer(d)  # we ignore sign of denom

    if length(n) > 1
        n_body = "($n_body)"
    end
    if length(d) > 1
        d_body = "($d_body)"
    end

    return n_sig, "$n_body/$d_body"
end

# --- Show ---
import Base: show

function show(io::IO, f::FFunction)
    sig, body = stringer(f)
    sig_str = sig ? "-" : ""
    print(io, sig_str * body)
end

# Generic fallback
function stringer(f::FFunction, vars::Vector{String}; do_latex::Bool=false)
    error("No stringer method for type $(typeof(f)) with variable names")
end

# --- FAtom with variable names ---
function stringer(a::FAtom, vars::Vector{String}; do_latex::Bool=false)
    exps = a.var_exponents
    @assert length(vars) == length(exps) "Number of symbols must match number of variables"

    if isnumeric(a)
        return sign_string(a.coeff) 
    else
        # build the variable part
        varparts = String[]
        for (i, e) in enumerate(exps)
            if e == 0
                continue
            elseif do_latex
                push!(varparts, e == 1 ? vars[i] : "$(vars[i])^{$e}")
            else
                push!(varparts, e == 1 ? vars[i] : "$(vars[i])"*str2sup(string(e)))
            end
        end

        var_str = isempty(varparts) ? "" : join(varparts, do_latex ? "\\cdot " : "")
        c = a.coeff
        sign, c_str = sign_string(c) 
        if is_abs_one(c)
            return (sign, vec)
        else
            return sign,  c_str*"*"*var_str
        end
    end
end

# --- FSum with variable names ---
function stringer(s::FSum, vars::Vector{String}; do_latex::Bool=false, braced::Bool=false)
    s = simplify(s)
    terms = s.terms
    if isempty(terms)
        return false, "0"
    end

    parts = String[]
    for (i, t) in enumerate(terms)
        sig, body = stringer(t, vars; do_latex=do_latex)
        if i == 1
            push!(parts, sig ? "-" * body : body)
        else
            push!(parts, sig ? "-" * body : "+" * body)
        end
    end

    out = join(parts, "")
    return false, out
end

# --- FRational with variable names ---
function stringer(r::FRational, vars::Vector{String}; do_latex::Bool=false)
    r = simplify(r)
    n = r.numer
    d = r.denom

    n_sig, n_str = stringer(n, vars; do_latex=do_latex, braced=true)
    _, d_str = stringer(d, vars; do_latex=do_latex, braced=true)

    if do_latex
        return n_sig, "\\frac{$n_str}{$d_str}"
    else
        n_wrapped = length(n.terms) > 1 ? "($n_str)" : n_str
        d_wrapped = length(d.terms) > 1 ? "($d_str)" : d_str
        return n_sig, "$n_wrapped/$d_wrapped"
    end
end

"""
    to_string(f::FFunction, vars::Vector{String};
              do_latex::Bool = false,
              braced::Bool = false,
              optional_sign::Bool = true) -> String

Converts an `FFunction` into a human‐readable string using variable names in `vars`.
- If `do_latex=true`, uses LaTeX syntax (e.g. `\frac{}` and superscripts).
- If `braced=true`, wraps sums in parentheses.
- If `optional_sign=false`, always prefixes a plus or minus sign.
"""
function to_stringer(f::FFunction, vars::Vector{String}; do_latex::Bool=false, braced::Bool=false)::Tuple{Bool, String}
    if braced && f isa FSum && length(f) > 1
        sig, body = stringer(f, vars; do_latex=do_latex, braced=braced)
    else
        sig, body = stringer(f, vars; do_latex=do_latex)
    end
    # Apply braces to sums if requested
    if braced && f isa FSum && length(f) > 1
        body = do_latex ? "\\left( $body \\right)" : "($body)"
    end
    return sig, body
end
function to_string(f::FFunction, vars::Vector{String}; do_latex::Bool=false, braced::Bool=false, optional_sign::Bool=true)::String
    sig, body = to_stringer(f, vars; do_latex=do_latex, braced=braced)

    if sig
        return "-" * body
    else
        if optional_sign
            return body
        end
        return "+" *body
    end
end

end # module FFunctions