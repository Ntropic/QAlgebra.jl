function simplify end

simplify(s::CAtom) = s
simplify(s::CAbstract) = s
simplify(s::CCustomType) = s
simplify(s::CSum)  = simplify_CSum(s.param_info, s.expr)
simplify(r::CRational) = simplify_CRational(r.param_info, r.numer, r.denom)
simplify(e::CExp)  = simplify_CExp(e.param_info, e.coeff, e.expr)
simplify(l::CLog)  = simplify_CLog(l.param_info, l.coeff, l.expr)
simplify(p::CProd) = simplify_CProd(p.param_info, p.coeff, p.terms)
simplify(v::CVector) = CVector(v.param_info, v.coeff, simplify.(v.expr); row=v.row)
simplify(M::CMatrix) = CMatrix(M.param_info, M.coeff, reshape(simplify.(M.expr[:]), size(M.expr)))
simplify(p::CPower)  = simplify_CPower(p.param_info, p.coeff, simplify(p.param_info, ), p.exponent)

# helpers
@inline _isint(q::Rational{Int}) = denominator(q) == 1
@inline _one_atom(param_info::ParameterInfo) = CAtom(param_info, ComplexRational(1,0,1), zeros(Int, param_info.dims))
@inline _const_atom(param_info::ParameterInfo, c::ComplexRational) = CAtom(param_info, c, zeros(Int, param_info.dims))

function simplify_CPower(param_info::ParameterInfo, coeff::ComplexRational, x::CFunction, q::Rational{Int})
    x = simplify(x)  # normalize inner first

    # x^0 → 1; x^1 → x; 1^q → 1
    if iszero(q)
        return _const_atom(param_info, coeff)
    elseif _isint(q) && q == 1//1
        return coeff * x
    elseif isone(x)
        return _const_atom(param_info, coeff)
    end

    # negative exponent → put in denominator
    if q < 0
        num = _const_atom(param_info, coeff)
        den = CPower(param_info, ComplexRational(1,0,1), x, -q, Val(:nosimp))
        return simplify_CRational(num, den)
    end

    # q > 0 cases
    if x isa CRational
        # (N/D)^q = (N^q)/(D^q); outer coeff multiplies numerator
        Nq = CPower(param_info, ComplexRational(1,0,1), x.numer, q)
        Dq = CPower(param_info, ComplexRational(1,0,1), x.denom, q)
        return simplify_CRational(param_info, coeff * Nq, Dq)

    elseif x isa CProd
        if _isint(q)
            n = Int(q)
            tP = [ t^n for t in x.terms ]
            return simplify( (coeff * (x.coeff^n)) * CProd(tP) )
        else
            cP = CPower(param_info, ComplexRational(1,0,1), _const_atom(param_info, x.coeff), q)
            tP = [ CPower(param_info, ComplexRational(1,0,1), t, q) for t in x.terms ]
            return simplify( coeff * (cP * CProd(tP)) )
        end

    elseif x isa CAtom
        if _isint(q)
            n = Int(q)
            return simplify( coeff * CAtom(x.coeff^n, x.var_exponents .* n) )
        else
            return CPower(param_info, coeff, x, q, Val(:nosimp))
        end

    elseif x isa CExp
        # (c*exp(y))^q = c^q * exp(q*y)
        if _isint(q)
            n = Int(q)
            return simplify( (x.coeff^n * coeff) * CExp(ComplexRational(1,0,1), n * x.x) )
        else
            cP = CPower(param_info, ComplexRational(1,0,1), _const_atom(param_info, x.coeff), q)
            return simplify( coeff * ( cP * CExp(param_info, ComplexRational(1,0,1), q * x.x) ) )
        end

    elseif x isa CPower
        # ( (y)^a )^b = y^(a*b); keep coefficient separately if b not integer
        ab = x.exponent * q
        if _isint(q)
            n = Int(q)
            return simplify( CPower(param_info, coeff * (x.coeff^n), x.x, ab) )
        else
            cP = CPower(param_info, ComplexRational(1,0,1), _const_atom(param_info, x.coeff), q)
            return simplify( coeff * ( cP * CPower(param_info, ComplexRational(1,0,1), x.x, ab) ) )
        end
    else
        return CPower(param_info, coeff, x, q, Val(:nosimp))
    end
end


function simplify_CSum(param_info::ParameterInfo, elements::AbstractVector{<:CFunction})
    sort!(elements)
    i = 1
    new_elements = CFunction[]
    curr_element = elements[1]
    for i in 2:length(elements)
        if typeof(curr_element) == typeof(elements[i]) 
            if addable(curr_element, elements[i])
                curr_element = unify_add(curr_element, elements[i])
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
        push!(new_elements, CAtom(0, zeros(Int, dims(elements[1]))))
    elseif length(new_elements) == 1
        return new_elements[1]
    end
    return _CSum(param_info, new_elements, Val(:nosimp))
end

function simplify_CRational(param_info::ParameterInfo, n::CFunction, d::CFunction)
    # remove fractions in coefficients in Rational
    if iszero(n)
        return n 
    end
    #println("n=$n, d=$d")
    curr_div = vcat(divisors(n), divisors(d))
    #println(curr_div)
    factor = lcm(curr_div...)
    #println(factor)
    n = factor * n
    d = factor * d
    # common denominator
    factor = gcd(n, d)
    n = n / factor
    d = d / factor

    min_n = min_exponents_simple(n)
    min_d = min_exponents_simple(d)
    min_vals = min.(min_n, min_d)
    if any(min_vals .> 0)
        n = vec_multiply(n, .-min_vals)
        d = vec_multiply(d, .-min_vals)
    end
    if allnegative(d) || (FLIP_IF_FIRST_TERM_NEGATIVE  && firstnegative(d))   # prefer negatives on numerator
        n = -n
        d = -d
    end
    # if denom now has exactly one term, collapse back to a sum
    if length(d) == 0 || iszero(d)
        error("Dividing by zero: n=$n, d=$d.")
    elseif isa(d, CAtom)
        if isnumeric(d)
            return n/d.coeff
        end
        return CRational(param_info, n, d, Val(:nosimp))
    elseif isa(d, CRational)
        return (n*d.denom) / d.numer
    else
        return CRational(param_info, n, d, Val(:nosimp))
    end
end

function simplify_CExp(param_info::ParameterInfo, coeff::ComplexRational, x::CFunction)
    # exp(log(y)) ⇒ y
    if x isa CLog
        return x.x * coeff 
    # exp(0) ⇒ 1
    elseif iszero(x)
        return CAtom(param_info, coeff, zeros(Int, dims(x)))
    end
    return CExp(param_info, coeff, x, Val(:nosimp))
end

function simplify_CLog(param_info::ParameterInfo, coeff::ComplexRational, x::CFunction)
    # log(exp(y)) ⇒ y
    if x isa CExp
        return x.x * coeff 
    end
    # log(1) ⇒ 0
    if isone(x)
        return CAtom(param_info, ComplexRational(0,0,1), zeros(Int, dims(x)))
    end
    return CLog(param_info, coeff, x, Val(:nosimp)) 
end

# CAtom, CSum, CProd, CRational, CExp, CLog

function simplify_CProd(param_info::ParameterInfo, coeff::ComplexRational, terms::AbstractVector{<:CFunction})
    coeff, atoms, sums, rats, exps, logs, pows = collect_prod_terms(param_info, coeff, terms)

    # 1) combine atoms back into one "base" factor (coeff pulled out already)
    has_term = false
    local term::CFunction
    if !isempty(atoms)
        var_ex = zeros(Int, length(atoms[1].var_exponents))
        @inbounds for a in atoms
            var_ex .+= a.var_exponents
        end
        term = CAtom(param_info, ComplexRational(1,0,1), var_ex)
        has_term = true
    end

    # 2) multiply sums/rationals (at most one of each after reduce)
    if !isempty(sums)
        s = length(sums) == 1 ? sums[1] : reduce(*, sums)
        term = has_term ? term * s : (has_term = true; s)
    end
    if !isempty(rats)
        r = length(rats) == 1 ? rats[1] : reduce(*, rats)
        term = has_term ? term * r : (has_term = true; r)
    end

    # 3) reattach exp/log/pow factors (keep a stable order: exp, log, power)
    final_terms = CFunction[]
    if has_term
        push!(final_terms, term)
    end
    if length(exps) > 1
        exps = [reduce(*, exps)]   # you already define *(::CExp, ::CExp)
    end
    append!(final_terms, exps)
    append!(final_terms, logs)
    append!(final_terms, pows)     # carry powers as-is (optionally merge same-base here)

    # 4) trivial cases
    if isempty(final_terms)
        # only scalar coeff remained; keep dimensions from first original term
        return CAtom(param_info, coeff, zeros(Int, dims(terms[1])))
    elseif length(final_terms) == 1
        return coeff * final_terms[1]
    end

    # 5) canonical sort of the tail (do not reorder the leading combined term)
    start_idx = has_term ? 2 : 1
    if start_idx <= length(final_terms)-1
        sort!(final_terms[start_idx:end])
    end

    return CProd(param_info, coeff, final_terms, Val(:nosimp))
end

# Helper 
function collect_prod_terms(param_info::ParameterInfo, coeff::ComplexRational, terms::AbstractVector{<:CFunction})
    atoms, sums, rats, exps, logs, pows = CAtom[], CSum[], CRational[], CExp[], CLog[], CPower[]
    for t in terms
        t = simplify(t)
        if t isa CProd
            c2, a2, s2, r2, e2, l2, p2 = collect_prod_terms(param_info, t.coeff, t.terms)
            coeff *= c2
            append!(atoms, a2); append!(sums, s2)
            append!(rats,  r2); append!(exps, e2); append!(logs, l2); append!(pows, p2)
        elseif t isa CAtom
            coeff *= t.coeff
            push!(atoms, CAtom(param_info, ComplexRational(1,0,1), copy(t.var_exponents)))
        elseif t isa CSum
            push!(sums, t)
        elseif t isa CRational
            push!(rats, t)
        elseif t isa CExp
            push!(exps, t)
        elseif t isa CLog
            push!(logs, t)
        elseif t isa CPower
            push!(pows, t)
        else
            error("unsupported term in CProd: $t")
        end
    end
    return coeff, atoms, sums, rats, exps, logs, pows
end


######## Helper function(s) ###################################################################################################

firstnegative(a::CAtom) = is_negative(a.coeff)
firstnegative(a::CAbstract) = false
firstnegative(s::CSum)  = firstnegative(s.terms[1])
firstnegative(p::CMultiComposite) = is_negative(p.coeff)
function firstnegative(r::CRational) 
    if allnegative(r.denom) || (FLIP_IF_FIRST_TERM_NEGATIVE  && firstnegative(r.denom))   # prefer negatives on numerator
        r.numer = -r.numer
        r.denom = -r.denom
        return true
    end
    return false
end
firstnegative(x::CComposite)  = is_negative(x.coeff)

#### can be added?!
function addable(a::CFunction, b::CFunction)
    error("Unimplemented addable for $typeof(a) and $typeof(b).")
end
function addable(a::CAtom, b::CAtom)::Bool
    return a.var_exponents == b.var_exponents
end
function addable(a::CSum, b::CSum)::Bool
    true
end
function addable(a::CProd, b::CProd)::Bool
    if length(a.terms) != length(b.terms) 
        return false
    end
    for (el_a, el_b) in zip(a.terms, b.terms)
        if el_a != el_b
            return false
        end
    end
    return true
end
function addable(a::CRational, b::CRational)::Bool  
    return a.denom == b.denom
end
function addable(a::CExp, b::CExp)::Bool
    if a.x == b.x
        return true
    end 
    return false
end
function addable(a::CLog, b::CLog)::Bool
    if a.x == b.x
        return true
    end 
    return false
end


# assume addable
function unify_add(a::CAtom, b::CAtom)::CFunction
    absum = a.coeff+b.coeff
    return CAtom(absum, copy(a.var_exponents))
end
function unify_add(a::CSum, b::CSum)::CFunction
    error("CSum shouldn't contain another CSum!")
    #return a+b
end
function unify_add(a::CProd, b::CProd)::CFunction
    absum = a.coeff+b.coeff
    return CProd(absum, copy(a.terms), Val(:nosimp))
end
function unify_add(a::CRational, b::CRational)::CFunction
    simple_numer = simplify(a.numer+b.numer)
    return CRational(simple_numer, copy(a.denom), Val(:nosimp))
end
function unify_add(a::CExp, b::CExp)::CFunction
    absum = a.coeff+b.coeff
    return CExp(absum, copy(a.x), Val(:nosimp))
end
function unify_add(a::CLog, b::CLog)::CFunction
    absum = a.coeff+b.coeff
    return CLog(absum, copy(a.x), Val(:nosimp))
end 


# divisors 
function divisors(a::T) where T <: CFunction
    error("Not implemented for type $(typeof(a)).")
end
function divisors(a::CAtom)::Vector{Int}
    return [a.coeff.c]
end
function divisors(a::CSum)::Vector{Int}
    return reduce(vcat, [divisors(t) for t in a.terms])
end
function divisors(a::CProd)::Vector{Int}
    return [a.coeff.c]
end
function divisors(a::CRational)::Vector{Int}
    return reduce(vcat, [divisors(t) for t in [a.numer, a.denom]])
end
function divisors(a::CExp)::Vector{Int}
    return [a.coeff.c]
end
function divisors(a::CLog)::Vector{Int}
    return [a.coeff.c]
end

function vec_multiply(x::CAtom, vector::Vector{Int})::CAtom
    return CAtom(x.coeff, x.var_exponents + vector)
end
function vec_multiply(x::T, vector::Vector{Int})::T where T <: CComposite
    return modify_exprs([vec_multiply(t, vector) for t in x.terms], Val(:nosimp))
end
function vec_multiply(x::CRational, vector::Vector{Int})::CRational
    return CRational(vec_multiply(x.numer), vec_multiply(x.denom), Val(:nosimp))
end
function vec_multiply(x::CProd, vector::Vector{Int})::CProd
    terms = x.terms
    terms[1] = vec_multiply(terms[1], vector)   
    return CProd(x.coeff, terms, Va(:nosimp)) 
end


# Assumes that the divisors are 1 
function numer(a::T) where T <: CFunction
    error("Not implemented for type $(typeof(a)).")
end
function numer(s::CAtom)::Vector{Int}
    # if subtype complex take real and imaginary parts separately 
    return [s.coeff.a, s.coeff.b]
end
function numer(s::CSum)::Vector{Int}
    ints::Vector{Int} = Int[]
    for t in s.terms 
        append!(ints, numer(t))
    end
    return ints 
end
function numer(a::CProd)::Vector{Int}
    return [a.coeff.a, a.coeff.b]
end
function numer(s::CRational)::Vector{Int}
    return vcat(numer(s.numer), numer(s.denom))
end
function numer(a::CExp)::Vector{Int}
    return [a.coeff.a, a.coeff.b]
end
function numer(a::CLog)::Vector{Int}
    return [a.coeff.a, a.coeff.b]
end

import Base: gcd
function gcd(s::CFunction)
    return gcd(numer(s))   
end
function gcd(s::CFunction, t::CFunction)
    return gcd(vcat(numer(s), numer(t)))
end
#gcd(a*2+b+a*4)