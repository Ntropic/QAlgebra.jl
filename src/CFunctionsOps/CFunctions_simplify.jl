function simplify end

function simplify(s::CAtom)
    return s 
end

function simplify_CSum(elements::AbstractVector{<:CFunction})
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
    end
    return CSum(new_elements, Val{:nosimp}())
end
function simplify(s::CSum)
    # first simplify the lower levels 
    return simplify_CSum(s.terms)
end

function simplify_CRational(n::CFunction, d::CFunction)
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
        error("Dividing by zero")
    elseif isa(d, CAtom)
        if isnumeric(d)
            return n/d.coeff
        end
        return CRational(n, d, Val{:nosimp}())
    elseif isa(d, CRational)
        return (n*d.denom) / d.numer
    else
        return CRational(n, d, Val{:nosimp}())
    end
end
function simplify(r::CRational)
    return simplify_CRational(r.numer, r.denom)
end

function simplify_CExp(coeff::ComplexRational, x::CFunction)
    # exp(log(y)) ⇒ y
    if x isa CLog
        return x.x * coeff 
    # exp(0) ⇒ 1
    elseif iszero(x)
        return CAtom(coeff, zeros(Int, dims(x)))
    end
    return CExp(coeff, x, Val{:nosimp}())
end
function simplify(e::CExp)
    return simplify_CExp(e.coeff, e.x)
end

function simplify_CLog(coeff::ComplexRational, x::CFunction)
    # log(exp(y)) ⇒ y
    if x isa CExp
        return x.x * coeff 
    end
    # log(1) ⇒ 0
    if isone(x)
        return CAtom(0, zeros(Int, dims(x)))
    end
    return CLog(coeff, x, Val{:nosimp}()) 
end
function simplify(l::CLog)
    return simplify_CLog(l.coeff, l.x)
end

# CAtom, CSum, CProd, CRational, CExp, CLog

function simplify_CProd(coeff::ComplexRational, terms::AbstractVector{<:CFunction})
    #    simplify(p::CProd) -> CFunction
    # Normalize a product by
    #  • flattening nested CProd  
    #  • merging all CAtom → single coeff & one “pure” CAtom  
    #  • merging at most one CSum and one CRational via `*`  
    #  • sorting all remaining factors (commutative)  
    #  • pulling out common divisors into coeff
    coeff, atoms, sums, rats, exps, logs = collect_prod_terms(coeff, terms)   # also simplifies recursively! 

    # 1) combine atoms back into one
    has_term = false
    if !isempty(atoms) 
        var_ex = zeros(Int, length(atoms[1].var_exponents))
        for a in atoms
            var_ex .+= a.var_exponents
        end
        base_atom = CAtom(ComplexRational(1,0,1), var_ex)
        # 2) multiply sums/rationals (there’s at most one of each)
        term = base_atom
        has_term = true
    end
    if !isempty(sums)
        if has_term
            term = term * reduce(*, sums)
        else
            if length(sums)==1
                term = sums[1]
            else
                term = reduce(*, sums)
            end
            has_term = true
        end
    end
    if !isempty(rats)
        if has_term
            term = term * reduce(*, rats)
        else
            if length(rats)==1
                term = rats[1]
            else
                term = reduce(*, rats)
            end
            has_term = true
        end
    end

    # 3) re‑attach exp/log factors
    #    keep order: exps, logs
    final_terms = CFunction[]
    if has_term
        push!(final_terms, term)
    end
    if length(exps) > 1 
        exps = [reduce(*, exps)] # * is defined for CExp and CExp
    end
    append!(final_terms, exps)
    append!(final_terms, logs)

    # 5) sort the tail factors canonically
    sort!(final_terms[2:end])

    # 6) if only one term, drop the CProd wrapper
    if length(final_terms) == 1
        return coeff*final_terms[1]
    end
    return CProd(coeff, final_terms, Val{:nosimp}())
end

function simplify(p::CProd)
    return simplify_CProd(p.coeff, p.terms)
end

# Helper 
function collect_prod_terms(coeff::ComplexRational, terms::AbstractVector{<:CFunction})
    atoms, sums, rats, exps, logs = CAtom[], CSum[], CRational[], CExp[], CLog[]
    for t in terms
        t = simplify(t)
        if t isa CProd
            c2, a2,s2,r2,e2,l2 = collect_prod_terms(t.coeff, t.terms)
            coeff *= c2
            append!(atoms, a2); append!(sums, s2)
            append!(rats,  r2); append!(exps, e2); append!(logs, l2)
        elseif t isa CAtom
            coeff *= t.coeff
            push!(atoms, CAtom(1, copy(t.var_exponents)))
        elseif t isa CSum
            push!(sums, t)
        elseif t isa CRational
            push!(rats, t)
        elseif t isa CExp
            push!(exps, t)
        elseif t isa CLog
            push!(logs, t)
        else
            error("unsupported term in CProd: $t")
        end
    end
    return coeff, atoms, sums, rats, exps, logs
end


######## Helper function(s) ###################################################################################################

firstnegative(a::CAtom) = is_negative(a.coeff)
firstnegative(s::CSum)  = firstnegative(s.terms[1])
firstnegative(p::CProd) = is_negative(p.coeff)
function firstnegative(r::CRational) 
    if allnegative(r.denom) || (FLIP_IF_FIRST_TERM_NEGATIVE  && firstnegative(r.denom))   # prefer negatives on numerator
        r.numer = -r.numer
        r.denom = -r.denom
        return true
    end
    return false
end
firstnegative(x::CExp)  = is_negative(x.coeff)
firstnegative(x::CLog)  = is_negative(x.coeff)

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
    return CProd(absum, copy(a.terms), Val{:nosimp}())
end
function unify_add(a::CRational, b::CRational)::CFunction
    simple_numer = simplify(a.numer+b.numer)
    return CRational(simple_numer, copy(a.denom), Val{:nosimp}())
end
function unify_add(a::CExp, b::CExp)::CFunction
    absum = a.coeff+b.coeff
    return CExp(absum, copy(a.x), Val{:nosimp}())
end
function unify_add(a::CLog, b::CLog)::CFunction
    absum = a.coeff+b.coeff
    return CLog(absum, copy(a.x), Val{:nosimp}())
end 

issimple(f::CFunction) = error("Not implemented for type $(typeof(f)).")
issimple(f::CSum) = all(t -> isa(t, CAtom), f.terms)
issimple(f::CProd) = all(t -> issimple(t), f.terms)
issimple(f::CRational) = issimple(f.numer) && issimple(f.denom)
issimple(f::CAtom) = true
issimple(f::CExp) = false 
issimple(f::CLog) = false 


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
function vec_multiply(x::CSum, vector::Vector{Int})::CSum
    return CSum([vec_multiply(t, vector) for t in x.terms], Val{:nosimp}())
end
function vec_multiply(x::CRational, vector::Vector{Int})::CRational
    return CRational(vec_multiply(x.numer), vec_multiply(x.denom), Val{:nosimp}())
end
function vec_multiply(x::CProd, vector::Vector{Int})::CProd
    terms = x.terms
    terms[1] = vec_multiply(terms[1], vector)   
    return CProd(x.coeff, terms, Val{:nosimp}()) 
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
    return vcat([numer(t) for t in s.terms]...)
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