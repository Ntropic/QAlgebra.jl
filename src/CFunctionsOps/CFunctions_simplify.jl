function simplify end

function simplify(s::CAtom)
    return s 
end

function simplify(e::CExp)
    y = simplify(e.x)
    # exp(log(y)) ⇒ y
    if y isa CLog
        return y.x
    # exp(0) ⇒ 1
    elseif iszero(simplify(y))
        return CAtom(1, zeros(Int, dims(y)))
    else
        return CExp(y)
    end
end

function simplify(l::CLog)
    y = simplify(l.x)
    # log(exp(y)) ⇒ y
    if y isa CExp
        return y
    end
    # log(1) ⇒ 0
    if isone(y)
        return CAtom(0, zeros(Int, dims(y)))
    end
    return CLog(y)
end

function simplify(r::CRational)
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
        error("Divinding by zero")
    elseif isa(d, CAtom)
        return n / d 
    elseif isa(d, CRational)
        return (n*d.denom) / d.numer
    else
        return CRational(n, d)
    end
end

function simplify(s::CSum)
    # first simplify the lower levels 
    elements = [simplify(e) for e in s.terms]
    s = CSum(elements)
    recursive_sort!(s)
    i = 1
    elements = s.terms
    new_elements = CFunction[]
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
        push!(new_elements, CAtom(0, zeros(Int, dims(s))))
    end
    return CSum(new_elements)
end



#    simplify(p::CProd) -> CFunction
# Normalize a product by
#  • flattening nested CProd  
#  • merging all CAtom → single coeff & one “pure” CAtom  
#  • merging at most one CSum and one CRational via `*`  
#  • sorting all remaining factors (commutative)  
#  • pulling out common divisors into coeff
function simplify(p::CProd)

    coeff, atoms, sums, rats, exps, logs = collect_prod_terms(p)   # also simplifies recursively! 

    # 1) combine atoms back into one
    has_term = false
    if !isempty(atoms) 
        var_ex = zeros(Int, length(atoms[1].var_exponents))
        for a in atoms
            var_ex .+= a.var_exponents
        end
        base_atom = CAtom(coeff, var_ex)
        # 2) multiply sums/rationals (there’s at most one of each)
        term = base_atom
        has_term = true
    end
    if !isempty(sums)
        if has_term
            term = simplify(term * reduce(*, sums))
        else
            if length(sums)==1
                term = sums
            else
                term = simplify(reduce(*, sums))
            end
            has_term = true
        end
    end
    if !isempty(rats)
        if has_term
            term = simplify(term * reduce(*, rats))
        else
            if length(rats)==1
                term = rats
            else
                term = simplify(reduce(*, rats))
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
    append!(final_terms, exps)
    append!(final_terms, logs)

    # 5) sort the tail factors canonically
    sort!(final_terms[2:end])

    # 6) if only one term, drop the CProd wrapper
    return length(final_terms)==1 ? final_terms[1] : CProd(final_terms...)
end
# Helper 
function collect_prod_terms(p::CProd)
    coeff = ComplexRational(1,0,1)
    atoms, sums, rats, exps, logs = CAtom[], CSum[], CRational[], CExp[], CLog[]
    for t in p.terms
        t = simplify(t)
        if t isa CProd
            c2, a2,s2,r2,e2,l2 = collect_prod_terms(t)
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


function unifiable(a::CAtom, b::CAtom)::Bool
    return a.var_exponents == b.var_exponents
end
function unifiable(a::CRational, b::CRational)::Bool
    return a.denom == b.denom
end
function unifiable(a::CSum, b::CSum)::Bool
    true
end
# assume inifiable
function unify(a::CAtom, b::CAtom)::CFunction
    absum = a.coeff+b.coeff
    if iszero(absum)
        var_zeros = zeros(Int, dims(a))
        return CAtom(0, var_zeros)
    end
    return CAtom(absum, copy(a.var_exponents))
end
function unify(a::CRational, b::CRational)::CFunction
    simple_numer = simplify(a.numer+b.numer)
    if iszero(simple_numer)
        var_zeros = zeros(Int, dims) 
        return CAtom(0, var_zeros)
    end
    return CRational(simple_numer, a.denom)
end
function unify(a::CSum, b::CSum)::CFunction
    return simplify(a+b)
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
    return CSum([vec_multiply(t, vector) for t in x.terms])
end
function vec_multiply(x::CRational, vector::Vector{Int})::CRational
    return CRational(vec_multiply(x.numer), vec_multiply(x.denom))
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
    return [s.coeff.a, s.coeff.b]
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