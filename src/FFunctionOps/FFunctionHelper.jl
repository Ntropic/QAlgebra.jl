function isonelike(c::ComplexRational)::Bool 
    if c.c == 1 && ( c.b == 0 && abs(c.a) == 1 )
        return true 
    end
    return false
end
function printnumeric(f::FAtom)::Bool
    if isnumeric(f) && isonelike(f.coeff) 
        return false 
    end
    return true 
end
function printnumeric(f::FSum)::Bool
    if length(f.terms) == 1 
        return printnumeric(f.terms[1])
    end 
    return true 
end
function printnumeric(f::FRational)::Bool
    return true 
end

function is_axis_multiple(x::ComplexRational, y::ComplexRational)::Tuple{Bool, ComplexRational}
    ratio = x / y   # y * ratio = x 
    # is the ratio either on.y real or only on the imaginary axis? 
    if ratio.b == 0 || ratio.a == 0
        return true, ratio
    else
        return false, ratio
    end
end

function has_common_multiple(v::Vector{ComplexRational})::Bool
    for i in 2:length(v) 
        if !is_axis_multiple(v[1], v[i])[1]'
            return false
        end
    end
    return true
end

function is_axis_multiple(v1::Vector{ComplexRational}, v2::Vector{ComplexRational})::Tuple{Bool, ComplexRational}
    @assert length(v1) == length(v2) "Vectors must have same length"

    # Compute the ratio for the first pair
    r = v1[1] / v2[1]

    for i in 2:length(v1)
        ri = v1[i] / v2[i]
        if ri.a * r.c != r.a * ri.c || ri.b * r.c != r.b * ri.c
            return false, r
        end
    end

    return true, r
end

function common_denominator_form(v::Vector{ComplexRational})::Tuple{ComplexRational, Vector{ComplexRational}}
    @assert length(v) ≥ 1 "You must supply at least one ComplexRational."
    c1 = v[1]   # ComplexRational(a1, b1, c1den)

    ratios = [ v[i] / c1  for i in 1:length(v) ]

    for (i, r) in enumerate(ratios)
        if !(r.a == 0 || r.b == 0)
            # If both a≠0 and b≠0, then v[i] is not a pure real‐ or pure imag‐multiple of c1.
            error("Entry #$(i) = $(v[i]) is not axis‐aligned with v[1] = $(c1).")
        end
    end

    dens = [ r.c for r in ratios ] 
    D    = foldl(lcm, dens)     # positive Int

    base = ComplexRational(c1.a, c1.b, c1.c * D)

    multiples = ComplexRational[]
    for i in 1:length(v)
        push!(multiples, v[i] / base)
    end
    return base, multiples
end

function get_coeffs_ops(f::FSum)::Tuple{Vector{ComplexRational}, Vector{Vector{Int}}}
    return [term.coeff for term in f.terms], [term.var_exponents for term in f.terms]
end

function common_exponent_offset(exponents::Vector{Vector{Int}})::Vector{Int}
    @assert !isempty(exponents)
    n = length(exponents[1])
    @assert all(length(e) == n for e in exponents)

    offset = Int[]
    for i in 1:n
        vals = [e[i] for e in exponents]
        if all(x -> x > 0, vals)
            push!(offset, minimum(vals))
        elseif all(x -> x < 0, vals)
            push!(offset, maximum(vals))  # smallest (least negative)
        else
            push!(offset, 0)  # neutral element 
        end
    end
    return offset
end

function isonelike(f::FAtom)::Bool
    return isonelike(f.coeff)
end
function simple_FSum(f::FSum) 
    return all([typeof(term)==FAtom for term in f.terms])
end
function separate_FSum(f::FSum )::Tuple{Bool, Union{FAtom, Nothing}, FSum}
    if simple_FSum(f) 
        coeffs, vs = get_coeffs_ops(f)
        if has_common_multiple(coeffs) 
            base, multiples = common_denominator_form(coeffs)
            # get exponent offsets 
            offset = common_exponent_offset(vs)
            new_f = FSum([FAtom(m, v.-offset) for (m, v) in zip(multiples, vs)])
            pre_f = FAtom(base, offset)
            return true, pre_f, new_f 
        end
    end
    return false, nothing, copy(f)
end