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
function simple_FSum(f::FSum)::Bool
    return all([typeof(term)==FAtom for term in f.terms])
end
function simple_FSum(f::FAtom)::Bool
    return true 
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

function simple_combinable_F(t1::Union{FAtom, FSum}, t2::Union{FAtom, FSum})::Tuple{Bool, ComplexRational}
    # Transform FAtom's to FSums 
    if t1 isa FAtom
        t1 = [t1]
    else 
        t1 = t1.terms
    end
    if t2 isa FAtom
        t2 = [t2]
    else
        t2 = t2.terms
    end

    #@assert length(t1) == length(t2) "Cannot separate pair terms with different lengths."
    if length(t1) != length(t2)
        return false, ComplexRational(0, 0, 1)
    end

    if all(a -> all(e == 0 for e in a.var_exponents), t1)
        return false, ComplexRational(0, 0, 1)
    end
    
    # check if corresponding terms have the similar exponents (i.e. the difference of exponents has to be the same for each term in FSum)
    exponents_diff = t1[1].var_exponents .- t2[1].var_exponents
    if !all([diff == 0 for diff in exponents_diff])
        return false, ComplexRational(0, 0, 1)
    end
    for (a1, a2) in zip(t1[2:end], t2[2:end])
        if a1.var_exponents .- a2.var_exponents != exponents_diff
            return false, ComplexRational(0, 0, 1)
        end
    end
    
    # check if the terms have a constant ratio in the coefficients 
    ratio = t1[1].coeff / t2[1].coeff
    if !(ratio.b == 0 || ratio.a == 0)
        return false, ratio
    end
    for (a1, a2) in zip(t1[2:end], t2[2:end])
        if a1.coeff / a2.coeff != ratio
            return false, ratio
        end
    end
    return true, ratio
end
# check for multiple elements in a Vector of FSum 
function simple_combinable_Fs(ts::AbstractVector{Union{FAtom,FSum}})::Tuple{Vector{Vector{Union{FAtom, FSum}}}, Vector{Vector{Int}}}
    groups = Vector{Vector{Union{FAtom,FSum}}}()
    indexes = Vector{Vector{Int}}()
    for (i, t) in enumerate(ts)
        placed = false
        for (inds, grp) in zip(indexes, groups)
            ok,_ = simple_combinable_F(grp[1], t)
            if ok 
                push!(grp, t)
                push!(inds, i)
                placed = true
                break
            end
        end
        if !placed
            push!(groups, [t])
            push!(indexes, [i])
        end
    end
    return groups, indexes
end

function ratios_Fs(ts::AbstractVector{Union{FAtom,FSum}})::Vector{ComplexRational}
    if ts[1] isa FAtom 
        c1 = ts[1].coeff
    else 
        c1 = ts[1].terms[1].coeff
    end
    ratios = ComplexRational[ComplexRational(1,0,1)]
    for t in ts[2:end] 
        if t isa FAtom 
            c2 = t.coeff
        else 
            c2 = t.terms[1].coeff
        end
        push!(ratios, c1/c2)
    end
    return ratios
end
function group_Fs(ts::AbstractVector{Union{FAtom,FSum}})::Union{FAtom, FSum, Tuple{Union{FAtom, FSum}, AbstractVector{Union{FAtom, FSum}}}}
    # find the correct way to group a group of Fs (they must be groupable, create the input vector with simple_combinable_Fs)
    if length(ts) == 1
        return ts[1]
    elseif length(ts) > 1
        # we need to find the best common ratio 
        ratios = ratios_Fs(ts)
        base, multiples = common_denominator_form(ratios)
        pre_F = ts[1]*base 

        if isa(pre_F, FAtom)
            curr_var_exponents = zeros(Int, length(ts[1].var_exponents))
        else
            curr_var_exponents = zeros(Int, length(ts[1].terms[1].var_exponents))
        end
        post_Fs = Union{FAtom, FSum}[ FAtom(m, curr_var_exponents) for m in multiples ]
        return (pre_F, post_Fs)
    end
end

"""
    how_to_combine_Fs(ts::Vector{Union{FAtom, FSum}}) :: Tuple{Vector{Union{FAtom, FSum, Tuple{FAtom, Vector{Union{FAtom, FSum}}}}}, Vector{Vector{Int}}}

Groups and combines `FAtom` and `FSum` objects in the input vector `ts` into composite structures that can be processed together. 
Returns a tuple containing both the groups as (FFunction) elements of a Vector and the indexs corresponding to the elements in the groups. 
The (FFunction) grouping is given either by a single FFunction element (either a single `FAtom` or a `FSum`) or by a tuple of an `FAtom` (F1) containing the shared factors, and a vector of FFunction elements (F2_i), so that 
together they represent a term of the form: F1 * (F2_1 + ... + F2_n). 
"""
function how_to_combine_Fs(ts::Vector{Union{FAtom,FSum}}) #::Tuple{Vector{Union{FAtom,FSum, Tuple{FAtom, Vector{Union{FAtom, FSum}}}}}, Vector{Vector{Int}}}
    if length(ts) == 1
        return [ts[1]], [[1]]
    elseif length(ts) > 1 
        groups, indexes = simple_combinable_Fs(ts)
        return [group_Fs(grp) for grp in groups], indexes
    else
        return [], []
    end
end


# Make 2 indexes equal => analogous to equally named function for qTerms in qExpressions.jl
function term_equal_indexes(atom::FAtom, coeff_inds1::Vector{Int}, coeff_inds2::Vector{Int})::Tuple{Bool,FAtom}
    changed_any::Bool = false
    new_exponents = copy(atom.var_exponents)
    for (i, j) in zip(coeff_inds1, coeff_inds2)
        if new_exponents[i] != 0
            changed_any = true
            new_exponents[j] += new_exponents[i]
            new_exponents[i] = 0
        end
    end
    return changed_any, FAtom(atom.coeff, new_exponents)
end
function term_equal_indexes(fsum::FSum, inds1::Vector{Int}, inds2::Vector{Int})::Tuple{Bool,FSum}
    changed_any::Bool = false
    new_terms::Vector{FFunction} = []
    for t in fsum.terms
        changed, new_term = term_equal_indexes(t, inds1, inds2)
        changed_any = changed_any || changed
        push!(new_terms, new_term)
    end
    return changed_any, FSum(fsum.index, new_terms)
end
function term_equal_indexes(frational::FRational, inds1::Vector{Int}, inds2::Vector{Int})::Tuple{Bool,FRational}
    changed_num, new_num = term_equal_indexes(frational.num, inds1, inds2)
    changed_den, new_den = term_equal_indexes(frational.den, inds1, inds2)
    return changed_num || changed_den, FRational(new_num, new_den)
end
