"""
    tree_iter(f::CFunction)

Iterator that yields a node `f` and then all of its children,
recursively (depth-first, pre-order).
"""
tree_iter(a::CAtom)     = (a,)                # leaf
tree_iter(a::CAbstract) = (a,)                # leaf
tree_iter(f::CComposite) = Iterators.flatten(((f,), tree_iter(f.expr)))
tree_iter(f::CMultiComposite) = Iterators.flatten(((f,), (tree_iter(ch) for ch in f.expr)))
tree_iter(r::CRational) = Iterators.flatten(((r,), tree_iter(r.numer), tree_iter(r.denom)))
tree_iter(M::CMatrix) = Iterators.flatten(((M,), (tree_iter(ch) for ch in M.expr[:])))


var_exponents_iter(a::CAtom) = (a.var_exponents,)  # 1-tuple, no alloc beyond the tuple
var_exponents_iter(q::CAbstract) = ((zeros(Int, dims(q))),)  # if you want to treat them as zero-terms
var_exponents_iter(f::CComposite) = var_exponents_iter(f.expr)
var_exponents_iter(f::CMultiComposite) = Iterators.flatten(var_exponents_iter.(f.expr))
var_exponents_iter(r::CRational) = Iterators.flatten((var_exponents_iter(r.numer), var_exponents_iter(r.denom)))
var_exponents_iter(q::CCustomType) = Iterators.flatten((Iterators.flatten(var_exponents_iter.(q.expr)),var_exponents_iter(q.ctype_def.fun)))
var_exponents_iter(v::CVector) = Iterators.flatten(var_exponents_iter.(v.expr))
var_exponents_iter(M::CMatrix) = Iterators.flatten(var_exponents_iter.(M.expr[:]))

"""
    var_exponents_iter_simple(f::CFunction)

Like `var_exponents_iter`, but only includes exponents from
`CAtom`, `CSum`, and `CProd`.
Other function types contribute nothing.
"""
var_exponents_iter_simple(a::CAtom) = (a.var_exponents,)
var_exponents_iter_simple(s::CSum)  = Iterators.flatten(var_exponents_iter_simple.(s.expr))
var_exponents_iter_simple(p::CProd) = Iterators.flatten(var_exponents_iter_simple.(p.expr))
var_exponents_iter_simple(f::CFunction) = (zeros(Int, dims(f)),)

"""
    isnumeric(f::CFunction)

True if all exponents in all terms are zero.
"""
isnumeric(f::CFunction) = all(exps -> all(==(0), exps), var_exponents_iter(f))

import Base: iszero, isempty, isone

"""
    iszero(a::CAtom)     -> Bool
    iszero(s::CSum)      -> Bool
    iszero(r::CRational) -> Bool

Returns `true` if the expression is identically zero:
"""
iszero(a::CAtom)        = iszero(a.coeff)
iszero(s::CSum)         = all(iszero, s.expr)
iszero(p::CComposite)        = iszero(p.coeff) || iszero(p.expr)
iszero(p::CMultiComposite)        = iszero(p.coeff) || any(iszero, p.expr)
iszero(r::CRational)    = iszero(r.numer)
function iszero(c::CCustomType)
    # 1. First check coeff
    iszero(c.coeff) && return true
    # 2. Substitute arguments into the base definition
    substituted = substitute_and_simplify(c.ctype_def.fun, c.ctype_def.abstract_parameters, c.expr)
    # 3. Check if the expanded form is zero
    return iszero(substituted)
end

isempty(s::CMultiComposite)        = isempty(s.expr)

isone(c::CFunction) = false 
isone(a::CAtom)     = isnumeric(a) && isone(a.coeff)

allnegative(a::CAtom) = is_negative(a.coeff)
allnegative(s::CSum)  = !isempty(s.expr) && all(allnegative, s.expr)
allnegative(p::CProd) = is_negative(p.coeff)
allnegative(r::CRational) = allnegative(r.numer)
allnegative(x::CExp)  = false
allnegative(x::CLog)  = false

"""
    min_exponents(f::CFunction) -> Vector{Int}

Component-wise minimum of all monomial exponent vectors appearing in `f`.
If `f` contains no atoms (e.g. empty containers), returns `Int[]`.

Relies on `var_exponents_iter(::CFunction)`.
"""
function min_exponents(f::CFunction)::Vector{Int}
    n = dims(f)                    # number of polynomial variables
    mins = zeros(Int, n)           # start with all zeros
    for exps in var_exponents_iter_simple(f)
        mins = min.(mins, exps)
    end
    return mins
end

"""
    contains_c_indexes(f::CFunction, idxs::Vector{Int})

True if any exponent at positions `idxs` is nonzero.
"""
function contains_c_indexes(f::CFunction, idxs::Vector{Int})
    for exps in var_exponents_iter(f)
        @inbounds for i in idxs
            if exps[i] != 0
                return true
            end
        end
    end
    return false
end


# Helper: does this expression contain any CVector or CMatrix?
contains_vec_or_mat(::CFunction) = false
contains_vec_or_mat(q::CCustomType) = any(contains_vec_or_mat, q.expr) || contains_vec_or_mat(q.ctype_def.fun)
contains_vec_or_mat(e::CComposite) = contains_vec_or_mat(e.expr)
contains_vec_or_mat(s::CMultiComposite) = any(contains_vec_or_mat, s.expr)
contains_vec_or_mat(r::CRational) = contains_vec_or_mat(r.numer) || contains_vec_or_mat(r.denom)
contains_vec_or_mat(::CVector) = true
contains_vec_or_mat(::CMatrix) = true

# --- numeric printing gate (specify if the coefficient is needed) -----------------------------------------------------
"""
    printnumeric(f::CFunction) -> Bool

Should `f` be printed in numeric form?
Default: true.
"""
printnumeric(::CFunction) = true

function printnumeric(f::CAtom)::Bool
    # Don’t print if numeric with coeff == 1
    isnumeric(f) && isonelike(f.coeff) ? false : true
end
function printnumeric(f::CMultiComposite)::Bool
    # If it’s a singleton, delegate to the single child
    length(f.expr) == 1 ? printnumeric(f.expr[1]) : true
end
printnumeric(::CRational) = true  # Rational: always print
printnumeric(v::CComposite) = !(isnumeric(v) && isonelike(v.coeff))



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

function get_coeffs_ops(f::CSum)::Tuple{Vector{ComplexRational}, Vector{Vector{Int}}}
    return vcat(coeff.(f.terms)...), vcat(var_exponents.(f.terms)...)
    #return [coeff(term) for term in f.terms], [term.var_exponents for term in f.terms]
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

function isonelike(f::CAtom)::Bool
    return isonelike(f.coeff) && isnumeric(f)
end
function simple_CSum(f::CSum)::Bool
    return !any(term -> typeof(term)==CAtom, f.terms)
end
function simple_CSum(f::CAtom)::Bool
    return true 
end
function simple_CSum(f::CFunction)::Bool
    return false 
end

function separate_CSum(f::CSum )::Tuple{Bool, Union{CAtom, Nothing}, CSum}
    if simple_CSum(f) 
        coeffs, vs = get_coeffs_ops(f)
        if has_common_multiple(coeffs) 
            base, multiples = common_denominator_form(coeffs)
            # get exponent offsets 
            offset = common_exponent_offset(vs)
            new_f = _CSum([CAtom(m, v.-offset) for (m, v) in zip(multiples, vs)])
            pre_f = CAtom(base, offset)
            return true, pre_f, new_f 
        end
    end
    return false, nothing, f
end

function simple_combinable_F(t1::T1, t2::T2)::Tuple{Bool, ComplexRational} where {T1 <: CFunction, T2 <: CFunction}
    # Transform CAtom's to CSums 
    if !isa(t1, CSum)
        t1 = [t1]
    else 
        t1 = t1.terms
    end
    if !isa(t2, CSum)
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
    
    # check if corresponding terms have the similar exponents (i.e. the difference of exponents has to be the same for each term in CSum)
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
# check for multiple elements in a Vector of CSum 
function simple_combinable_Fs(ts::AbstractVector{<:CFunction})::Tuple{Vector{Vector{<:CFunction}}, Vector{Vector{Int}}}
    groups = Vector{Vector{Union{CAtom,CSum}}}()
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

function ratios_Fs(ts::AbstractVector{Union{CAtom,CSum}})::Vector{ComplexRational}
    if ts[1] isa CAtom 
        c1 = ts[1].coeff
    else 
        c1 = ts[1].terms[1].coeff
    end
    ratios = ComplexRational[ComplexRational(1,0,1)]
    for t in ts[2:end] 
        if t isa CAtom 
            c2 = t.coeff
        else 
            c2 = t.terms[1].coeff
        end
        push!(ratios, c1/c2)
    end
    return ratios
end
function group_Fs(ts::AbstractVector{Union{CAtom,CSum}})::Union{CAtom, CSum, Tuple{Union{CAtom, CSum}, AbstractVector{Union{CAtom, CSum}}}}
    # find the correct way to group a group of Fs (they must be groupable, create the input vector with simple_combinable_Fs)
    if length(ts) == 1
        return ts[1]
    elseif length(ts) > 1
        # we need to find the best common ratio 
        ratios = ratios_Fs(ts)
        base, multiples = common_denominator_form(ratios)
        pre_F = ts[1]*base 

        if isa(pre_F, CAtom)
            curr_var_exponents = zeros(Int, length(ts[1].var_exponents))
        else
            curr_var_exponents = zeros(Int, length(ts[1].terms[1].var_exponents))
        end
        post_Fs = Union{CAtom, CSum}[ CAtom(m, curr_var_exponents) for m in multiples ]
        return (pre_F, post_Fs)
    end
end

"""
    how_to_combine_Fs(ts::Vector{Union{CAtom, CSum}}) :: Tuple{Vector{Union{CAtom, CSum, Tuple{CAtom, Vector{Union{CAtom, CSum}}}}}, Vector{Vector{Int}}}

Groups and combines `CAtom` and `CSum` objects in the input vector `ts` into composite structures that can be processed together. 
Returns a tuple containing both the groups as (CFunction) elements of a Vector and the indexs corresponding to the elements in the groups. 
The (CFunction) grouping is given either by a single CFunction element (either a single `CAtom` or a `CSum`) or by a tuple of an `CAtom` (F1) containing the shared factors, and a vector of CFunction elements (F2_i), so that 
together they represent a term of the form: F1 * (F2_1 + ... + F2_n). 
"""
function how_to_combine_Fs(ts::Vector{Union{CAtom,CSum}}) #::Tuple{Vector{Union{CAtom,CSum, Tuple{CAtom, Vector{Union{CAtom, CSum}}}}}, Vector{Vector{Int}}}
    if length(ts) == 1
        return [ts[1]], [[1]]
    elseif length(ts) > 1 
        groups, indexes = simple_combinable_Fs(ts)
        return [group_Fs(grp) for grp in groups], indexes
    else
        return [], []
    end
end

# --- make two indexes equal (variable reindexing/merging) ---------------------

# Moves exponent from j -> i for each pair (i,j) in coeff_ind_order.
# Returns (changed_any::Bool, transformed_expression)

function term_equal_indexes(atom::CAtom, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CAtom}
    new_exponents = copy(atom.var_exponents)
    changed_any = false
    @inbounds for (i, j) in coeff_ind_order
        ei = new_exponents[j]
        changed_any |= (ei != 0)
        new_exponents[i] += ei
        new_exponents[j] = 0
    end
    return changed_any, CAtom(atom.coeff, new_exponents)
end

function term_equal_indexes(fsum::CSum, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CSum}
    changed_any = false
    new_terms = Vector{CFunction}(undef, length(fsum.terms))
    @inbounds for k in eachindex(fsum.terms)
        changed, new_term = term_equal_indexes(fsum.terms[k], coeff_ind_order)
        changed_any |= changed
        new_terms[k] = new_term
    end
    return changed_any, _CSum(new_terms)
end

function term_equal_indexes(frational::CRational, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CRational}
    changed_num, new_num = term_equal_indexes(frational.numer, coeff_ind_order)
    changed_den, new_den = term_equal_indexes(frational.denom, coeff_ind_order)
    return (changed_num || changed_den), CRational(new_num, new_den, Val(:nosimp))
end

# NEW
function term_equal_indexes(fprod::CProd, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CFunction}
    changed_any = false
    new_terms = Vector{CFunction}(undef, length(fprod.terms))
    @inbounds for k in eachindex(fprod.terms)
        changed, new_term = term_equal_indexes(fprod.terms[k], coeff_ind_order)
        changed_any |= changed
        new_terms[k] = new_term
    end
    return changed_any, CProd(fprod.coeff, new_terms, Val(:nosimp))
end

function term_equal_indexes(fexp::CExp, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CFunction}
    changed, nx = term_equal_indexes(fexp.x, coeff_ind_order)
    return changed, CExp(fexp.coeff, nx, Val(:nosimp))
end

function term_equal_indexes(flog::CLog, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CFunction}
    changed, nx = term_equal_indexes(flog.x, coeff_ind_order)
    return changed, CLog(flog.coeff, nx, Val(:nosimp))
end

function term_equal_indexes(fpwr::CPower, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CFunction}
    changed, nx = term_equal_indexes(fpwr.x, coeff_ind_order)
    return changed, CPower(fpwr.coeff, nx, fpwr.exponent, Val(:nosimp))
end

function term_equal_indexes(v::CVector, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CVector}
    changed_any = false
    new_entries = Vector{CFunction}(undef, length(v.expr))
    @inbounds for k in eachindex(v.expr)
        changed, e = term_equal_indexes(v.expr[k], coeff_ind_order)
        changed_any |= changed
        new_entries[k] = e
    end
    return changed_any, CVector(v.coeff, new_entries; row=v.row)
end

function term_equal_indexes(M::CMatrix, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CMatrix}
    changed_any = false
    flat = M.expr[:]
    new_flat = Vector{CFunction}(undef, length(flat))
    @inbounds for k in eachindex(flat)
        changed, e = term_equal_indexes(flat[k], coeff_ind_order)
        changed_any |= changed
        new_flat[k] = e
    end
    new_mat = reshape(new_flat, size(M.expr))
    return changed_any, CMatrix(M.coeff, new_mat)
end