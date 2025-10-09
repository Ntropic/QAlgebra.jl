export separate_by_cond

"""
    separate_by_cond(c::CFunction, group1::BitVector, group2::BitVector)
        -> Vector{Tuple{CFunction,CFunction}}

Attempt to split the coefficient function `c` into factors whose parameters lie
inside two disjoint groups. Entries with `group1[i] == true` are required to be
separated from those with `group2[i] == true`. Parameters that are false in both
masks may appear in either output, but are preferentially assigned to the second
group.

Returns a vector of pairs `(left, right)` such that summing their products
reconstructs `c`. For every pair:
- `left` only depends on parameters from `group1`.
- `right` depends on parameters from `group2` or those unassigned in either
  mask.
- `left` always has coefficient `CR_ONE`; any scalar factors are carried by
  `right` (including the neutral `c_one` term when appropriate).

Throws an error if `c` still contains `CAbstract` leaves or the expression
cannot be separated without mixing both parameter sets (e.g. a power of a mixed
sum).
"""
function separate_by_cond(c::T, group1::BitVector, group2::BitVector) where T <: CFunction
    dims = c.param_info.dims
    length(group1) == dims || throw(DimensionMismatch("group1 mask must match parameter dimension"))
    length(group2) == dims || throw(DimensionMismatch("group2 mask must match parameter dimension"))
    any(group1 .& group2) && throw(ArgumentError("group masks must be disjoint"))
    flex = .!(group1 .| group2)
    _ensure_no_abstracts(c)

    pairs = _separate_pairs(c, group1, group2, flex)
    true_terms = CFunction[]
    false_terms = CFunction[]
    mixed_terms = Tuple{CFunction,CFunction}[]

    for (tp, fp) in pairs
        if iszero(tp) || iszero(fp)
            continue
        elseif isone(fp)
            push!(true_terms, tp)
        elseif isone(tp)
            push!(false_terms, fp)
        else
            push!(mixed_terms, (tp, fp))
        end
    end

    true_part = _combine_sum(c.param_info, true_terms)
    false_part = _combine_sum(c.param_info, false_terms)
    results = Tuple{CFunction,CFunction}[]
    if !iszero(true_part)
        push!(results, _normalize_pair(true_part, _one_atom(c.param_info)))
    end
    if !iszero(false_part)
        push!(results, _normalize_pair(_one_atom(c.param_info), false_part))
    end
    for pair in mixed_terms
        push!(results, _normalize_pair(pair[1], pair[2]))
    end
    return results
end

function separate_by_cond(c::T, cond::BitVector)  where T <: CFunction
    length(cond) == c.param_info.dims || throw(DimensionMismatch("condition mask must match parameter dimension"))
    return separate_by_cond(c, cond, .!cond)
end

# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

function _ensure_no_abstracts(f::CFunction)
    for leaf in leaf_iter(f)
        leaf isa CAbstract && throw(ArgumentError("separate_by_cond requires expressions without CAbstract leaves"))
    end
end

function _combine_sum(param_info::ParameterInfo, terms::Vector{CFunction})
    filtered = CFunction[]
    for term in terms
        iszero(term) && continue
        push!(filtered, term)
    end
    if isempty(filtered)
        return zero_catom(param_info)
    elseif length(filtered) == 1
        return filtered[1]
    else
        return _CSum(param_info, filtered)
    end
end

function _separate_pairs(f::CFunction, group1::BitVector, group2::BitVector, flex::BitVector)
    acting = where_acting(f)
    has_g1 = any(acting .& group1)
    has_other = any(acting .& (group2 .| flex))

    if has_g1 && has_other
        return _split_mixed(f, group1, group2, flex)
    elseif has_g1
        return [(f, _one_atom(f.param_info))]
    else
        return [(_one_atom(f.param_info), f)]
    end
end

function _split_mixed(a::CAtom, group1::BitVector, group2::BitVector, flex::BitVector)
    dims = a.param_info.dims
    true_exps = spzeros(Int, dims)
    false_exps = spzeros(Int, dims)
    nz_ind = a.var_exponents.nzind
    nz_val = a.var_exponents.nzval
    has_true = false
    has_false = false
    for (pos, idx) in enumerate(nz_ind)
        val = nz_val[pos]
        if group1[idx]
            true_exps[idx] = val
            has_true = true
        else
            false_exps[idx] = val
            has_false = true
        end
    end
    if has_true && has_false
        true_atom = CAtom(a.param_info, a.coeff, true_exps)
        false_atom = CAtom(a.param_info, CR_ONE, false_exps)
        return [(true_atom, false_atom)]
    elseif has_true
        return [(a, _one_atom(a.param_info))]
    elseif has_false
        false_atom = CAtom(a.param_info, a.coeff, false_exps)
        return [(_one_atom(a.param_info), false_atom)]
    else
        return [(a, _one_atom(a.param_info))]
    end
end

function _split_mixed(s::CSum, group1::BitVector, group2::BitVector, flex::BitVector)
    result = Tuple{CFunction,CFunction}[]
    for term in s.expr
        append!(result, _separate_pairs(term, group1, group2, flex))
    end
    return result
end

function _split_mixed(p::CProd, group1::BitVector, group2::BitVector, flex::BitVector)
    iszero(p.coeff) && return [(zero_catom(p.param_info), _one_atom(p.param_info))]
    blocks = Vector{Vector{Tuple{CFunction,CFunction}}}()
    if !isone(p.coeff)
        push!(blocks, [(_const_atom(p.param_info, p.coeff), _one_atom(p.param_info))])
    end
    for factor in p.expr
        push!(blocks, _separate_pairs(factor, group1, group2, flex))
    end
    return _multiply_pair_blocks(blocks, p.param_info)
end

function _split_mixed(e::CExp, group1::BitVector, group2::BitVector, flex::BitVector)
    base_pairs = _separate_pairs(e.expr, group1, group2, flex)
    if any(! (isone(tp) || isone(fp)) for (tp, fp) in base_pairs)
        throw(ArgumentError("cannot separate exp of an expression mixing parameter groups"))
    end
    true_terms = CFunction[]
    false_terms = CFunction[]
    for (tp, fp) in base_pairs
        if isone(fp)
            push!(true_terms, tp)
        elseif isone(tp)
            push!(false_terms, fp)
        end
    end
    true_inner = _combine_sum(e.param_info, true_terms)
    false_inner = _combine_sum(e.param_info, false_terms)

    true_factor = iszero(true_inner) ? _const_atom(e.param_info, e.coeff) : CExp(e.param_info, e.coeff, true_inner, Val(:nosimp))
    false_factor = iszero(false_inner) ? _one_atom(e.param_info) : CExp(e.param_info, CR_ONE, false_inner, Val(:nosimp))
    return [(true_factor, false_factor)]
end

function _split_mixed(p::CPower, group1::BitVector, group2::BitVector, flex::BitVector)
    base_pairs = _separate_pairs(p.expr, group1, group2, flex)
    length(base_pairs) == 1 || throw(ArgumentError("cannot separate power whose base expands into multiple mixed terms"))
    bt, bf = base_pairs[1]
    if isone(bf)
        true_factor = CPower(p.param_info, p.coeff, bt, p.exponent, Val(:nosimp))
        return [(true_factor, _one_atom(p.param_info))]
    elseif isone(bt)
        false_factor = CPower(p.param_info, CR_ONE, bf, p.exponent, Val(:nosimp))
        coeff_atom = _const_atom(p.param_info, p.coeff)
        return [(coeff_atom, false_factor)]
    elseif denominator(p.exponent) == 1
        true_factor = CPower(p.param_info, p.coeff, bt, p.exponent, Val(:nosimp))
        false_factor = CPower(p.param_info, CR_ONE, bf, p.exponent, Val(:nosimp))
        return [(true_factor, false_factor)]
    else
        throw(ArgumentError("cannot distribute non-integer power across mixed parameter groups"))
    end
end

function _split_mixed(f::CFunction, ::BitVector, ::BitVector, ::BitVector)
    throw(ArgumentError("separation for expressions of type $(typeof(f)) is not implemented"))
end

function _multiply_pair_blocks(blocks::Vector{Vector{Tuple{CFunction,CFunction}}}, param_info::ParameterInfo)
    isempty(blocks) && return [(_one_atom(param_info), _one_atom(param_info))]
    acc = [(_one_atom(param_info), _one_atom(param_info))]
    for block in blocks
        acc = _multiply_pair_lists(acc, block)
        isempty(acc) && break
    end
    return acc
end

function _multiply_pair_lists(lhs::Vector{Tuple{CFunction,CFunction}}, rhs::Vector{Tuple{CFunction,CFunction}})
    buffer = Tuple{CFunction,CFunction}[]
    for (lt, lf) in lhs
        for (rt, rf) in rhs
            new_t = lt * rt
            new_f = lf * rf
            if iszero(new_t) || iszero(new_f)
                continue
            end
            push!(buffer, (new_t, new_f))
        end
    end
    return buffer
end

function _normalize_pair(left::CFunction, right::CFunction)
    coeffs = coeff(left)
    if length(coeffs) == 1
        c = coeffs[1]
        if c == CR_ONE
            return left, right
        elseif iszero(c)
            return left, right * c
        else
            return modify_coeff(left, CR_ONE), right * c
        end
    end
    return left, right
end
