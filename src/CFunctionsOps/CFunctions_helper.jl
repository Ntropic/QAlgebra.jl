# tree_iter_sub iterates for a CCustomtype exclusively 
_tree_iter_sub(n::CAtom,      ::CCustomType) =  (n,)
function _tree_iter_sub(n::CAbstract,  c::CCustomType)
    argpos = c.ctype_def.index_map[n.index]
    return tree_iter(c.expr[argpos])
end
_tree_iter_sub(n::CComposite, c::CCustomType) = Iterators.flatten(((n,), _tree_iter_sub(n.expr, c)))
_tree_iter_sub(n::CMultiComposite, c::CCustomType) = Iterators.flatten(((n,), (_tree_iter_sub(ch, c) for ch in n.expr)))
_tree_iter_sub(n::CRational,  c::CCustomType) = Iterators.flatten(((n,), _tree_iter_sub(n.numer, c), _tree_iter_sub(n.denom, c)))
_tree_iter_sub(v::CVector,    c::CCustomType) = Iterators.flatten(((v,), (_tree_iter_sub(ch, c) for ch in v.expr)))
_tree_iter_sub(M::CMatrix,    c::CCustomType) = Iterators.flatten(((M,), (_tree_iter_sub(ch, c) for ch in M.entries[:])))


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
function tree_iter(c::CCustomType) 
    if c.ctype_def.has_abstract 
        return Iterators.flatten(((c,), _tree_iter_sub(c.ctype_def.fun, c))) 
    else
        return Iterators.flatten(((c,), (tree_iter(ch) for ch in c.expr)))
    end
end

_leaf_iter_sub(n::CAtom,      ::CCustomType) = (n,)
function _leaf_iter_sub(n::CAbstract,  c::CCustomType)
    argpos = c.ctype_def.index_map[n.index]
    leaf_iter(c.expr[argpos])
end
_leaf_iter_sub(n::CComposite, c::CCustomType) = _leaf_iter_sub(n.expr, c)
_leaf_iter_sub(n::CMultiComposite, c::CCustomType) = Iterators.flatten((_leaf_iter_sub(ch, c) for ch in n.expr))
_leaf_iter_sub(n::CRational,  c::CCustomType) = Iterators.flatten((_leaf_iter_sub(n.numer, c), _leaf_iter_sub(n.denom, c)))
_leaf_iter_sub(v::CVector,    c::CCustomType) = Iterators.flatten((_leaf_iter_sub(ch, c) for ch in v.expr))
_leaf_iter_sub(M::CMatrix,    c::CCustomType) = Iterators.flatten((_leaf_iter_sub(ch, c) for ch in M.entries[:]))

"""
    leaf_iter(f::CFunction)

Iterator over the leaves of an expression tree (depth-first).
A "leaf" is either a `CAtom` or a `CAbstract`.
All composite/container nodes are skipped.
"""
leaf_iter(a::CAtom)     = (a,)       # leaf
leaf_iter(a::CAbstract) = (a,)       # leaf
leaf_iter(f::CComposite) = leaf_iter(f.expr)
leaf_iter(f::CMultiComposite) = Iterators.flatten((leaf_iter(ch) for ch in f.expr))
leaf_iter(r::CRational) = Iterators.flatten((leaf_iter(r.numer), leaf_iter(r.denom)))
leaf_iter(v::CVector) = Iterators.flatten((leaf_iter(ch) for ch in v.expr))
leaf_iter(M::CMatrix) = Iterators.flatten((leaf_iter(ch) for ch in M.entries[:]))
function leaf_iter(c::CCustomType)
    if c.ctype_def.has_abstract 
        return _leaf_iter_sub(c.ctype_def.fun, c) 
    else
        return Iterators.flatten((leaf_iter(ch) for ch in c.expr))
    end
end

import ..QAlgebra: unique_sorted!
"""
    contains_which_abstracts(f::CFunction)::Vector{CAbstract}

Returns a Vector with unique and sorted CAbstracts present in the expression tree. 
"""
function contains_which_abstracts(f::CFunction)::Vector{CAbstractDefinition}
    indexes = contains_which_abstract_indexes(f)
    return [f.param_info.abstract_definitions[i] for i in indexes]
end
function contains_which_abstract_indexes(f::CFunction)::Vector{Int}
    all_indexes::Vector{Int} = []
    for leaf in leaf_iter(f) 
        if isa(leaf, CAbstract)
            push!(all_indexes, leaf.index)
        end
    end
    unique_sorted!(all_indexes)
    return all_indexes
end
function abstract_from_abstractdef(defs::CAbstractDefinition)::CAbstract
    return CAbstract(defs.param_info, ComplexRational(1,0,1), defs.index)
end

""" 
    which_ensemble_acting(f::CFunction)::Vector{BitVector}

Returns a vector of vectors of booleans. Each inner vector specifies which of its subsystem indexes are acted upon by the QObj. 
This includes actions from CFunctions. Th function should only be applied after substituting all QAbstract terms. 
Their present can be checked via `contains_abstract(q)`.
"""
function which_ensemble_acting(f::CFunction)::Vector{BitVector}
    where_non_trivial::Vector{BitVector} = [falses(n) for n in f.param_info.how_many_by_ensemble]
    return which_ensemble_acting(f, where_non_trivial)
end
function which_ensemble_acting(f::CFunction, where_non_trivial::Vector{BitVector})::Vector{BitVector}
    for leaf in leaf_iter(f) 
        which_ensemble_acting_atom!(leaf, where_non_trivial)
    end
    return where_non_trivial
end
which_ensemble_acting_atom!(f::CAbstract, where_non_trivial::Vector{BitVector}) = error("Cannot determine the acting ensembles for a CAbstract. Use abstracts only in CType definitions.") 
function which_ensemble_acting_atom!(f::CAtom, where_non_trivial::Vector{BitVector})::Vector{BitVector}
    for (param_ind, where_acting) in zip(f.param_info.indexed_parameter_indexes, f.param_info.where_acting_by_parameter)
        if any(!=(0), f.var_exponents[param_ind])
            vecvec_or!(where_non_trivial, where_acting)
        end
    end
    return where_non_trivial
end

""" 
    where_acting(f::CFunction)::BitVector

Returns a vector of booleans specifying if the associated parameter is present in the expression.
"""
function where_acting(f::CFunction)::BitVector
    acting::BitVector = falses(dims(f))
    where_acting!(f, acting)
    return acting
end
function where_acting!(f::CFunction, acting::BitVector )::BitVector
    for leaf in leaf_iter(f) 
        where_acting_atom!(leaf, acting)
    end
    return acting
end
where_acting_atom!(f::CAbstract, acting::BitVector =[]) = error("Cannot determine where acting for a CAbstract. Use abstracts only in CType definitions. ")
function where_acting_atom!(f::CAtom, acting::BitVector)::BitVector
    # or operation between acting and f.var_exponents being overwritten on acting 
    acting .|= (f.var_exponents .!= 0)
    return acting
end

# Returns index strings, and time strings
function where_acting_to_index_strings(param_indexes::ParameterIndexes, acting::BitVector; do_latex::Bool=false)::Tuple{Vector{String}, Vector{String}}
    current_indexes::Vector{String} = []
    for (str, inds) in zip(param_indexes.labels, param_indexes.label_parameter_indexes)
        if any(acting[inds])
            push!(current_indexes, str)
        end
    end
    current_t_indexes::Vector{String} = []
    t_str_elements = do_latex ? param_indexes.t_labels_latex : param_indexes.t_labels 
    for (str, inds) in zip(t_str_elements, param_indexes.label_parameter_t_indexes)
        if any(acting[inds])
            push!(current_t_indexes, str)
        end
    end
    return current_indexes , current_t_indexes
end
where_acting_to_index_strings(f::CFunction; do_latex::Bool=false)::Tuple{Vector{String}, Vector{String}} = where_acting_to_index_strings(f.param_info.param_indexes, where_acting(f), do_latex=do_latex)

has_indexes(param_indexes::ParameterIndexes, acting::BitVector)::Bool = any(acting[param_indexes.all_indexes])
has_indexes(f::CFunction) = contains_c_indexes(f, f.param_info.param_indexes.all_indexes)


"""
    var_exponents_iter(f::CFunction)

Iterates over the `var_exponents` vectors of all leaves in the expression tree.
- For `CAtom` leaves: yields the actual `var_exponents`.
- For `CAbstract` leaves: yields a zero vector (same length as number of dims).
"""
function var_exponents_iter(f::CFunction)
    Iterators.map(_leaf2exps, leaf_iter(f))
end
_leaf2exps(a::CAtom)     = a.var_exponents
_leaf2exps(a::CAbstract) = zeros(Int, dims(a))

"""
    var_exponents_iter_simple(f::CFunction)

Like `var_exponents_iter`, but only includes exponents from
`CAtom`, `CSum`, and `CProd`.
Other function types contribute trivial exponents.
"""
var_exponents_iter_simple(a::CAtom) = (a.var_exponents,)
var_exponents_iter_simple(s::CSum)  = Iterators.flatten(var_exponents_iter_simple.(s.expr))
var_exponents_iter_simple(p::CProd) = Iterators.flatten(var_exponents_iter_simple.(p.expr))
var_exponents_iter_simple(f::CFunction) = (zeros(Int, dims(f)),)

import Base: isnumeric
"""
    isnumeric(f::CFunction)

True if all exponents in all terms are zero.
"""
isnumeric(f::CAtom) = iszero(f.coeff) || all(==(0), f.var_exponents)
isnumeric(f::CAbstract) = false
isnumeric(f::CFunction) = all(leaf -> isnumeric(leaf) , leaf_iter(f))

import Base: iszero, isempty, isone

"""
    iszero(a::CAtom)     -> Bool
    iszero(s::CSum)      -> Bool
    iszero(r::CRational) -> Bool

Returns `true` if the expression is identically zero:
"""
iszero(a::CAtomic)        = iszero(a.coeff)
iszero(s::CSum)         = all(iszero, s.expr)
iszero(p::CComposite)        = iszero(p.coeff) || iszero(p.expr)
iszero(p::CMultiComposite)        = iszero(p.coeff) || any(iszero, p.expr)
iszero(r::CRational)    = iszero(r.numer)
function iszero(c::CCustomType)
    # 1. First check coeff
    iszero(c.coeff) && return true
    # 2. Substitute arguments into the base definition
    substituted = simplify(substitute(c.ctype_def.fun, c.ctype_def.abstract_parameters, c.expr))
    # 3. Check if the expanded form is zero
    return iszero(substituted)
end

isempty(s::CMultiComposite)        = isempty(s.expr)

isone(c::CFunction) = false 
isone(a::CAtom)     = isnumeric(a) && isone(a.coeff)

# only used for printing the signs! doesn't mean evaluated function is negative! 
allnegative(a::CAtom) = is_negative(a.coeff)
allnegative(s::CSum)  = !isempty(s.expr) && all(allnegative, s.expr)
allnegative(p::CProd) = is_negative(p.coeff)
allnegative(r::CRational) = allnegative(r.numer)
allnegative(x::CFunction)  = is_negative(x.coeff)

"""
    min_exponents(f::CFunction) -> Vector{Int}

Component-wise minimum of all monomial exponent vectors appearing in `f`.
If `f` contains no atoms (e.g. empty containers), returns `Int[]`.

Relies on `var_exponents_iter(::CFunction)`.
"""
function min_exponents(f::CFunction)::Vector{Int}
    mins = zeros(Int, dims(f))           # start with all zeros
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


contains_vec_or_mat(c::CFunction) = any(element -> isa(element, Union{CMatrix, CVector}), tree_iter(c))

# --- numeric printing gate (specify if the coefficient is needed) -----------------------------------------------------
"""
    printnumeric(f::CFunction) -> Bool

Decides if a numeric coefficient should be printed for `f`.
"""
function printnumeric(f::CFunction)::Bool
    if f isa CRational
        return true
    elseif isnumeric(f) && isonelike(coeff(f)[1])
        return false
    elseif f isa CMultiComposite && length(f.expr) == 1
        return printnumeric(f.expr[1])
    else
        return true
    end
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
        if !is_axis_multiple(v[1], v[i])[1]
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
    return !any(term -> typeof(term)==CAtom, f.expr)
end
function simple_CSum(f::CAtom)::Bool
    return true 
end
function simple_CSum(f::CFunction)::Bool
    return false 
end

# factor a simple sum into:  (pre_F) * (new_f)
# Only runs when the sum is of "monomial-like" terms (e.g., atoms) so we can safely
# factor a common scalar and a common exponent offset.
function separate_CSum(f::CSum)::Tuple{Bool, Union{CAtom,Nothing}, CSum}
    # treat as "simple" if every summand is a monomial-like term (e.g. CAtom)
    # (adjust this predicate if you later support more term kinds)
    if all(t -> t isa CAtom, f.expr)
        # scalar coeffs of each term
        coeffs = ComplexRational[first(coeff(t)) for t in f.expr]

        if has_common_multiple(coeffs)
            base, multiples = common_denominator_form(coeffs)

            # exponent vectors of each term
            vs = var_exponents.(f.expr)
            offset = common_exponent_offset(vs)

            # rebuild normalized terms with adjusted coeffs/exponents
            new_terms = [CAtom(f.param_info, m, v .- offset) for (m, v) in zip(multiples, vs)]

            # shared factor in front
            pre_f = CAtom(f.param_info, base, offset)

            # always return a CSum (avoid _CSum so we don't collapse to a single term)
            return true, pre_f, CSum(f.param_info, new_terms)
        end
    end
    return false, nothing, f
end


function simple_combinable_F(t1::T1, t2::T2)::Tuple{Bool, ComplexRational} where {T1<:CFunction, T2<:CFunction}
    v1 = t1 isa CSum ? t1.expr : [t1]
    v2 = t2 isa CSum ? t2.expr : [t2]

    if length(v1) != length(v2)
        return false, ComplexRational(0,0,1)
    end

    # if all terms in v1 are constant in the polynomial variables, skip
    if all(x -> all(e == 0 for e in var_exponents(x)), v1)
        return false, ComplexRational(0,0,1)
    end

    # exponent pattern must match up to a constant offset across corresponding terms
    exponents_diff = var_exponents(v1[1]) .- var_exponents(v2[1])
    if any(!=(0), exponents_diff)
        return false, ComplexRational(0,0,1)
    end
    for (a1, a2) in zip(v1[2:end], v2[2:end])
        if var_exponents(a1) .- var_exponents(a2) != exponents_diff
            return false, ComplexRational(0,0,1)
        end
    end

    # constant scalar ratio of the coefficients across corresponding terms
    ratio = first(coeff(v1[1])) / first(coeff(v2[1]))
    # keep your original requirement that the ratio is purely real or purely imaginary
    if !(ratio.b == 0 || ratio.a == 0)
        return false, ratio
    end
    for (a1, a2) in zip(v1[2:end], v2[2:end])
        if first(coeff(a1)) / first(coeff(a2)) != ratio
            return false, ratio
        end
    end
    return true, ratio
end


# group inputs into buckets of pairwise-combinable terms
function simple_combinable_Fs(ts::AbstractVector{<:CFunction})
    groups  = Vector{Vector{CFunction}}()
    indexes = Vector{Vector{Int}}()
    for (i, t) in enumerate(ts)
        placed = false
        for (inds, grp) in zip(indexes, groups)
            ok, _ = simple_combinable_F(grp[1], t)
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


# ratios between the leading coefficients of each element in a group
function ratios_Fs(ts::AbstractVector{<:CFunction})::Vector{ComplexRational}
    # take the first summand if it's a CSum; otherwise the term itself
    first_term = ts[1] isa CSum ? ts[1].expr[1] : ts[1]
    c1 = first(coeff(first_term))

    ratios = ComplexRational[ComplexRational(1,0,1)]
    for t in ts[2:end]
        inner = t isa CSum ? t.expr[1] : t
        c2 = first(coeff(inner))
        push!(ratios, c1 / c2)
    end
    return ratios
end


"""
    group_Fs(ts)

Given a combinable group `ts`, decide how to express it as either:
  - a single `CFunction`, or
  - a tuple `(F1, Fs)` meaning `F1 * (Fs[1] + ... + Fs[n])`.

`Fs[i]` are returned as scalar monomials (atoms with zero exponents) so that
they act purely as numeric multipliers in the grouped sum.
"""
function group_Fs(ts::AbstractVector{<:CFunction})::Union{CFunction, Tuple{CFunction, Vector{CFunction}}}
    if length(ts) == 1
        return ts[1]
    else
        ratios = ratios_Fs(ts)
        base, multiples = common_denominator_form(ratios)

        # pre-factor: first element scaled by the 'base' ratio
        pre_F = ts[1] * base

        # zero-exponent vector for creating pure scalar atoms
        # use the first representative's exponent length as dimension
        rep = ts[1] isa CSum ? ts[1].expr[1] : ts[1]
        dim = length(var_exponents(rep))
        zexp = zeros(Int, dim)

        post_Fs = CFunction[ CAtom(rep.param_info, m, zexp) for m in multiples ]
        return (pre_F, post_Fs)
    end
end


"""
    how_to_combine_Fs(ts::Vector{CFunction})

Returns:
  - `groups_as_functions::Vector{Union{CFunction, Tuple{CFunction, Vector{CFunction}}}}`
  - `indexes::Vector{Vector{Int}}` (original positions per group)
"""
function how_to_combine_Fs(ts::Vector{CFunction})
    if isempty(ts)
        return CFunction[], Vector{Vector{Int}}()
    elseif length(ts) == 1
        return [ts[1]], [[1]]
    else
        groups, indexes = simple_combinable_Fs(ts)
        return [group_Fs(grp) for grp in groups], indexes
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
    return changed_any, CAtom(atom.param_info, atom.coeff, new_exponents)
end

function term_equal_indexes(fsum::CSum, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CSum}
    changed_any = false
    new_terms = Vector{CFunction}(undef, length(fsum.expr))
    @inbounds for k in eachindex(fsum.expr)
        changed, new_term = term_equal_indexes(fsum.expr[k], coeff_ind_order)
        changed_any |= changed
        new_terms[k] = new_term
    end
    return changed_any, _CSum(new_terms)
end

function term_equal_indexes(fractional::CRational, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CRational}
    changed_num, new_num = term_equal_indexes(fractional.numer, coeff_ind_order)
    changed_den, new_den = term_equal_indexes(fractional.denom, coeff_ind_order)
    return (changed_num || changed_den), CRational(fractional.param_info, new_num, new_den, Val(:nosimp))
end

# NEW
function term_equal_indexes(fprod::CProd, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CFunction}
    changed_any = false
    new_terms = Vector{CFunction}(undef, length(fprod.expr))
    @inbounds for k in eachindex(fprod.expr)
        changed, new_term = term_equal_indexes(fprod.expr[k], coeff_ind_order)
        changed_any |= changed
        new_terms[k] = new_term
    end
    return changed_any, CProd(fprod.param_info, fprod.coeff, new_terms, Val(:nosimp))
end

function term_equal_indexes(fexp::CExp, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CFunction}
    changed, nx = term_equal_indexes(fexp.expr, coeff_ind_order)
    return changed, CExp(fexp.param_info, fexp.coeff, nx, Val(:nosimp))
end

function term_equal_indexes(flog::CLog, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CFunction}
    changed, nx = term_equal_indexes(flog.expr, coeff_ind_order)
    return changed, CLog(flog.param_info, flog.coeff, nx, Val(:nosimp))
end

function term_equal_indexes(fpwr::CPower, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CFunction}
    changed, nx = term_equal_indexes(fpwr.expr, coeff_ind_order)
    return changed, CPower(fpwr.param_info, fpwr.coeff, nx, fpwr.exponent, Val(:nosimp))
end

function term_equal_indexes(v::CVector, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CVector}
    changed_any = false
    new_entries = Vector{CFunction}(undef, length(v.expr))
    @inbounds for k in eachindex(v.expr)
        changed, e = term_equal_indexes(v.expr[k], coeff_ind_order)
        changed_any |= changed
        new_entries[k] = e
    end
    return changed_any, CVector(v.param_info, v.coeff, new_entries; row=v.row)
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
    return changed_any, CMatrix(M.param_info, M.coeff, new_mat)
end

function term_equal_indexes(A::CAbstract, coeff_ind_order::Vector{Tuple{Int, Int}}) 
    error("Cannot substitute indexes in Abstract expressions. ")
end

function term_equal_indexes(C::CCustomType, coeff_ind_order::Vector{Tuple{Int, Int}})::Tuple{Bool, CCustomType}
    changed_any = false
    new_parameters::Vector{CFunction} = []
    for x in C.expr
        c, new_x = term_equal_indexes(x, var_tuples)
        push!(new_parameters, new_x)
        changed_any |= c
    end
    return changed_any, modify_expr(f, new_parameters)
end