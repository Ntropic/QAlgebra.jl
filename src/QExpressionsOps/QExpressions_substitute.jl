export Substitution, Substitution_t, Substitution_index, -->

import ..CFunctions: substitute
using ..QSpaces: map_by_subspace, map_by_tindex



# Handy alias for the concrete mapping used below
const _IndexToken = Union{SubSpaceIndex, Symbol, AbstractString}
"""
    Substitution_index(from, to)

Subspace-index substitution mapping `from` → `to`.

Tokens can be `SubSpaceIndex`, `Symbol` or `String` corresponding to index labels
in the ambient `QSpace`.
"""
struct Substitution_index
    from::_IndexToken
    to::_IndexToken
end


# --- typed Substitution + ASCII operator ------------------------------------------
"""
    Substitution(from_expr::QExpr, to_expr::QExpr)

Create a substitution struct from `from_expr` → `to_expr`.  
Requires both to be in the same `QSpace`.
"""
struct Substitution
    from::QAbstract
    to::QAtom
    qspace::QSpace
end

"""
    Substitution_t(from, to)

Time-index substitution mapping `from` → `to`.

Accepted tokens: integers (≥ 0) and symbols/strings such as `:t`, `:t1`, `:t_2`.
"""
struct Substitution_t
    from::Int
    to::Int
    function Substitution_t(from::Union{Int,Symbol,AbstractString},
                            to::Union{Int,Symbol,AbstractString})
        from_idx = _parse_time_token(from)
        to_idx   = _parse_time_token(to)
        return new(from_idx, to_idx)
    end
end

# Primary constructor: take two QExprs, extract (QAbstract -> QAtom), check qspace
function Substitution(from_expr::QExpr, to_expr::QExpr)
    if is_t(from_expr) && is_t(to_expr)
        from_idx = _time_index_from_expr(from)
        from_idx === nothing && error("Cannot infer time index from expression $from_expr.")
        to_idx = _time_index_from_expr(to)
        to_idx === nothing && error("Cannot infer time index from expression $to_expr.")
        return Substitution_t(from_idx, to_idx)
    end
    from_abs = extract_qabstract(from_expr)
    to_atom  = extract_qatom(to_expr)
    ss_from  = from_expr.qspace
    ss_to    = to_expr.qspace
    ss_from === ss_to || error("Statespace mismatch between 'from' and 'to'.")
    return Substitution(from_abs, to_atom, ss_from)
end
function -->(from::Union{QExpr,QAtomProduct}, to::Union{QExpr,QAtomProduct})
    from_idx = _time_index_from_expr(from)
    to_idx = _time_index_from_expr(to)
    if !(from_idx === nothing || to_idx === nothing)
        return Substitution_t(from_idx, to_idx)
    end
    return Substitution(from, to)
end


struct _ResolvedIndexSubstitution
    from::SubSpaceIndex
    to::SubSpaceIndex
end


function -->(from::Union{Symbol,String,SubSpaceIndex}, to::Union{Symbol,String,SubSpaceIndex})
    if _is_time_token(from) && _is_time_token(to)
        return Substitution_t(from, to)
    end
    return Substitution_index(from, to)   # placeholder, gets resolved when we have param info from a quantum object 
end


# decide which type it is.
@inline function _is_time_token(token::AbstractString)::Bool
    str = String(token)
    return str == "t" || occursin(r"^t_?\d+$", str)  # accepts t, t0, t1, t_0, t_1 
end
@inline _is_time_token(token::Symbol)::Bool = _is_time_token(String(token))
@inline _is_time_token(::SubSpaceIndex)::Bool = false

@inline _parse_time_token(token::Symbol) = _parse_time_token(String(token))
function _parse_time_token(token::AbstractString)::Int
    str = String(token)
    str == "t" && return 0
    m = match(r"^t_?(\d+)$", str)
    m === nothing && error("Invalid time index token '$token'. Expected forms like :t, :t1, :t_2.")
    return parse(Int, m.captures[1])
end
@inline function _parse_time_token(token::Integer)::Int
    token < 0 && error("Time indexes must be non-negative, got $token.")
    return Int(token)
    return _ResolvedIndexSubstitution(from, to)
end
_resolve_index_token(token::SubSpaceIndex, subspace_info::SubSpaceInfo) = token
_resolve_index_token(token::Symbol, subspace_info::SubSpaceInfo) = SubSpaceIndex(token, subspace_info)

@inline function _resolve_index_token(token::AbstractString, subspace_info::SubSpaceInfo)
    return SubSpaceIndex(token, subspace_info)
end

function _resolve_index_substitution(sp::Substitution_index, subspace_info::SubSpaceInfo)
    from = _resolve_index_token(sp.from, subspace_info)
    to   = _resolve_index_token(sp.to, subspace_info)
    from.outer == to.outer || error("Index substitution must stay within the same subspace (got $from → $to).")
    return _ResolvedIndexSubstitution(from, to)
end

@inline function _changed_pairs(mapvec::Vector{Int})::Vector{Tuple{Int,Int}}
    pairs = Tuple{Int,Int}[]
    for i in eachindex(mapvec)
        val = mapvec[i]
        val == i && continue
        push!(pairs, (i, val))
    end
    return pairs
end

@inline _coeff_pairs(sp::Substitution_t, qspace::QSpace) = _changed_pairs(map_by_tindex(sp.from, sp.to, qspace.param_info))
@inline _coeff_pairs(sp::_ResolvedIndexSubstitution, qspace::QSpace) = _changed_pairs(map_by_subspace(sp.from, sp.to, qspace.param_info))

@inline function _reindex_coeff(coeff::CFunctions.CFunction, pairs::Vector{Tuple{Int,Int}})
    isempty(pairs) && return (false, coeff)
    return CFunctions.term_equal_indexes(coeff, pairs)
end

@inline function _time_index_from_coeff(coeff::CFunctions.CFunction, ::QSpace)
    return nothing
end

@inline function _time_index_from_coeff(coeff::CFunctions.CAtom, qspace::QSpace)
    isone(coeff.coeff) || return nothing
    exps = coeff.var_exponents
    idxs = findall(!iszero, exps)
    length(idxs) == 1 || return nothing
    idx = idxs[1]
    exps[idx] == 1 || return nothing
    param = qspace.params[idx]
    param.is_t || return nothing
    return param.t_index
end

@inline _time_index_from_atomproduct(prod::QAtomProduct) =
    (is_t(prod) ? _time_index_from_coeff(prod.coeff_fun, prod.qspace) : nothing)

@inline _time_index_from_expr(prod::QAtomProduct) = _time_index_from_atomproduct(prod)

function _time_index_from_expr(expr::QExpr)
    length(expr) == 1 || return nothing
    term = expr.terms[1]
    term isa QAtomProduct || return nothing
    return _time_index_from_atomproduct(term)
end
@inline _identity_op(qspace::QSpace, expanded::Int) = qspace.I_op[expanded]

function extract_qabstract(q::QExpr)::QAbstract
    if length(q) > 1 
        error("Cannot substitute composites. abstract_op must contain only an abstract operator.")
    end
    term = q.terms[1]
    if !isa(term, QAtomProduct)
        error("QExpr must contain only a QAtomProduct")
    end
    return extract_qabstract(term)
end
function extract_qabstract(term::QAtomProduct)::QAbstract
    if length(term.expr) != 1 || !isa(term.expr[1], QAbstract)
        error("abstract_op must contain exactly one QAbstract")
    end
    return term.expr[1]
end
#now same for QAtom 
function extract_qatom(q::QExpr)::QAtom
    if length(q) > 1 
        error("Cannot substitute composites. abstract_op must contain only an abstract operator.")
    end
    term = q.terms[1]
    if !isa(term, QAtomProduct)
        error("QExpr must contain only a QAtomProduct")
    end
    return extract_qatom(term)
end
function extract_qatom(term::QAtomProduct)::QAtom
    if length(term.expr) != 1 || !isa(term.expr[1], QAtom)
        error("abstract_op must contain exactly one QAbstract")
    end
    return term.expr[1]
end

# --- helpers ------------------------------------------------------------------
# Deals with index_map mapping one index to another within QAbstract
function qAtom_index_flip(q::QAtom, index_map::Vector{Tuple{SubSpaceIndex,SubSpaceIndex}}, qspace::QSpace)::Vector{QAtomProduct}
    qs::Vector{QAtom} = [q]
    cs::Vector{ComplexRational} = [ComplexRational(1,0,1)]
    for (index1, index2) in index_map
        new_qs::Vector{QAtom} = []    
        new_cs::Vector{ComplexRational} = []
        for (qi, ci) in zip(qs, cs)
            _, new_terms, new_coeffs = term_equal_indexes( qi, index1.expanded, index2.expanded, qspace.subspaces[index1.inner])
            append!(new_qs, new_terms)
            append!(new_cs, new_coeffs*ci)
        end
        qs = new_qs
        cs = new_cs
    end
    return [QAtomProduct(qspace, qspace.c_one*c, [q]) for (q, c) in zip(qs, cs)]
end

# --- substitution on single terms (new order: target first, then Substitution) -----
# QTerm that is NOT the abstract operator: no change
function substitute_qAtom(target::QTerm, sp::Substitution)::Vector{QAtomProduct}
    ss = sp.qspace
    return [QAtomProduct(ss, ss.c_one, [target])]
end

# QAbstract: do the replacement when it matches the 'from'
function substitute_qAtom(target::QAbstract, sp::Substitution)::Vector{QAtomProduct}
    a = sp.from
    r = sp.to
    ss = sp.qspace
    if target.key_index == a.key_index && target.sub_index == a.sub_index
        qs = QExpr(ss, qAtom_index_flip(r, target.index_map, ss))
        if target.exponent != 1
            qs = qs^target.exponent
        end
        if target.dag
            qs = qs'
        end
        return qs.terms
    else
        return [QAtomProduct(ss, ss.c_one, [target])]
    end
end

# General unresolved substitution type
function substitute(target::T, sp::Substitution_index)::T where T <: Union{QExpr, diffQEq}
    sp_neq = _resolve_index_substitution(sp, target.qspace.subspace_info) # figure out subspace_info along the way 
    return substitute(target, sp_neq) 
end
# -------------------------  Substitution of QAbstract with QAtom -----------------------------------------------

substitution_properties_fulfilled(sub::Substitution)::Bool = substitution_properties_fulfilled(sub.from , QExpr(sub.qspace, sub.to))

"""
    substitute(T, sp::Substitution; checks=false) where T <: QObj

Apply a substitution `sp::Substitution` (alias `Substitution`) to different target types:

Optionally, property checks can be enabled via `checks=true`.
"""
function substitute(target::T, sp::Substitution)::T where T <: Union{QExpr, diffQEq}
    substitution_properties_fulfilled(sp)
    target.qspace === sp.qspace || error("Statespace mismatch between target and Substitution.")
    return _substitute(target, sp) 
end
function _substitute(target::QExpr, sp::Substitution)
    new_terms = QComposite[]
    for term in target.terms
        append!(new_terms, _substitute(term, sp))
    end
    return QExpr(target.qspace, new_terms)
end

function _substitute(target::QAtomProduct, sp::Substitution)::Vector{QComposite}
    ss = target.qspace
    expr = target.expr
    coeff_fun = target.coeff_fun

    new_expr::QExpr = QExpr(ss, substitute_qAtom(expr[1], sp))
    for t in expr[2:end]
        new_terms = substitute_qAtom(t, sp)
        new_new_expr = new_expr * new_terms[1]
        for tt in new_terms[2:end]
            new_new_expr += new_expr + tt
        end
        new_expr = new_new_expr
    end
    terms = new_expr.terms
    return [QAtomProduct(ss, coeff_fun*t.coeff_fun, t.expr) for t in terms]
end


# Any single composite holding one expr
function _substitute(targ::T, sp::Substitution) where {T<:QComposite}
    # Note: QMultiComposite is <: QComposite; a more specific method follows below.
    return modify_expr(targ, _substitute(targ.expr, sp))
end

function _substitute(target::QSum, sp::Substitution)
    new_expr = _substitute(target.expr, sp)
    new_expr === target.expr && return QComposite[target]
    return modify_expr(target, new_expr, Val(:nodecollision))
end

# Any multi-composite holding many sub-expressions
function _substitute(targ::T, sp::Substitution) where {T<:QMultiComposite}
    substituted = map(x -> _substitute(x, sp), targ.expr)
    return modify_expr(targ, map(only, substituted))
end

# diffQEq → diffQEq
function _substitute(target::diffQEq, sp::Substitution)::diffQEq
    lhs = _substitute(target.left_hand_side, sp) 
    if length(lhs) != 1
        error("Substitution of $(sp.from) with $(sp.to) in $target did not result in a single term.")
    end
    rhs = _substitute(target.expr, sp)
    return diffQEq(target.qspace, lhs[1], rhs)
end

# --- substitution contexts --------------------------------------------------------

abstract type AbstractSubContext end

# Cache for one substitution over time indices; reused by all nested calls.
struct TimeSubContext <: AbstractSubContext
    from::Int
    to::Int
    coeff_pairs::Vector{Tuple{Int,Int}}
end

# Cache for one substitution over ensemble indexes, including expanded slots.
struct IndexSubContext <: AbstractSubContext
    qspace::QSpace
    from::SubSpaceIndex
    to::SubSpaceIndex
    from_exp::Int
    to_exp::Int
    coeff_pairs::Vector{Tuple{Int,Int}}
end

@inline function _check_time_bounds(sp::Substitution_t, qspace::QSpace)
    max_t = qspace.max_t_ind
    sp.from <= max_t || error("Time index $(sp.from) out of bounds for QSpace with max_t_ind=$(max_t).")
    sp.to   <= max_t || error("Time index $(sp.to) out of bounds for QSpace with max_t_ind=$(max_t).")
end

function _time_context(qspace::QSpace, sp::Substitution_t, checks::Bool)::TimeSubContext
    checks && _check_time_bounds(sp, qspace)
    return TimeSubContext( sp.from, sp.to, _coeff_pairs(sp, qspace))
end

function _index_context(qspace::QSpace, sp::_ResolvedIndexSubstitution)::IndexSubContext
    return IndexSubContext(qspace, sp.from, sp.to, sp.from.expanded, sp.to.expanded, _coeff_pairs(sp, qspace))
end

function _index_context(qspace::QSpace, sp::Substitution_index)::IndexSubContext
    resolved = _resolve_index_substitution(sp, qspace.subspace_info)
    return _index_context(qspace, resolved)
end

# --- time index substitution -------------------------------------------------------

function substitute(target::T, sp::Substitution_t; checks::Bool=false)::T where T <: Union{QExpr, diffQEq}
    ctx = _time_context(target.qspace, sp, checks)
    return _substitute(target, ctx)
end

function _substitute(target::QExpr, ctx::TimeSubContext, )::QExpr
    terms = target.terms
    changed = false
    new_terms = Vector{QComposite}(undef, length(terms))
    for (i, term) in enumerate(terms)
        new_term = _substitute(term, ctx)
        changed |= (new_term !== term)
        new_terms[i] = new_term
    end
    return changed ? QExpr(target.qspace, new_terms) : target
end

function _substitute(target::diffQEq, ctx::TimeSubContext)::diffQEq
    new_lhs = _substitute( target.left_hand_side, ctx)
    new_rhs = _substitute( target.expr, ctx)
    return diffQEq(target.qspace, new_lhs, new_rhs)
end

function _substitute(term::QTerm, ctx::TimeSubContext)::QTerm
    term.time_index == ctx.from || return term
    ctx.checks && term.time_index == -1 && error("Cannot change time index of a time-independent QTerm.")
    return modify_time_index(term, ctx.to)
end

function _substitute(op::QAbstract, ctx::TimeSubContext)::QAbstract
    op.time_index == ctx.from || return op
    ctx.checks && !of_time(op) && error("Cannot change time index of time-independent operator $(op.operator_type.name).")
    return modify_time_index(op, ctx.to)
end

function _substitute(target::QAtomProduct, ctx::TimeSubContext)::QAtomProduct
    expr = target.expr
    new_expr = Vector{QAtom}(undef, length(expr))
    for (i, atom) in enumerate(expr)
        new_atom = _substitute(atom, ctx)
        new_expr[i] = new_atom
    end
    _, new_coeff = _reindex_coeff(target.coeff_fun, ctx.coeff_pairs)
    return modify_coeff_expr(target, new_coeff, new_expr)
end

function _substitute(target::QSum, ctx::TimeSubContext)::QSum 
    new_expr = _substitute( target.expr, ctx)
    return only(modify_expr(target, new_expr))
end

function _substitute(target::T, ctx::TimeSubContext)::T where T <: QComposite
    new_expr = _substitute( target.expr, ctx)
    _, new_coeff = _reindex_coeff(target.coeff_fun, ctx.coeff_pairs)
    return modify_coeff_expr(target, new_coeff, new_expr)
end
#function build_QCompositeProduct(terms::Vector{QComposite}, qspace::QSpace)::QCompositeProduct # continue here
#function _substitute(target::QCompositeProduct, ctx::TimeSubContext)::QCompositeProduct
#    new_expr = _substitute( target.expr, ctx)
#    _, new_coeff = _reindex_coeff(target.coeff_fun, ctx.coeff_pairs)
#    @assert length(new_expr) > 1 "QCompositeProducts with length 1 shouldn't exist."
#end

function _substitute(target::T, ctx::TimeSubContext)::T where T <: QMultiComposite
    new_expr = Vector{QComposite}(undef, length(target.expr))
    for (i, comp) in enumerate(target.expr)
        new_comp = _substitute( comp, ctx)
        new_expr[i] = new_comp
    end
    _, new_coeff = _reindex_coeff(target.coeff_fun, ctx.coeff_pairs)
    return modify_coeff_expr(target, new_coeff,  new_expr)
end


# --- index substitution -----------------------------------------------------------

function substitute(target::T, sp::_ResolvedIndexSubstitution)::T where T <: Union{QExpr, diffQEq}
    ctx = _index_context(target.qspace, sp)
    return _substitute(target, ctx)
end

function _substitute(target::QExpr, ctx::IndexSubContext, )::QExpr
    terms = target.terms
    new_terms = Vector{QComposite}(undef, length(terms))
    for (i, term) in enumerate(terms)
        new_terms[i] = _substitute(term, ctx)
    end
    return QExpr(target.qspace, new_terms)
end
function _substitute(target::diffQEq, ctx::IndexSubContext)::diffQEq
    new_lhs = _substitute(target.left_hand_side, ctx)
    new_rhs = _substitute(target.expr, ctx)
    return diffQEq(target.qspace, new_lhs, new_rhs)
end

function _substitute(term::QTerm, ctx::IndexSubContext)::QTerm
    from_exp = ctx.from_exp
    to_exp = ctx.to_exp
    from_exp == to_exp && return term
    op_indices = term.op_indices
    op_from = op_indices[from_exp]
    id_from = _identity_op(ctx.qspace, from_exp)
    op_from == id_from && return term
    op_to = op_indices[to_exp]
    id_to = _identity_op(ctx.qspace, to_exp)
    if op_to != id_to && op_to != op_from
        error("Target index $(ctx.to) already carries a different operator; cannot substitute $(ctx.from) → $(ctx.to).")
    end
    new_indices = copy(op_indices)
    new_indices[to_exp] = op_from
    new_indices[from_exp] = id_from
    return QTerm(new_indices, term.time_index)
end

function _substitute(op::QAbstract, ctx::IndexSubContext)::QAbstract
    ctx.from == ctx.to && return op
    pair = (ctx.from, ctx.to)
    if pair in op.index_map
        return op
    end
    return add_to_index_map(op, pair)
end

function _substitute(target::QAtomProduct, ctx::IndexSubContext)::QAtomProduct
    expr = target.expr
    changed = false
    new_expr = Vector{QAtom}(undef, length(expr))
    for (i, atom) in enumerate(expr)
        new_atom = _substitute(atom, ctx)
        changed |= (new_atom !== atom)
        new_expr[i] = new_atom
    end
    coeff_changed, new_coeff = _reindex_coeff(target.coeff_fun, ctx.coeff_pairs)
    changed || coeff_changed || return target
    expr_out = changed ? new_expr : expr
    coeff_out = coeff_changed ? new_coeff : target.coeff_fun
    return QAtomProduct(target.qspace, coeff_out, expr_out, target.separate_expectation_values)
end


function _substitute(target::QSum, ctx::IndexSubContext)::QSum
    new_expr = _substitute(target.expr, ctx)
    qspace = target.qspace
    info = qspace.subspace_info
    ensemble = info.ensemble_index_by_outer_index[ctx.from.outer]
    # ensemble != 0 || error("Index $(Index2String(ctx.from, info)) does not belong to an ensemble subspace.")
    target_block = target.blocks[ensemble]
    positions = findall(idx -> idx.expanded == ctx.from_exp, target_block.indexes)
    isempty(positions) && return QSum(qspace, new_expr, target.blocks)

    info.ensemble_index_by_outer_index[ctx.to.outer] == ensemble || error("Cannot substitute index $(Index2String(ctx.from, info)) with $(Index2String(ctx.to, info)): different ensemble blocks.")

    blocks = clone_blocks(target.blocks)
    block = blocks[ensemble]

    old_inner = ctx.from.inner
    new_inner = ctx.to.inner

    if old_inner != new_inner
        _swap_constraint_columns!(block.constraints, old_inner, new_inner)
    end

    for pos in positions
        block.indexes[pos] = ctx.to
        row = block.constraints[pos]
        row[new_inner] = true
    end

    perm = sortperm(block.indexes; by=expanded)
    if !isempty(perm)
        sorted_indexes = block.indexes[perm]
        sorted_constraints = block.constraints[perm]
        block.indexes[:] = sorted_indexes
        block.constraints[:] = sorted_constraints
    end

    return QSum(qspace, new_expr, blocks)
end

function _substitute(target::T, ctx::IndexSubContext)::T where T <: QComposite 
    new_expr = _substitute(target.expr, ctx)
    _, new_coeff = _reindex_coeff(target.coeff_fun, ctx.coeff_pairs)
    return modify_coeff_expr(target, new_coeff, new_expr) 
end


function _substitute(target::T, ctx::IndexSubContext)::T where T <: QMultiComposite
    new_expr = _substitute.(target.expr, Ref(ctx))
    _, new_coeff = _reindex_coeff(target.coeff_fun, ctx.coeff_pairs)
    return QCumulant(target.qspace, coeff_out, new_atom, new_expr, target.order, target.where_acting)
end
