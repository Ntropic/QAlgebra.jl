export Dag, Commutator


const _CHK   = Val(true)
const _NOCHK = Val(false)
# Default-to-checked helper entry points
@inline _add(a, b) = _add(a, b, _CHK)
@inline _sub(a, b) = _sub(a, b, _CHK)
@inline _mul(a, b) = _mul(a, b, _CHK)
@inline _comm(a, b) = _comm(a, b, _CHK)

### Basic Operations: 

####### Unary Minus #############################################################
function -(t::QAtomProduct)::QAtomProduct
    return QAtomProduct(t.statespace, -t.coeff_fun, t.expr) # removed copy
end
function -(t::QExpr)::QExpr
    return QExpr(t.statespace, .-t.terms)
end
function -(t::T)::T where T<:QComposite
    return modify_expr(t, -t.expr)
end

#### Binary + ####################################################################
# Public APIs (checked by default)
+(Q1::QExpr, Q2::QExpr)::QExpr = _add(Q1, Q2)
(+(Q1::QExpr, Q2::T)::QExpr) where {T<:QComposite} = _add(Q1, Q2)
+(Q1::QExpr, N::Number)::QExpr = _add(Q1, N)
+(N::Number, Q1::QExpr)::QExpr = _add(N, Q1)

# Helpers (Val-driven)
@inline function _add(Q1::QExpr, Q2::QExpr, ::Val{C})::QExpr where {C}
    statespace_check_if(Val(C), Q1, Q2)
    new_terms = simplify_QExpr(vcat(Q1.terms, Q2.terms))
    return QExpr(Q1.statespace, new_terms)
end
@inline function _add(Q1::QExpr, Q2::T, ::Val{C})::QExpr where {T<:QComposite,C}
    statespace_check_if(Val(C), Q1, Q2)
    return QExpr(Q1.statespace, vcat(Q1.terms, Q2))
end
# Number interactions: no check
@inline function _add(Q1::QExpr, N::Number, ::Val{C})::QExpr where {C}
    new_terms = vcat(Q1.terms, QAtomProduct(Q1.statespace, Q1.statespace.c_one * N, QTerm(Q1.statespace.I_op)))
    return QExpr(Q1.statespace, new_terms)
end
@inline _add(N::Number, Q1::QExpr, ::Val{C}) where {C} = _add(Q1, N, _NOCHK)


#### Binary - ####################################################################
-(Q1::QExpr, Q2::QExpr)::QExpr = _sub(Q1, Q2)
(-(Q1::QExpr, Q2::T)::QExpr) where {T<:QComposite} = _sub(Q1, Q2)
(-(Q2::T, Q1::QExpr)::QExpr) where {T<:QComposite} = -Q1 + Q2
-(Q1::QExpr, N::Number)::QExpr = _sub(Q1, N)
-(N::Number, Q1::QExpr)::QExpr = _sub(N, Q1)

@inline function _sub(Q1::QExpr, Q2::QExpr, ::Val{C})::QExpr where {C}
    statespace_check_if(Val(C), Q1, Q2)
    new_terms = simplify_QExpr(vcat(Q1.terms, .-Q2.terms))
    return QExpr(Q1.statespace, new_terms)
end
@inline function _sub(Q1::QExpr, Q2::T, ::Val{C})::QExpr where {T<:QComposite,C}
    statespace_check_if(Val(C), Q1, Q2)
    new_terms = vcat(Q1.terms, -Q2)
    return QExpr(Q1.statespace, new_terms)
end
@inline function _sub(Q1::QExpr, N::Number, ::Val{C})::QExpr where {C}
    new_terms = vcat(Q1.terms, QAtomProduct(Q1.statespace, -Q1.statespace.c_one * N, QTerm(Q1.statespace.I_op)))
    return QExpr(Q1.statespace, new_terms)
end
@inline _sub(N::Number, Q1::QExpr, ::Val{C}) where {C} = _add(-N, Q1, _NOCHK)

#### Multiply ####################################################################
# Multiplies two QTerm’s from the same statespace. (explicit StateSpace argument)
function multiply_qterm(t1::Vector{Vector{Int}}, t2::Vector{Vector{Int}}, ss::StateSpace)
    i = 0
    results::Vector{Vector{Tuple{ComplexRational,Is}}} = Vector{Tuple{ComplexRational,Is}}[]
    for s in ss.subspaces, _ in eachindex(s.ss_inner_ind)
        i += 1
        push!(results, s.op_set.op_product(t1[i], t2[i]))
    end
    new_inds = Vector{Is}[]
    new_coeffs = ComplexRational[]
    for combo in Base.Iterators.product(results...)
        push!(new_coeffs, foldl(*, (f for (f, _) in combo); init=ComplexRational(1, 0, 1)))
        push!(new_inds, [idx for (_, idx) in combo])
    end
    return new_inds, new_coeffs
end
function multiply_qterm(t1::QTerm, t2::QTerm, ss::StateSpace)
    new_inds, new_coeffs = multiply_qterm(t1.op_indices, t2.op_indices, ss)
    return QTerm[QTerm(inds) for inds in new_inds], new_coeffs
end
function *(t1::QTerm, t2::QTerm, statespace::StateSpace)  # not user-facing
    return multiply_qterm(t1, t2, statespace)
end

# ==============================
# Main Multiplication Functions
# ==============================
# Public APIs (checked by default)
*(p1::QAtomProduct, p2::QAtomProduct)::Vector{QComposite} = _mul(p1, p2)
(*(p1::T1, p2::T2)::Vector{QComposite}) where {T1<:QComposite,T2<:QComposite} = _mul(p1, p2)
(*(p1::QSum, p2::T2)::Vector{QSum}) where {T2<:QComposite} = _mul(p1, p2)
(*(p1::QCompositeProduct, p2::T2)::Vector{QComposite}) where {T2<:QComposite} = _mul(p1, p2)
(*(p2::T2, p1::QCompositeProduct)::Vector{QComposite}) where {T2<:QComposite} = _mul(p2, p1)
*(p1::QCompositeProduct, p2::QCompositeProduct)::Vector{QComposite} = _mul(p1, p2)

*(Q1::QExpr, Q2::QExpr)::QExpr = _mul(Q1, Q2)
(*(Q1::QExpr, Q2::T)::QExpr) where {T<:QComposite} = _mul(Q1, Q2)
(*(Q1::T, Q2::QExpr)::QExpr) where {T<:QComposite} = _mul(Q1, Q2)

# Number interactions: public APIs (no check)
function *(num::Number, p1::QAtomProduct)::Vector{QComposite}
    if iszero(num)
        return QComposite[]
    else
        return [modify_coeff(p1, p1.coeff_fun * num)]
    end
end
*(p1::QAtomProduct, num::Number)::Vector{QComposite} = _mul(num, p1, _NOCHK)

function *(Q1::QExpr, num::Number)::QExpr
    if num == 0
        return QExpr(Q1.statespace, QComposite[])
    end
    terms::Vector{QComposite} = QComposite[]
    for q in Q1.terms
        append!(terms, _mul(q, num, _NOCHK))
    end
    return QExpr(Q1.statespace, terms)
end
*(num::Number, Q1::QExpr)::QExpr = Q1 * num

function *(Q1::T, num::Number)::Vector{QComposite} where {T<:QComposite}
    return [modify_expr(Q1, Q1.expr * num)]
end
(*(num::Number, Q2::T)::Vector{QComposite}) where {T<:QComposite} = Q2 * num

# Helpers (Val-driven)
@inline function _mul(Q1::QExpr, Q2::QExpr, ::Val{C})::QExpr where {C}
    statespace_check_if(Val(C), Q1, Q2)
    new_terms::Vector{QComposite} = QComposite[]
    for t1 in Q1.terms
        for t2 in Q2.terms
            append!(new_terms, _mul(t1, t2, _NOCHK))
        end
    end
    return QExpr(Q1.statespace, new_terms)
end

@inline function _mul(Q1::QExpr, Q2::T, ::Val{C})::QExpr where {T<:QComposite,C}
    statespace_check_if(Val(C), Q1, Q2)
    terms::Vector{QComposite} = QComposite[]
    for q in Q1.terms
        append!(terms, _mul(q, Q2, _NOCHK))
    end
    return QExpr(Q1.statespace, terms)
end

@inline function _mul(Q1::T, Q2::QExpr, ::Val{C})::QExpr where {T<:QComposite,C}
    statespace_check_if(Val(C), Q1, Q2)
    terms::Vector{QComposite} = QComposite[]
    for q in Q2.terms
        append!(terms, _mul(Q1, q, _NOCHK))
    end
    return QExpr(Q1.statespace, terms)
end

@inline function _mul(p1::QAtomProduct, p2::QAtomProduct, ::Val{C}) where {C}
    statespace_check_if(Val(C), p1, p2)
    return multiply_QAtomProducts(p1, p2)
end

@inline function _mul(num::Number, p1::QAtomProduct, ::Val{C}) where {C}
    if iszero(num)
        return QComposite[]
    else
        return [modify_coeff(p1, p1.coeff_fun * num)]
    end
end
@inline _mul(p1::QAtomProduct, num::Number, ::Val{C}) where {C} = _mul(num, p1, _NOCHK)

@inline function _mul(p1::T1, p2::T2, ::Val{C}) where {T1<:QComposite,T2<:QComposite,C}
    statespace_check_if(Val(C), p1, p2)
    ss = p1.statespace  # same after check
    if is_numeric(p1)
        return [modify_coeff(p2, get_coeff(p1) * get_coeff(p2))]
    end
    p1_coeff, p1_new = separate_coeff_qcomposite(p1)
    p2_coeff, p2_new = separate_coeff_qcomposite(p2)
    coeff = p1_coeff * p2_coeff
    c, t = add_QComposite_to_QCompositeProduct([p1_new], p2_new, ss)
    return [QCompositeProductCleanup(ss, c * coeff, t, Val(:nosimp))]
end

@inline function _mul(p1::QSum, p2::T2, ::Val{C})::Vector{QSum} where {T2<:QComposite,C}
    statespace_check_if(Val(C), p1, p2)
    new_expr::Vector{QComposite} = []
    for t in p1.expr
        append!(new_expr, _mul(t, p2, _NOCHK))
    end
    return [modify_expr(p1, new_expr)]
end

@inline function _mul(p1::QCompositeProduct, p2::T2, ::Val{C})::Vector{QComposite} where {T2<:QComposite,C}
    statespace_check_if(Val(C), p1, p2)
    p2_coeff, p2_new = separate_coeff_qcomposite(p2)
    coeff = p1.coeff_fun * p2_coeff
    return multiply_QCompositeProducts(coeff, p1.expr, [p2_new], Val(:nosimp))
end

@inline function _mul(p2::T2, p1::QCompositeProduct, ::Val{C})::Vector{QComposite} where {T2<:QComposite,C}
    statespace_check_if(Val(C), p2, p1)
    p2_coeff, p2_new = separate_coeff_qcomposite(p2)
    coeff = p1.coeff_fun * p2_coeff
    return multiply_QCompositeProducts(coeff, [p2_new], p1.expr, Val(:nosimp))
end

@inline function _mul(p1::QCompositeProduct, p2::QCompositeProduct, ::Val{C})::Vector{QComposite} where {C}
    statespace_check_if(Val(C), p1, p2)
    return multiply_QCompositeProducts(p1.coeff_fun * p2.coeff_fun, p1.expr, p2.expr, Val(:nosimp))
end

# ==============================
# Exponentiation
# ==============================
function ^(Q::QExpr, n::Integer)::QExpr
    if n < 0
        error("Negative exponent not defined for subtype of QComposite or QExpr.")
    elseif n == 0
        return Identity(Q.statespace)
    end
    result = Identity(Q.statespace)
    base = Q
    exp = n
    while exp > 0
        if isodd(exp)
            result = _mul(result, base, _NOCHK)
        end
        base = _mul(base, base, _NOCHK)
        exp = exp ÷ 2
    end
    return result
end

function ^(Q::T, n::Integer) where {T<:QComposite}
    if n < 0
        error("Negative exponent not defined for subtype of QComposite or QExpr.")
    elseif n == 0
        return [IdentityQAtomProduct(Q.statespace)]
    end
    result = [Identity(Q.statespace)]
    base = Q
    exp = n
    while exp > 0
        if isodd(exp)
            result = result .* base  # elementwise uses your existing scalar/vec rules
        end
        base = base * base
        exp = exp ÷ 2
    end
    return result
end

# ==============================
# Commutators
# ==============================
"""
    Commutator(Q1::Union{QExpr, S}, Q2::Union{QExpr, T}) where {S<:QComposite, T<:QComposite}

Computes the commutator [Q1, Q2] = Q1 * Q2 - Q2 * Q1.
Outermost call checks that both share the same statespace.
"""
# Public APIs
(Commutator(Q1::QExpr, Q2::QExpr)::QExpr) = _comm(Q1, Q2)
(Commutator(Q1::QExpr, Q2::T)::QExpr) where {T<:QComposite} = _comm(Q1, Q2)
(Commutator(Q1::T, Q2::QExpr)::QExpr) where {T<:QComposite} = _comm(Q1, Q2)
function Commutator(Q1::S, Q2::T)::Vector{QComposite} where {S<:QComposite,T<:QComposite}
    return _comm(Q1, Q2)
end

# Helpers
@inline function _comm(Q1::QExpr, Q2::QExpr, ::Val{C})::QExpr where {C}
    statespace_check_if(Val(C), Q1, Q2)
    return _sub(_mul(Q1, Q2, _NOCHK), _mul(Q2, Q1, _NOCHK), _NOCHK)
end
@inline function _comm(Q1::QExpr, Q2::T, ::Val{C})::QExpr where {T<:QComposite,C}
    statespace_check_if(Val(C), Q1, Q2)
    return _sub(_mul(Q1, Q2, _NOCHK), _mul(Q2, Q1, _NOCHK), _NOCHK)
end
@inline function _comm(Q1::T, Q2::QExpr, ::Val{C})::QExpr where {T<:QComposite,C}
    statespace_check_if(Val(C), Q1, Q2)
    return _sub(_mul(Q1, Q2, _NOCHK), _mul(Q2, Q1, _NOCHK), _NOCHK)
end
@inline function _comm(Q1::S, Q2::T, ::Val{C})::Vector{QComposite} where {S<:QComposite,T<:QComposite,C}
    statespace_check_if(Val(C), Q1, Q2)
    return vcat(_mul(Q1, Q2, _NOCHK), .-_mul(Q2, Q1, _NOCHK))
end

# ==============================
# Commutator sugar with Vector{QExpr}
# ==============================
function +(Q::QExpr, exprs::Vector{QExpr})::QExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return c + Q
end
function +(exprs::Vector{QExpr}, Q::QExpr)::QExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return c + Q
end
function -(Q::QExpr, exprs::Vector{QExpr})::QExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return Q - c
end
function -(exprs::Vector{QExpr}, Q::QExpr)::QExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return c - Q
end
function *(Q::QExpr, exprs::Vector{QExpr})::QExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return Q * c
end
function *(exprs::Vector{QExpr}, Q::QExpr)::QExpr
    if length(exprs) != 2
        error("Only vectors of length 2 can be interpreted as Commutator arguments.")
    end
    c = Commutator(exprs[1], exprs[2])
    return c * Q
end

# ==============================
# Identity helpers
# ==============================
function Identity(qspace::StateSpace)::QExpr
    return QExpr(qspace, QAtomProduct(qspace, qspace.c_one, QAtom[QTerm(qspace.I_op)]))
end
function Identity(qspace::StateSpace, coeff_fun::CFunction)::QExpr
    return QExpr(qspace, QAtomProduct(qspace, coeff_fun, QAtom[QTerm(qspace.I_op)]))
end
function IdentityQAtomProduct(qspace::StateSpace)::QAtomProduct
    return QAtomProduct(qspace, qspace.c_one, QAtom[QTerm(qspace.I_op)])
end
function IdentityQAtomProduct(qspace::StateSpace, coeff_fun::CFunction)::QAtomProduct
    return QAtomProduct(qspace, coeff_fun, QAtom[QTerm(qspace.I_op)])
end

# ==============================
# Daggers (single-arg; no pairwise checks)
# ==============================
function Dag(qspace::StateSpace, t::QTerm)::Vector{Tuple{QTerm,ComplexRational}}
    new_op_inds::Vector{Vector{Tuple{ComplexRational,Is}}} = []
    curr_op_inds = t.op_indices
    i = 1
    for subspace in qspace.subspaces
        for _ in 1:length(subspace.keys)
            op = curr_op_inds[i]
            push!(new_op_inds, subspace.op_set.op_dag(op))
            i += 1
        end
    end
    terms::Vector{QTerm} = QTerm[]
    coeffs::Vector{ComplexRational} = ComplexRational[]
    for combo in product(new_op_inds...)
        curr_inds::Vector{Vector{Int}} = []
        curr_coeff::ComplexRational = ComplexRational(1, 0, 1)
        for i in 1:length(combo)
            curr_coeff *= combo[i][1]
            push!(curr_inds, combo[i][2])
        end
        push!(terms, QTerm(curr_inds))
        push!(coeffs, curr_coeff)
    end
    return collect(zip(terms, coeffs))
end

function Dag(qspace::StateSpace, t::QAbstract)::Vector{Tuple{QAbstract,ComplexRational}}
    if t.operator_type.hermitian
        return Tuple{QAbstract,ComplexRational}[(t, one(ComplexRational))]
    end
    new_t = dag_copy(t)
    return Tuple{QAbstract,ComplexRational}[(new_t, one(ComplexRational))]
end

function Dag(p::QAtomProduct)::Vector{QAtomProduct}
    terms::Vector{Vector{Tuple{QAtom,ComplexRational}}} = []
    for t in reverse(p.expr)
        push!(terms, Dag(p.statespace, t))
    end
    coeff_fun = p.coeff_fun
    new_atom_products::Vector{QAtomProduct} = []
    for combo in Iterators.product(terms...)
        curr_atoms = [a[1] for a in combo]
        curr_coeff = [a[2] for a in combo]
        coeff = reduce(*, curr_coeff)
        if !iszero(coeff)
            push!(new_atom_products, QAtomProduct(p.statespace, coeff_fun * coeff, curr_atoms, p.separate_expectation_values))
        end
    end
    return new_atom_products
end

function Dag(Q::QExpr)::QExpr
    terms::Vector{QComposite} = QComposite[]
    sizehint!(terms, length(Q.terms))
    for term in Q.terms
        append!(terms, Dag(term))
    end
    return QExpr(Q.statespace, terms)
end

function Dag(t::T)::Vector{QComposite} where {T<:QComposite}
    t_new = copy(t)
    t_new.expr = Dag(t.expr)
    return QComposite[t_new]
end

function Dag(t::QMultiComposite)::Vector{QMultiComposite}
    dag_exprs::Vector{QExpr} = []
    for expr in reverse(t.expr)
        append!(dag_exprs, Dag(expr))
    end
    return [modify_expr(t, dag_exprs)]
end

adjoint(Q::QExpr) = Dag(Q)
adjoint(Q::T) where {T<:QComposite} = Dag(Q)