export substitute, Substitution, -->


# --- typed Substitution + ASCII operator ------------------------------------------
struct Substitution
    from::QAbstract
    to::QAtom
    statespace::StateSpace
end

# Primary constructor: take two QExprs, extract (QAbstract -> QAtom), check statespace
function Substitution(from_expr::QExpr, to_expr::QExpr)
    from_abs = extract_qabstract(from_expr)
    to_atom  = extract_qatom(to_expr)
    ss_from  = from_expr.statespace
    ss_to    = to_expr.statespace
    ss_from === ss_to || error("Statespace mismatch between 'from' and 'to'.")
    return Substitution(from_abs, to_atom, ss_from)
end

# ASCII-friendly constructor:  q_from --> q_to
const --> = Substitution

# Handy alias for the concrete mapping used below
const AtoP = Substitution



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
function qAtom_index_flip(q::QAtom, index_map::Vector{Tuple{SubSpaceIndex,SubSpaceIndex}}, statespace::StateSpace)::Vector{QAtomProduct}
    qs::Vector{QAtom} = [q]
    cs::Vector{ComplexRational} = [ComplexRational(1,0,1)]
    for (index1, index2) in index_map
        new_qs::Vector{QAtom} = []    
        new_cs::Vector{ComplexRational} = []
        for (qi, ci) in zip(qs, cs)
            _, new_terms, new_coeffs = term_equal_indexes( qi, index1.expanded, index2.expanded, statespace.subspaces[index1.inner])
            append!(new_qs, new_terms)
            append!(new_cs, new_coeffs*ci)
        end
        qs = new_qs
        cs = new_cs
    end
    return [QAtomProduct(statespace, statespace.c_one*c, [q]) for (q, c) in zip(qs, cs)]
end

# --- substitution on single terms (new order: target first, then Substitution) -----
# QTerm that is NOT the abstract operator: no change
function substitute_qAtom(target::QTerm, sp::AtoP)::Vector{QAtomProduct}
    ss = sp.statespace
    return [QAtomProduct(ss, ss.c_one, [target])]
end

# QAbstract: do the replacement when it matches the 'from'
function substitute_qAtom(target::QAbstract, sp::AtoP)::Vector{QAtomProduct}
    a = sp.from
    r = sp.to
    ss = sp.statespace
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

# --- top-level substitute (new order everywhere) ------------------------------
# QAtomProduct → Vector{QComposite}
function substitute(target::QAtomProduct, sp::AtoP)::Vector{QComposite}
    # sanity: statespace compatibility
    target.statespace === sp.statespace || error("Statespace mismatch between target and Substitution.")
    ss = target.statespace
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

substitution_properties_fulfilled(sub::Substitution)::Bool = substitution_properties_fulfilled(sub.from , QExpr(sub.statespace, sub.to))
# QExpr → QExpr
function substitute(target::QExpr, sp::AtoP, checks::Bool=false)::QExpr
    if checks
        substitution_properties_fulfilled(sp)
    end
    target.statespace === sp.statespace || error("Statespace mismatch between target and Substitution.")
    new_terms = QComposite[]
    for term in target.terms
        append!(new_terms, substitute(term, sp))
    end
    return QExpr(target.statespace, new_terms)
end

# Any single composite holding one expr
function substitute(targ::T, sp::AtoP) where {T<:QComposite}
    # Note: QMultiComposite is <: QComposite; a more specific method follows below.
    return [modify_expr(targ, substitute(targ.expr, sp))]
end

# Any multi-composite holding many sub-expressions
function substitute(targ::T, sp::AtoP) where {T<:QMultiComposite}
    return [modify_expr(targ, map(x -> substitute(x, sp), targ.expr))]
end

# diff_QEq → diff_QEq
function substitute(target::diff_QEq, sp::AtoP; checks::Bool=true)::diff_QEq
    if checks
        substitution_properties_fulfilled(sp)
    end
    target.statespace === sp.statespace || error("Statespace mismatch between target and Substitution.")
    lhs = substitute(target.left_hand_side, sp) 
    if length(lhs) != 1
        error("Substitution of $(sp.from) with $(sp.to) in $target did not result in a single term.")
    end
    rhs = substitute(target.expr, sp, false)
    return diff_QEq(target.statespace, lhs[1], rhs, do_braket=target.do_braket)
end

