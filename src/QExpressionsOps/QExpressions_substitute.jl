export substitute 

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

# substitute abstract operator 
# input is abstract_op, replacement and target
simpleQ = Union{QExpr, QAtomProduct}
# Deals with index_map mapping one index to another ithin QAbstract
function qAtom_index_flip(q::QAtom, index_map::Vector{Tuple{Int,Int}}, statespace::StateSpace)::Vector{QAtomProduct}
    qs::Vector{QAtom} = [q]
    cs::Vector{ComplexRational} = [ComplexRational(1,0,1)]
    for (index1, index2) in index_map
        new_qs::Vector{QAtom} = []    
        new_cs::Vector{ComplexRational} = []
        for (qi, ci) in zip(qs, cs)
            _, new_terms, new_coeffs = term_equal_indexes(qi, index1, index2, statespace)
            append!(new_qs, new_terms)
            append!(new_cs, new_coeffs*ci)
        end
        qs = new_qs
        cs = new_cs
    end
    return [QAtomProduct(statespace, statespace.fone*c, [q]) for (q, c) in zip(qs, cs)]
end
function substitute_qAtom(abstract_op::QAbstract, replacement::QAtom, target::QTerm, statespace::StateSpace)::Vector{QAtomProduct}
    return [QAtomProduct(statespace, statespace.fone, [target])]
end
function substitute_qAtom(abstract_op::QAbstract, replacement::QAtom, target::QAbstract, statespace::StateSpace)::Vector{QAtomProduct}
    # check if its the same QAbstract operator 
    if target.key_index == abstract_op.key_index && target.sub_index == abstract_op.sub_index 
        qs = qAtom_index_flip(replacement, target.index_map, statespace)

        if target.exponent != 1
            qs = qs.^target.exponent
        end
        if target.dag
            qs = qs'
        end
        return qs
    else
        return [QAtomProduct(statespace, statespace.fone, [target])]
    end
end

""" 
    substitute(abstract_op::Union{QExpr, QAtomProduct, QAbstract}, replacement::Union{QExpr, QAtomProduct, QAtom}, target::diff_QEq, statespace::StateSpace) -> diff_q
    substitute(abstract_op::Union{QExpr, QAtomProduct, QAbstract}, replacement::Union{QExpr, QAtomProduct, QAtom}, target::QExpr) -> QExpr

Substitutes `abstract_op` with `replacement` in `target`. Keeps track of index changes due to for example `neq`. 
"""
function substitute(abstract_op::QAbstract, replacement::QAtom, target::QAtomProduct)::Vector{QComposite}
    # recursively navigate expression, and substitue
    statespace = target.statespace
    expr = target.expr
    coeff_fun = target.coeff_fun
    new_expr::QExpr = QExpr(target.statespace, substitute_qAtom(abstract_op, replacement, expr[1], statespace))
    for t in expr[2:end]
        new_terms = substitute_qAtom(abstract_op, replacement, t, statespace)
        new_new_expr = new_expr * new_terms[1]
        for t in new_terms[2:end]
            new_new_expr += new_expr + t
        end
        new_expr = new_new_expr
    end
    terms = simplify(new_expr).terms
    return [QAtomProduct(statespace, coeff_fun*t.coeff_fun, t.expr) for t in terms]  
end


function substitute(abstract_op::Union{simpleQ, QAbstract}, replacement::Union{simpleQ, QAtom}, target::QExpr)::QExpr
    if !isa(abstract_op, QAbstract)
        abstract_op = extract_qabstract(abstract_op)
    end
    if !isa(replacement, QAtom)
        replacement = extract_qatom(replacement)
    end
    return substitute(abstract_op, replacement, target)
end
function substitute(abstract_op::QAbstract, replacement::QAtom, target::QExpr)::QExpr
    # recursively navigate expression, and substitue
    new_terms = QComposite[]
    for term in target.terms
        append!(new_terms, substitute(abstract_op, replacement, term))
    end
    return QExpr(target.statespace, new_terms)
end
function substitute(a::QAbstract, r::QAtom, targ::T) where T<:QComposite
    # T<:QMultiComposite is *also* <:QComposite, 
    # so we need the QMultiComposite method to be more specific
    return [modify_expr(targ, substitute(a, r, targ.expr))]
end

# For anything that holds *many* sub‑expressions
function substitute(a::QAbstract, r::QAtom, targ::T) where T<:QMultiComposite
    return [modify_expr(targ, map(x -> substitute(a, r, x), targ.expr))]
end

function substitute(abstract_op::Union{simpleQ, QAbstract}, replacement::Union{simpleQ, QAtom}, target::diff_QEq)::diff_QEq 
    if !isa(abstract_op, QAbstract)
        abstract_op = extract_qabstract(abstract_op)
    end
    if !isa(replacement, QAtom)
        replacement = extract_qatom(replacement)
    end
    return substitute(abstract_op, replacement, target)
end
function substitute(abstract_op::QAbstract, replacement::QAtom, target::diff_QEq)::diff_QEq
    lhs = substitute(abstract_op, replacement, target.left_hand_side)
    if length(lhs) != 1
        error("Substitution of $abstract_op with $replacement in $target did not result in a single term.")
    end
    lhs = lhs[1]
    rhs = substitute(abstract_op, replacement, target.expr)
    return diff_QEq(target.statespace, lhs, rhs, target.braket)
end