export substitute 

function extract_qabstract(q::qExpr)::QAbstract
    if length(q) > 1 
        error("Cannot substitute composites. abstract_op must contain only an abstract operator.")
    end
    term = q.terms[1]
    if !isa(term, qAtomProduct)
        error("qExpr must contain only a qAtomProduct")
    end
    return extract_qabstract(term)
end
function extract_qabstract(term::qAtomProduct)::QAbstract
    if length(term.expr) != 1 || !isa(term.expr[1], QAbstract)
        error("abstract_op must contain exactly one QAbstract")
    end
    return term.expr[1]
end
#now same for QAtom 
function extract_qatom(q::qExpr)::QAtom
    if length(q) > 1 
        error("Cannot substitute composites. abstract_op must contain only an abstract operator.")
    end
    term = q.terms[1]
    if !isa(term, qAtomProduct)
        error("qExpr must contain only a qAtomProduct")
    end
    return extract_qatom(term)
end
function extract_qatom(term::qAtomProduct)::QAtom
    if length(term.expr) != 1 || !isa(term.expr[1], QAtom)
        error("abstract_op must contain exactly one QAbstract")
    end
    return term.expr[1]
end

# substitute abstract operator 
# input is abstract_op, replacement and target
simpleQ = Union{qExpr, qAtomProduct}
# Deals with index_map mapping one index to another ithin QAbstract
function qAtom_index_flip(q::QAtom, index_map::Vector{Tuple{Int,Int}}, statespace::StateSpace)::Vector{qAtomProduct}
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
    return [qAtomProduct(statespace, statespace.fone*c, [q]) for (q, c) in zip(qs, cs)]
end
function substitute_qAtom(abstract_op::QAbstract, replacement::QAtom, target::QTerm, statespace::StateSpace)::Vector{qAtomProduct}
    return [qAtomProduct(statespace, statespace.fone, [target])]
end
function substitute_qAtom(abstract_op::QAbstract, replacement::QAtom, target::QAbstract, statespace::StateSpace)::Vector{qAtomProduct}
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
        return [qAtomProduct(statespace, statespace.fone, [target])]
    end
end

""" 
    substitute(abstract_op::Union{qExpr, qAtomProduct, QAbstract}, replacement::Union{qExpr, qAtomProduct, QAtom}, target::Diff_qEQ, statespace::StateSpace) -> diff_q
    substitute(abstract_op::Union{qExpr, qAtomProduct, QAbstract}, replacement::Union{qExpr, qAtomProduct, QAtom}, target::qExpr) -> qExpr

Substitutes `abstract_op` with `replacement` in `target`. Keeps track of index changes due to for example `neq`. 
"""
function substitute(abstract_op::QAbstract, replacement::QAtom, target::qAtomProduct)::Vector{QComposite}
    # recursively navigate expression, and substitue
    statespace = target.statespace
    expr = target.expr
    coeff_fun = target.coeff_fun
    new_expr::qExpr = qExpr(target.statespace, substitute_qAtom(abstract_op, replacement, expr[1], statespace))
    for t in expr[2:end]
        new_terms = substitute_qAtom(abstract_op, replacement, t, statespace)
        new_new_expr = new_expr * new_terms[1]
        for t in new_terms[2:end]
            new_new_expr += new_expr + t
        end
        new_expr = new_new_expr
    end
    terms = simplify(new_expr).terms
    return [qAtomProduct(statespace, coeff_fun*t.coeff_fun, t.expr) for t in terms]  
end


function substitute(abstract_op::Union{simpleQ, QAbstract}, replacement::Union{simpleQ, QAtom}, target::qExpr)::qExpr
    if !isa(abstract_op, QAbstract)
        abstract_op = extract_qabstract(abstract_op)
    end
    if !isa(replacement, QAtom)
        replacement = extract_qatom(replacement)
    end
    return substitute(abstract_op, replacement, target)
end
function substitute(abstract_op::QAbstract, replacement::QAtom, target::qExpr)::qExpr
    # recursively navigate expression, and substitue
    new_terms = QComposite[]
    for term in target.terms
        append!(new_terms, substitute(abstract_op, replacement, term))
    end
    return qExpr(target.statespace, new_terms)
end
function substitute(a::QAbstract, r::QAtom, targ::T) where T<:QComposite
    # T<:QMultiComposite is *also* <:QComposite, 
    # so we need the QMultiComposite method to be more specific
    cp = copy(targ)
    cp.expr = substitute(a, r, targ.expr)
    return [cp]  # Tuple or Vector, depending on your convention
end

# For anything that holds *many* sub‑expressions
function substitute(a::QAbstract, r::QAtom, targ::T) where T<:QMultiComposite
    cp = copy(targ)
    cp.expr = map(x -> substitute(a, r, x), targ.expr)
    return [cp]
end

function substitute(abstract_op::Union{simpleQ, QAbstract}, replacement::Union{simpleQ, QAtom}, target::Diff_qEQ)::Diff_qEQ 
    if !isa(abstract_op, QAbstract)
        abstract_op = extract_qabstract(abstract_op)
    end
    if !isa(replacement, QAtom)
        replacement = extract_qatom(replacement)
    end
    return substitute(abstract_op, replacement, target)
end
function substitute(abstract_op::QAbstract, replacement::QAtom, target::Diff_qEQ)::Diff_qEQ
    lhs = substitute(abstract_op, replacement, target.left_hand_side)
    if length(lhs) != 1
        error("Substitution of $abstract_op with $replacement in $target did not result in a single term.")
    end
    lhs = lhs[1]
    rhs = substitute(abstract_op, replacement, target.expr)
    return Diff_qEQ(target.statespace, lhs, rhs, target.braket)
end