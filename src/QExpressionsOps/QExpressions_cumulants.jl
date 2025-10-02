using ..Cumulants: ReducedIndexedCumulant
export Cumulant

struct QCumulant <: QComposite
    qspace::QSpace
    coeff_fun::CFunction
    atom::QTerm
    expr::QExpr
    order::Int
    where_acting::Vector{Int}
end

function modify_expr(q::QCumulant, expr::QExpr)::Vector{QComposite}
    return QComposite[QCumulant(q.qspace, q.coeff_fun, q.atom, expr, q.order, q.where_acting)]
end

function modify_coeff(q::QCumulant, coeff::CFunction)::QCumulant
    old = q.coeff_fun
    if coeff === old
        return q
    end
    scale = iszero(old) ? coeff : coeff / old
    return QCumulant(q.qspace, coeff, q.atom, multiply_coeff(q.expr, scale), q.order, q.where_acting)
end

each_term(q::QCumulant) = q.expr.terms
each_coeff(q::QCumulant)::Vector{CFunction} = [q.coeff_fun; each_coeff(q.expr)]
get_coeff(q::QCumulant) = q.coeff_fun

function replace_indexes(I_op::Vector{Vector{Int}}, curr_op_indexes::Vector{Vector{Int}}, indexes::Vector{Int})
    new_op = copy(I_op)
    for ind in indexes
        new_op[ind] = copy(curr_op_indexes[ind])
    end
    return new_op
end

copy_qterm(term::QTerm) = QTerm(term.op_indices, term.time_index)

function _build_approximation(atom::QTerm, qspace::QSpace, where_acting::Vector{Int}, red_cum::ReducedIndexedCumulant)
    curr_op_indexes = atom.op_indices
    I_op = qspace.I_op
    approx = QComposite[]
    for app in red_cum.approximation
        curr_atoms = QTerm[QTerm(replace_indexes(I_op, curr_op_indexes, where_acting[ind]), atom.time_index) for ind in app.indices]
        term_coeff = qspace.c_one * app.coeff
        prod = QAtomProduct(qspace, term_coeff, curr_atoms, true)
        push!(approx, prod)
    end
    return QExpr(qspace, approx)
end

function Cumulant(qprod::QAtomProduct)::QCumulant
    length(qprod.expr) == 1 || error("Cumulant requires a single QTerm inside the QAtomProduct.")
    atom = qprod.expr[1]
    atom isa QTerm || error("Cumulant only defined for QAtomProduct containing a QTerm; got $(typeof(atom)).")
    qspace = qprod.qspace
    where_acting = where_acting_index(atom, qspace)
    order = length(where_acting) # don'T use order here, since where_acting is also needed 
    red_cum = qspace.cumulant_cache(order)
    return _Cumulant(qprod.coeff_fun, atom, qspace, order, where_acting, red_cum)
end

function Cumulant(qexpr::QExpr)::QCumulant
    length(qexpr.terms) == 1 || error("Cumulant expects a QExpr with exactly one term; got $(length(qexpr.terms)).")
    term = qexpr.terms[1]
    term isa QAtomProduct || error("Cumulant only defined for QExpr whose single term is a QAtomProduct; got $(typeof(term)).")
    return Cumulant(term)
end

function Cumulant(::QComposite)
    error("Cumulant only defined for QAtomProduct (with a single QTerm) or QExpr containing exactly one such product.")
end

function Cumulant(::QAtom)
    error("Cumulant only defined for QAtomProduct (with a single QTerm) or QExpr containing exactly one such product.")
end

@inline function _Cumulant(coeff_fun::CFunction, atom::QTerm, qspace::QSpace, order::Int, where_acting::Vector{Int}, red_cum::ReducedIndexedCumulant)::QCumulant
    approx_expr = _build_approximation(atom, qspace, where_acting, red_cum)
    return QCumulant(qspace, coeff_fun, copy_qterm(atom), approx_expr, order, copy(where_acting))
end
