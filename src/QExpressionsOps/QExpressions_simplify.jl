function same_term_type(t1::QAtom, t2::QAtom) # generic case 
    return false
end
function same_term_type(t1::QTerm, t2::QTerm)
    return t1.op_indices == t2.op_indices 
end
function same_term_type(t1::QAbstract, t2::QAbstract) 
    return t1.key_index == t2.key_index && t1.sub_index == t2.sub_index && t1.index_map == t2.index_map
end
function same_term_type(t1::QAtomProduct, t2::QAtomProduct)::Bool
    # check operators individually -> must all be the same!
    if length(t1.expr) != length(t2.expr) || t1.separate_expectation_values == t2.separate_expectation_values 
        return false
    end
    for (a,b) in zip(t1.expr, t2.expr) 
        if !same_term_type(a,b)
            return false 
        end
    end
    return true
end
function same_term_type(s1::QSum, s2::QSum)::Bool
    return s1.subsystem_index == s2.subsystem_index && s1.element_indexes == s2.element_indexes
end
function same_term_type(t1::S, t2::T)::Bool where {S<:QComposite, T<:QComposite}   # term grouping either not implemented for this object or objects aren'T the same 
    return false 
end

function combine_term(t1::QAtomProduct, t2::QAtomProduct)::QAtomProduct
    return QAtomProduct(t1.statespace, FFunctions.simplify(t1.coeff_fun + t2.coeff_fun), t1.expr)
end
function combine_term(s1::QSum, s2::QSum)::QSum
    return QSum(s1.statespace, simplify(s1.expr + s2.expr), s1.indexes, s1.subsystem_index, s1.element_indexes, s1.neq)
end

function simplifyqAtomProduct(p::QAtomProduct)::Vector{QComposite}
    # — pre‑allocate two empty buffers of the correct element‑type —
    buf1 = Vector{Tuple{ComplexRational,Vector{QAtom}}}()
    buf2 = Vector{Tuple{ComplexRational,Vector{QAtom}}}()
    current, nextbuf = buf1, buf2

    # seed the first buffer
    empty!(current)
    push!(current, (one(ComplexRational), copy(p.expr)))

    if !p.separate_expectation_values 
        while true
            empty!(nextbuf)
            did_any = false

            for (coeff, terms) in current
                outs, changed = simplify_pairs(terms, p.statespace)
                
                did_any |= changed
                for (dc, t) in outs
                    if !iszero(dc)
                        push!(nextbuf, (coeff * dc, t))
                    end
                end
            end

            # if nothing changed, we're done
            if !did_any
                break
            end

            # swap buffers for the next iteration
            current, nextbuf = nextbuf, current
        end
    end

    # wrap the final products
    return [ QAtomProduct(p.statespace, c * p.coeff_fun, t) for (c,t) in current ]
end


"""
    simplify(q) -> simplified

Simplifies `qExpr`-based symbolic quantum expressions by recursively reducing internal structures:

- `QAtomProduct`: Applies pairwise simplifications repeatedly and merges results into canonical `QAtomProduct`s.
- `QMultiComposite`: Simplifies each expression element-wise.
- `QComposite`: Simplifies its internal expression and returns a new `QComposite`.
- `qExpr`: Flattens and simplifies terms, then combines like terms where possible.
- `QSum`: Simplifies its expression array and returns a new `QSum`.
- `diff_QEq`: Replaces its right-hand side with a simplified version.

Returns either a single simplified object or a list of canonical components depending on input type.
"""
function simplify(q::QAtomProduct)::Vector{QAtomProduct}
    return [q]
end
function simplify(qcomp::T)::Vector{T} where T <: QMultiComposite
    new_exprs::Vector{qExpr} = [simplify(x) for x in qcomp.expr]
    q = copy(qcomp)
    q.expr = new_exprs
    return [q]
end
function simplify(p::T)::Vector{T} where T <: QComposite
    new_p = copy(p) 
    new_p.expr = simplify(p.expr)
    return [new_p]
end 

function simplify(q::qExpr)::qExpr
    # If there are no terms, return an empty qExpr.
    if isempty(q.terms)
        return qExpr(q.statespace, QComposite[])
    end

    expr::Vector{QComposite} = QComposite[]
    
    for t in q.terms 
        append!(expr, simplify(t))
    end
    q = qExpr(q.statespace, expr)
    
    # First, sort qExpr without modifying the original.
    sorted_terms = _sort(expr)

    combined_terms = QComposite[]
    i = 1
    curr_term = sorted_terms[1]
    while i < length(sorted_terms)
        # Combine adjacent like terms.
        next_term = sorted_terms[i+1]
        if same_term_type(curr_term, next_term)
            curr_term = combine_term(curr_term, next_term)
        else
            if !iszero(curr_term)
                push!(combined_terms, copy(curr_term))
            end
            curr_term = next_term
        end
        i += 1
    end
    if !iszero(curr_term)
        if isa(curr_term, QSum)
            simplified_curr_term = simplify(curr_term)
            if !iszero(simplified_curr_term)
                append!(combined_terms, simplified_curr_term)
            end
        else
            push!(combined_terms, copy(curr_term))
        end
    end
    return qExpr(q.statespace, combined_terms)
end
function simplify(q::diff_QEq)::diff_QEq
    simp_rhs = simplify(q.expr)
    return diff_QEq(q.left_hand_side, simp_rhs, q.statespace, q.braket)
end