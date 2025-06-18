function same_term_type(t1::qAtom, t2::qAtom) # generic case 
    return false
end
function same_term_type(t1::qTerm, t2::qTerm)
    return t1.op_indices == t2.op_indices 
end
function same_term_type(t1::qAbstract, t2::qAbstract) 
    return t1.key_index == t2.key_index && t1.sub_index == t2.sub_index && t1.index_map == t2.index_map
end
function same_term_type(t1::qAtomProduct, t2::qAtomProduct)::Bool
    # check operators individually -> must all be the same!
    if length(t1.expr) != length(t2.expr) 
        same = false
    end
    for (a,b) in zip(t1.expr, t2.expr) 
        if !same_term_type(a,b)
            return false 
        end
    end
    return true
end
function same_term_type(s1::qSum, s2::qSum)::Bool
    return s1.subsystem_index == s2.subsystem_index && s1.element_indexes == s2.element_indexes
end
function same_term_type(t1::qComposite, t2::qComposite)::Bool   # term grouping either not implemented for this object or objects aren'T the same 
    return false 
end

function combine_term(t1::qAtomProduct, t2::qAtomProduct)::qAtomProduct
    return qAtomProduct(t1.statespace, simplify(t1.coeff_fun + t2.coeff_fun), t1.expr)
end
function combine_term(s1::qSum, s2::qSum)::qSum
    return qSum(s1.statespace, simplify(s1.expr + s2.expr), s1.indexes, s1.subsystem_index, s1.element_indexes, s1.neq)
end

import ..FFunctions: simplify
function simplify(p::qAtomProduct)::Vector{qComposite}
    p_new = qComposite[pad_before_qAbstracts(p)]
    did_any = true 
    while did_any
        left_terms = qComposite[] 
        for curr_p in p_new 
            append!(left_terms, qTerms2left(curr_p)) 
        end 
        did_any = false 
        p_new = qComposite[] 
        
        for i in 1:length(left_terms)
            new_term, did_it = reduce_qabstractpairs(left_terms[i])
            push!(p_new, new_term)
            if did_it
                did_any = true
            end
        end
    end
    return p_new
    #return qExpr(p.statespace, p_new)
end  
function simplify(p::qCompositeProduct)::Vector{qCompositeProduct}
    return [copy(p)]
end         
function simplify(qcomp::qMultiComposite; kwargs...)::Vector{qComposite}
    new_exprs::Vector{qExpr} = [simplify(x; kwargs...) for x in qcomp.expr]
    q = copy(qcomp)
    q.expr = new_exprs
    return [q]
end
function simplify(p::qComposite)::Vector{qComposite}
    new_p = copy(p) 
    new_p.expr = simplify(p.expr)
    return [new_p]
end 

function simplify(q::qExpr)::qExpr
    # If there are no terms, return an empty qExpr.
    if isempty(q.terms)
        return qExpr(q.statespace, qComposite[])
    end

    expr::Vector{qComposite} = qComposite[]
    
    for t in q.terms 
        append!(expr, simplify(t))
    end
    q = qExpr(q.statespace, expr)
    
    # First, sort qExpr without modifying the original.
    sorted_q = sort(q)
    sorted_terms = copy(sorted_q.terms)

    combined_terms = qComposite[]
    i = 1
    curr_term = sorted_terms[1]
    curr_i = 1
    while i < length(sorted_terms)
        # Combine adjacent like terms.
        next_term = sorted_terms[i+1]
        if same_term_type(curr_term, next_term)
            curr_term = combine_term(curr_term, next_term)
        else
            if !iszero(curr_term)
                if isa(curr_term, qSum)
                    simplified_curr_term = simplify(curr_term)
                    if !iszero(simplified_curr_term)
                        append!(combined_terms, simplified_curr_term)
                    end
                else
                    push!(combined_terms, copy(curr_term))
                end
                curr_term = next_term
                curr_i = i + 1
            end
        end
        i += 1
    end
    if !iszero(curr_term)
        if isa(curr_term, qSum)
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
function simplify(s::qSum)::Vector{qComposite}
    simplified_expr = simplify(s.expr)
    return [qSum(s.statespace, simplified_expr, s.indexes, s.subsystem_index, s.element_indexes, s.neq)]
end


function simplify(q::diff_qEQ)::diff_qEQ
    simp_rhs = simplify(q.expr)
    return diff_qEQ(q.left_hand_side, simp_rhs, q.statespace, q.braket, q.do_sigma)
end