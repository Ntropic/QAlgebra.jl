function reduce_qabstract_pair(a::QAbstract, b::QAbstract, statespace::StateSpace)::Tuple{Vector{QAbstract}, Bool}
    # must be exactly the same operator subtype
    if !same_term_type(a, b)
        if commutes_QAbstract(a, b, statespace)
            if b < a
                return QAbstract[b, a], false
            end
        end
        return QAbstract[a,b], false
    end

    # try hermitian/unitary or plain accumulation
    if xor(a.dag, b.dag)
        if a.operator_type.hermitian
            t = copy(a)
            t.dag = false
            t.exponent += b.exponent
            if t.exponent == 0
                return QAbstract[], true
            end
            return QAbstract[t], true
        elseif a.operator_type.unitary
            t = copy(a)
            t.dag = false
            # note: exponent^dag flips sign
            t.exponent = a.exponent * (-1)^Int(a.dag) + b.exponent * (-1)^Int(b.dag)
            if t.exponent == 0
                return QAbstract[], true
            end
            return QAbstract[t], true
        end
        return QAbstract[a,b], false
    else #if !xor(a.dag, b.dag)  # bath have same dag
        t = copy(a)
        t.dag = a.dag
        t.exponent += b.exponent
        if a.operator_type.hermitian
            t.dag = false 
        elseif a.operator_type.unitary
            t.dag = false
            t.exponent = .-t.exponent
        end
        if t.exponent == 0
            return QAbstract[], true
        end
        return QAbstract[t], true
    end

    return QAbstract[a,b], false
end

#### Main simplify step function #################################################
@inline function simplify_pair(x::QAbstract, y::QAbstract, ss::StateSpace)
    new_abs, did = reduce_qabstract_pair(x, y, ss)
    if did
        out = Vector{Tuple{ComplexRational, Vector{QAtom}}}()
        for t in new_abs
            push!(out, (one(ComplexRational), Vector{QAtom}([t])))
        end
        return out, true
    end
    return Tuple{ComplexRational, Vector{QAtom}}[(one(ComplexRational), Vector{QAtom}([x,y]))], false
end

@inline function simplify_pair(x::QAbstract, y::QTerm, ss::StateSpace)
    commuting_qterm = copy(ss.neutral_op)
    non_commuting_qterm = copy(ss.neutral_op)
    c_abstract = where_acting(x, ss)
    curr_combo = Vector{QAtom}()
    for (i, (b, q_ind)) in enumerate(zip(c_abstract, y.op_indices))
        if b
            non_commuting_qterm[i] = q_ind
        else
            commuting_qterm[i] = q_ind
        end
    end
    if !is_numeric(commuting_qterm, ss)
        push!(curr_combo, QTerm(commuting_qterm))
    else
        return [(one(ComplexRational), [x, y])], false
    end
    push!(curr_combo, x)
    if !is_numeric(non_commuting_qterm, ss)
        push!(curr_combo, QTerm(non_commuting_qterm))
    end
    return Tuple{ComplexRational, Vector{QAtom}}[(one(ComplexRational), curr_combo)], true
end

@inline function simplify_pair(x::QTerm, y::QAbstract, ss::StateSpace)
    # just reverse order, logic is symmetric
    return [(one(ComplexRational), Vector{QAtom}([x, y]))], false
end

@inline function simplify_pair(x::QTerm, y::QTerm, ss::StateSpace)
    Ts, Cs = multiply_qterm(x, y, ss)
    out = Vector{Tuple{ComplexRational, Vector{QAtom}}}()
    for (t, c) in zip(Ts, Cs)
        if !iszero(c) 
            if is_numeric(t, ss) 
                push!(out, (c, Vector{QAtom}([])))
            else
                push!(out, (c, Vector{QAtom}([t])))
            end
        end
    end
    return out, true
end


function simplify_pairs(expr::Vector{QAtom}, ss::StateSpace)::Tuple{Vector{Tuple{ComplexRational, Vector{QAtom}}}, Bool}
    # We're going to track triples (coeff, terms, i)
    # where i is the next index at which to apply simplify_pair.
    results = Tuple{ComplexRational, Vector{QAtom}, Int}[(one(ComplexRational), copy(expr), length(expr)-1)]
    changed = false
    # Keep going until every triple has i == 0
    while any(triple -> triple[3] > 0, results)
        new_results = Vector{Tuple{ComplexRational, Vector{QAtom}, Int}}()
        for (c, terms, i) in results
            if i > 0
                # try a single-pair rewrite at position i
                x, y = terms[i], terms[i+1]
                outs, did = simplify_pair(x, y, ss)
                
                if did
                    changed = true
                    # for each rewrite, build the new terms and decrement i
                    for (dc, pair) in outs
                        new_terms = Vector{QAtom}()
                        append!(new_terms, terms[1:i-1])
                        append!(new_terms, pair)
                        append!(new_terms, terms[i+2:end])
                        if length(pair) == 0
                            push!(new_results, (c * dc, new_terms, i-2))
                        else
                            push!(new_results, (c * dc, new_terms, i-1))
                        end
                    end
                else
                    # no change here, just move on
                    push!(new_results, (c, terms, i-1))
                end
            else
                # already at i==0, carry forward unchanged
                push!(new_results, (c, terms, 0))
            end
        end
        
        results = new_results
    end
    
    # strip off the index and return just (coeff, terms)
    return [(c, terms) for (c, terms, _) in results], changed
end
