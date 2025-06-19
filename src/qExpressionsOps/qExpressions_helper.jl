function reduce_qabstract_pair(a::qAbstract, b::qAbstract)::Tuple{Vector{qAbstract}, Bool}
    # must be exactly the same operator subtype
    if !same_term_types(a, b)
        return qAbstract[a,b], false
    end

    # try hermitian/unitary or plain accumulation
    if xor(a.dag, b.dag)
        if a.operator_type.hermitian
            t = copy(a)
            t.dag = false
            t.exponent += b.exponent
            return qAbstract[t], true
        elseif a.operator_type.unitary
            t = copy(a)
            t.dag = false
            # note: exponent^dag flips sign
            t.exponent = a.exponent * (-1)^Int(a.dag) + b.exponent * (-1)^Int(b.dag)
            return qAbstract[t], true
        end
    elseif !(a.dag || b.dag)   # neither is dag
        t = copy(a)
        t.dag = false
        t.exponent += b.exponent
        return qAbstract[t], true
    end

    return qAbstract[a,b], false
end

#### Main simplify step function #################################################
function simplify_pair(x::qAtom, y::qAtom, ss::StateSpace)::Tuple{Vector{Tuple{ComplexRational, Vector{qAtom}}}, Bool}
    if isa(x, qAbstract) 
        if isa(y, qAbstract) # 1) two abstracts → try reduction
            new_abs, did = reduce_qabstract_pair(x, y)
            if did
                out = Vector{Tuple{ComplexRational, Vector{qAtom}}}()
                for t in new_abs
                    push!(out, (one(ComplexRational), Vector{qAtom}([t])))
                end
                return out, true
            end
        else  # 1) qAbstract * qTerm → try (partial) commutation
            commuting_qterm = copy(ss.neutral_op)
            non_commuting_qterm = copy(ss.neutral_op)
            c_abstract = where_acting(x, ss)
            curr_combo = Vector{qAtom}()
            for (i, (b, q_ind)) in enumerate(zip(c_abstract, y.op_indices))
                if b
                    non_commuting_qterm[i] = q_ind
                else
                    commuting_qterm[i] = q_ind
                end
            end
            if !is_numeric(commuting_qterm, ss)
                push!(curr_combo, qTerm(commuting_qterm))
            else
                return [(one(ComplexRational), [x, y])], false
            end
            push!(curr_combo, x)
            if !is_numeric(non_commuting_qterm, ss)
                push!(curr_combo, qTerm(non_commuting_qterm))
            end
            return [(one(ComplexRational), curr_combo)], true
        end 
    else
        if isa(y, qTerm) # 3) two qTerms → multiply/unify
            Ts, Cs = multiply_qterm(x, y, ss)
            out = Vector{Tuple{ComplexRational, Vector{qAtom}}}()
            for (t, c) in zip(Ts, Cs)
                if !iszero(c) 
                    if is_numeric(t, ss) 
                        push!(out, (c, Vector{qAtom}([])))
                    else
                        push!(out, (c, Vector{qAtom}([t])))
                    end
                end
            end
            return out, true
        else  # 4) nothing changed
            return [(one(ComplexRational), Vector{qAtom}([x, y]))], false
        end
    end 
end

function simplify_pairs(expr::Vector{qAtom}, ss::StateSpace)::Tuple{Vector{Tuple{ComplexRational, Vector{qAtom}}}, Bool}
    results = [(one(ComplexRational), copy(expr))]
    changed = false
    for i in (length(expr)-1):-1:1
        new_results = Vector{Tuple{ComplexRational, Vector{qAtom}}}()
        for (c, terms) in results
            x, y = terms[i], terms[i+1]
            outs, did = simplify_pair(x, y, ss)
            if did
                changed = true
                for (dc, pair) in outs
                    new_expr = Vector{qAtom}()
                    append!(new_expr, terms[1:i-1])
                    append!(new_expr, pair)
                    append!(new_expr, terms[i+2:end])
                    push!(new_results, (c * dc, new_expr))
                end
            else
                push!(new_results, (c, terms))
            end
        end
        results = new_results
    end
    return results, changed 
end