#### First output is a vector of tuples each a coefficient and a QAtom product, the second indicates whether we changed the order, the third indicates whether we should stop
# Functions return: (new_atoms, changed, go_left)

@inline function simplify_pair(a::QAbstract, b::QAbstract, ss::StateSpace)
    # Different operator subtypes: maybe repartition
    if !same_term_type(a, b)
        if commutes_QAtom(a, b, ss) && b < a
            return [(one(ComplexRational), QAtom[b, a])], true, true
        else
            return [(one(ComplexRational), QAtom[a, b])], false, false
        end
    end

    # Same subtype ⇒ accumulation/cancellation rules
    if xor(a.dag, b.dag)
        if a.operator_type.hermitian
            exp = a.exponent + b.exponent
            if exp == 0
                return [(one(ComplexRational), QAtom[])], true, false
            else
                return [(one(ComplexRational), QAtom[modify_exp_dag(a, exp, false)])], true, false
            end
        elseif a.operator_type.unitary
            # dagger means inverse ⇒ signed exponent sum
            exp = (a.dag ? -a.exponent : a.exponent) + (b.dag ? -b.exponent : b.exponent)
            if exp == 0
                return [(one(ComplexRational), QAtom[])], true, false
            else
                dag = exp < 0
                exp = abs(exp)
                return [(one(ComplexRational), QAtom[modify_exp_dag(a, exp, dag)])], true, false
            end
        else
            # No special rule; keep order
            return [(one(ComplexRational), QAtom[a, b])], false, false
        end
    else
        # Same dag
        exp = a.exponent + b.exponent
        dag = a.dag
        if a.operator_type.hermitian
            dag = false
        elseif a.operator_type.unitary && exp < 0
            dag = !dag
            exp = -exp
        end
        if exp == 0
            return [(one(ComplexRational), QAtom[])], true, false
        else
            return [(one(ComplexRational), QAtom[modify_exp_dag(a, exp, dag)])], true, false
        end
    end
end


# QTerm × QTerm  → multiply (may branch)
@inline function simplify_pair(x::QTerm, y::QTerm, ss::StateSpace)::Tuple{Vector{Tuple{ComplexRational, Vector{QAtom}}}, Bool, Bool}
    if x.time_index == y.time_index 
        Ts, Cs = multiply_qterm(x, y, ss)
        out = Vector{Tuple{ComplexRational, Vector{QAtom}}}()
        @inbounds for (t, c) in zip(Ts, Cs)
            if !iszero(c) 
                if isnumeric(t, ss) 
                    push!(out, (c, Vector{QAtom}([])))
                else
                    push!(out, (c, Vector{QAtom}([t])))
                end
            end
        end
        return out, true, true
    end
    return [(one(ComplexRational), QAtom[x, y])], false, false
end

@inline function simplify_pair(x::QAbstract, y::QTerm, ss::StateSpace)::Tuple{Vector{Tuple{ComplexRational, Vector{QAtom}}}, Bool, Bool}
    if commutes_QAtom(x, y, ss)
        return [(one(ComplexRational), QAtom[y, x])], true, true
    else
        return [(one(ComplexRational), QAtom[x, y])], false, false
    end
end
@inline function simplify_pair(x::QTerm, y::QAbstract, ss::StateSpace)::Tuple{Vector{Tuple{ComplexRational, Vector{QAtom}}}, Bool, Bool}
    return [(one(ComplexRational), QAtom[x, y])], false, false
end

"""
    add_QAtom_to_QAtomProduct(terms, a, ss) -> Vector{(coeff, Vector{QAtom})}

Append `a` to `terms`, then bubble it left:
- keep going while swaps/mults/cancellations happen,
- always step back one after a change so new junctions can simplify.
Branches are preserved (sum of products).
"""
function add_QAtom_to_QAtomProduct(terms::Vector{QAtom}, a::QAtom, ss::StateSpace)::Vector{Tuple{ComplexRational, Vector{QAtom}}}
    seed = copy(terms)
    push!(seed, a)

    # index of the rightmost pair (i,i+1) to examine; 0 means "no pair"
    start_i = length(seed) > 1 ? length(seed) - 1 : 0
    states = Tuple{ComplexRational,Vector{QAtom},Int}[(one(ComplexRational), seed, start_i)]
    results = Tuple{ComplexRational,Vector{QAtom}}[]

    while !isempty(states)
        next = Tuple{ComplexRational,Vector{QAtom},Int}[]
        for (c, t, i) in states
            if i > 0
                @inbounds outs, changed, go_left = simplify_pair(t[i], t[i+1], ss)  # no where_acting
                if changed
                    @inbounds for (dc, pair) in outs
                        nt = Vector{QAtom}()
                        append!(nt, t[1:i-1])
                        append!(nt, pair)          # pair may be length 0, 1, or 2
                        append!(nt, t[i+2:end])

                        if go_left 
                            new_index = i-1
                            if new_index < 1 
                                new_index = 2
                            end
                        else
                            new_index = length(pair) == 2 ? i+1 : i 
                        end
                        if new_index == 0 # can't go further left 
                            new_index = 2 # hence try going right
                        end

                        # can'T we just add another case here for when it gets smaller than 1 through a go left operation?
                        if new_index >= length(nt) 
                            push!(results, (c*dc, nt))
                        else
                            push!(next, (c*dc, nt, new_index))
                        end
                    end
                else
                    # no local change → move right
                    new_i = i+1
                    if new_i >= length(t)   
                        push!(results, (c, t))  # i did change this block!
                    else
                        push!(next, (c, t, new_i))
                    end
                end
            else
                push!(results, (c, t))
            end
        end
        states = next
    end
    # Outer vector = sum, inner Vector{QAtom} = product
    return results
end


function multiply_QAtomProducts_terms(p1::Vector{QAtom}, p2::Vector{QAtom}, statespace::StateSpace)
    # start from simplified p1 (sum of products)
    states = [(ComplexRational(1,0,1), p1) ] # (coeff, Vector{QAtom})

    # insert each atom from p2 into every current branch
    for a in p2
        new_states = Tuple{ComplexRational,Vector{QAtom}}[]
        for (c, t) in states
            for (dc, nt) in add_QAtom_to_QAtomProduct(t, a, statespace)
                push!(new_states, (c*dc, nt))
            end
        end
        states = new_states
    end

    return states
end

"""
    multiply_QAtomProducts(p1::QAtomProduct, p2::QAtomProduct) → Vector{QComposite}

Wraps `multiply_QAtomProducts_terms` into your `QComposite` type.
Assumes scalar-like `coeff_fun` fields multiply.
"""
function multiply_QAtomProducts(p1::QAtomProduct, p2::QAtomProduct)::Vector{QComposite}
    ss = p1.statespace
    if p1.separate_expectation_values != p2.separate_expectation_values
        error("Cannot multiply QAtomProducts with different `separate_expectation_values`")
    end
    new_coeff_fun = p1.coeff_fun * p2.coeff_fun
    if p1.separate_expectation_values   # don't simplify the term
        return [ QAtomProduct(ss, new_coeff_fun, vcat(p1.expr, p2.expr), p1.separate_expectation_values) ]
    end
    termsums = multiply_QAtomProducts_terms(p1.expr, p2.expr, ss)
    return [ QAtomProduct(ss, c * new_coeff_fun, t, p1.separate_expectation_values) for (c, t) in termsums ]
end
