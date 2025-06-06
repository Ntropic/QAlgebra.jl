function pad_before_qAbstracts(p::qAtomProduct)::qAtomProduct
    ss = p.statespace
    neutral = qTerm(ss.neutral_op)

    # Find the last qAtom position
    last_atom_index = findlast(t -> isa(t, qTerm), p.expr)

    # If there are no qAtoms, no padding is necessary
    if last_atom_index === nothing
        return copy(p)
    end

    terms = copy(p.expr)
    # Iterate in reverse to keep indices stable while inserting
    for i in last_atom_index-1:-1:2
        both_qAbstract = isa(terms[i], qAbstract) && isa(terms[i-1], qAbstract)
        if both_qAbstract
            insert!(terms, i, copy(neutral))
        end
    end
    # Special case: if the first element is a qAbstract
    if isa(terms[1], qAbstract)
        insert!(terms, 1, copy(neutral))
    end
    return qAtomProduct(ss, copy(p.coeff_fun), terms)
end

function qTerms2left(p::qAtomProduct)::Vector{qAtomProduct}
    statespace = p.statespace
    coeff_fun = p.coeff_fun
    terms = copy(p.expr)
    last_atom_index = findlast(t -> isa(t, qTerm), terms)
    all_terms::Vector{Vector{qAtom}} = [copy(p.expr)]
    all_coeffs::Vector{ComplexRational} = [one(ComplexRational)]
    # for storing intermediate creations 
    new_all_terms::Vector{Vector{qAtom}} = []
    new_all_coeffs::Vector{ComplexRational} = []

    # Find the last qAtom position
    last_atom_index = findlast(t -> isa(t, qTerm), terms)
    for i in last_atom_index:-2:3
        for k in 1:length(all_terms)
            terms = all_terms[k]
            coeff = all_coeffs[k]
            curr_qterm = terms[i]
            curr_qabstract = terms[i-1]
            c_abstract = where_acting(curr_qabstract, statespace) # boolean vector
            c_term = where_acting(curr_qterm, statespace)  # boolean vector
            if !any(c_term) 
                deleteat!(terms, i)
                push!(new_all_terms, terms)
                push!(new_all_coeffs, copy(coeff))
            elseif all([nand(a,b) for (a,b) in zip(c_abstract, c_term)])   # commutes 
                prev_qterm = term[i-2]
                new_terms, new_coeffs = multiply_qterm(prev_qterm, curr_qterm, statespace)
                # construct new terms 
                deleteat!(terms, i)
                for (t,c) in zip(new_terms, new_coeffs)
                    new_term = copy(terms)
                    new_term[i-2] = t 
                    push!(new_all_terms, new_term)
                    push!(new_all_coeffs, coeff*c)
                end
            else 
                # split curr_qterm into 2, one that commutes with curr_qabstract and one that doesn't
                commuting_qterm = copy(statespace.neutral_op)
                non_commuting_qterm = copy(statespace.neutral_op)
                # add at indexes that aren't in c_abstract
                for (i, b, q_ind) in enumerate(c_abstract, curr_qterm.op_indices)
                    if b
                        non_commuting_qterm[i] = q_ind
                    else
                        commuting_qterm[i] = q_ind
                    end
                end
                prev_qterm = term[i-2]
                new_terms, new_coeffs = multiply_qterm(qTerm(prev_qterm), qTerm(commuting_qterm), statespace)
                terms[i] = qTerm(non_commuting_qterm)
                for (t, c) in zip(new_terms, new_coeffs)
                    new_term = copy(terms)
                    new_terms[i-2] = t 
                    push!(new_all_terms, new_term)
                    push!(new_all_coeffs, coeff*c)
                end
            end
        end
        all_terms = copy(new_all_terms)
        all_coeffs = copy(new_all_coeffs)
        new_terms = Vector{qAtom}[]
        new_coeffs = Vector{ComplexRational}[]
    end
    if length(p.expr) > 0
        for t in all_terms
            if is_numeric(t[1], statespace)
                deleteat!(t, 1)
            end
        end
    end
    # create qAtomProduct for each 
    return qAtomProduct[qAtomProduct(statespace, c, t) for (c,t) in zip(all_coeffs, all_terms)]
end

function reduce_qabstractpairs(p::qAtomProduct)::Tuple{qAtomProduct, Bool}
    # remove pairs of qAbstract if possible 
    ss = p.statespace
    terms = copy(p.expr)
    i = 1
    did_any = false
    while i < length(terms)
        did_this = false
        if is_qAbstract(terms[i]) && is_qAbstract(terms[i+1])
            # check if they are of the same subtypes
            # could probably reduce this to same_term_type function call!
            if same_term_types(terms[i], terms[i+1]) 
                if xor(terms[i].dag, terms[i+1].dag) 
                    if terms[i].operator_type.hermitian 
                        new_term = copy(terms[i])
                        new_term.dag = false 
                        new_term.exponent += terms[i+1].exponent
                        did_this = true
                    elseif terms[i].operator_type.unitary
                        new_term = copy(terms[i])
                        new_term.dag = false 
                        new_term.exponent = terms[i].exponent*(-1)^terms[i].dag + terms[i+1].exponent*(-1)^terms[i+1].dag
                        did_this = true
                    end
                elseif nor(terms[i].dag, terms[i+1].dag)  # neither dag 
                    new_term = copy(terms[i])
                    new_term.dag = false 
                    new_term.exponent += terms[i+1].exponent
                    did_this = true
                end
            end
            if did_this
                deleteat!(terms, i+1)
                terms[i] = new_term
                did_any = true
            else
                i += 1
            end
        end
    end
    return qAtomProduct(ss, p.coeff_fun, terms), did_any
end