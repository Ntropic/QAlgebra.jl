# Helper Functions for QComposite Constructions ==> To keep all coefficients in front  
function separate_coeff_qcomposite(q::T)::Tuple{CFunction, QComposite} where {T <: QComposite}
    return q.coeff_fun, modify_coeff(q, q.qspace.c_one)
end
function separate_coeff_qcomposite(q::QSum)::Tuple{CFunction, QComposite} 
    return q.qspace.c_one, q 
end
function separate_coeff_qcomposites(qs::Vector{QComposite}, qspace::QSpace)::Tuple{CFunction, Vector{QComposite}} 
    c_one = qspace.c_one 
    coeff_fun = c_one  
    new_vector::Vector{QComposite} = QComposite[]
    sizehint!(new_vector, length(qs))
    for q in qs 
        if isa(q, QSum)
            push!(new_vector, q)
        else
            coeff_fun *= q.coeff_fun 
            push!(new_vector, modify_coeff(q, c_one))
        end
    end
    return coeff_fun, new_vector 
end

# split replacements into (coeff, stripped factors) while dropping numeric identities
const CompositeBranch = Tuple{CFunction, Vector{QComposite}}

function _make_branch(replacements::AbstractVector{<:QComposite}, qspace::QSpace)::CompositeBranch
    coeff = qspace.c_one
    factors = QComposite[]
    sizehint!(factors, length(replacements))
    for rep in replacements
        part_coeff, stripped = separate_coeff_qcomposite(rep)
        coeff *= part_coeff
        if !(stripped isa QAtomProduct && isnumeric(stripped))
            push!(factors, stripped)
        end
    end
    return coeff, factors
end

# contains QComposite Vector (for product), changed, go_left, new_coeff_fun
# ======> Pair Sorting & Simplifications <===========================================================================================================
function simplify_pair_composite(a::S, b::T, qspace::QSpace)::Tuple{Vector{CompositeBranch}, Bool, Bool} where {S <: QComposite, T <: QComposite}
    if b < a && commutes(a,b)
        return CompositeBranch[_make_branch(QComposite[b, a], qspace)], true, true
    else
        return CompositeBranch[], false, true
    end
end
function simplify_pair_composite(a::QAtomProduct, b::QAtomProduct, qspace::QSpace)::Tuple{Vector{CompositeBranch}, Bool, Bool}
    branches = CompositeBranch[]
    for prod in multiply_QAtomProducts(a, b)
        push!(branches, _make_branch(QComposite[prod], qspace))
    end
    return branches, true, true
end 
function simplify_pair_composite(a::QExp, b::QExp, qspace::QSpace)::Tuple{Vector{CompositeBranch}, Bool, Bool} 
    if commutes(a,b)
        new_expr = a.expr+b.expr 
        if isnumeric(new_expr)
            if length(new_expr) == 0 
                coeff_fun = qspace.c_one 
            elseif length(new_expr) == 1
                coeff_fun = exp(new_expr.terms[1].coeff_fun)
            else
                error("Numeric terms should be of length 0 or 1. ")
            end
            return CompositeBranch[(coeff_fun, QComposite[])], true, false
        end
        return CompositeBranch[_make_branch(modify_expr(a, new_expr), qspace)], true, false
    else
        return CompositeBranch[], false, false
    end
end
function simplify_pair_composite(a::QSum, b::QSum, qspace::QSpace)::Tuple{Vector{CompositeBranch}, Bool, Bool}
    return CompositeBranch[_make_branch(a*b, qspace)], true, false
end
function simplify_pair_composite(a::QSum, b::QAtomProduct, qspace::QSpace)::Tuple{Vector{CompositeBranch}, Bool, Bool}
    return CompositeBranch[_make_branch(modify_expr(a, a.expr * b), qspace)], true, false
end
function simplify_pair_composite(b::QAtomProduct, a::QSum, qspace::QSpace)::Tuple{Vector{CompositeBranch}, Bool, Bool}
    return CompositeBranch[_make_branch(modify_expr(a, b * a.expr), qspace)], true, false
end


"""
    add_QComposite_to_QCompositeProduct(terms::Vector{QComposite}, a::QComposite, ss::QSpace) → Vector{Tuple{CFunction, Vector{QComposite}}}  # sum of products

Append `a` and bubble it left:
- if a pair changes (swap / multiply / unify), splice the replacement and step left
- if no change, step left
- branching is preserved
"""
function add_QComposite_to_QCompositeProduct(terms::AbstractVector{<:QComposite}, a::T, ss::QSpace)::Vector{CompositeBranch} where {T<:QComposite}
    seed = QComposite[]
    append!(seed, terms)
    push!(seed, a)

    start_i = length(seed) > 1 ? length(seed) - 1 : 0
    states = Vector{Tuple{CFunction,Vector{QComposite},Int}}()
    push!(states, (ss.c_one, seed, start_i))
    results = CompositeBranch[]

    while !isempty(states)
        next_states = Tuple{CFunction,Vector{QComposite},Int}[]
        for (c, t, i) in states
            if i > 0
                branches, changed, go_left = simplify_pair_composite(t[i], t[i+1], ss)  
                if changed
                    for (dc, pair_terms) in branches
                        nt = Vector{QComposite}()
                        append!(nt, t[1:i-1])
                        append!(nt, pair_terms)
                        append!(nt, t[i+2:end])

                        new_coeff = c * dc
                        if length(nt) <= 1
                            push!(results, (new_coeff, nt))
                        else
                            new_index = go_left ? max(1, i-1) : min(i, length(nt)-1)
                            if new_index >= length(nt)
                                push!(results, (new_coeff, nt))
                            else
                                push!(next_states, (new_coeff, nt, new_index))
                            end
                        end
                    end
                else
                    new_i = i + 1
                    if new_i >= length(t)
                        push!(results, (c, t))
                    else
                        push!(next_states, (c, t, new_i))
                    end
                end
            else
                push!(results, (c, t))
            end
        end
        states = next_states
    end

    return results
end

function multiply_QCompositeProduct_terms(p1::AbstractVector{<:QComposite}, p2::AbstractVector{<:QComposite})
    isempty(p1) && isempty(p2) && error("Cannot multiply two empty composite products.")
    qspace = isempty(p1) ? p2[1].qspace : p1[1].qspace

    seed = QComposite[]
    append!(seed, p1)
    states = CompositeBranch[]
    push!(states, (qspace.c_one, seed))

    for a in p2 
        new_states = CompositeBranch[]
        for (c, t) in states
            for branch in add_QComposite_to_QCompositeProduct(t, a, qspace)
                dc, nt = branch
                push!(new_states, (c * dc, nt))
            end
        end
        states = new_states
    end
    return states
end

function multiply_QCompositeProducts(coeff::CFunction, p1::AbstractVector{<:QComposite}, p2::AbstractVector{<:QComposite})
    states = multiply_QCompositeProduct_terms(p1, p2)
    qspace = isempty(p1) ? p2[1].qspace : p1[1].qspace
    results = QComposite[]
    for (c, terms) in states
        append!(results, _QCompositeProduct(qspace, coeff * c, terms))
    end
    return results
end
function multiply_QCompositeProducts(coeff::CFunction, p1::AbstractVector{<:QComposite}, p2::AbstractVector{<:QComposite}, ::Val{:nosimp})
    states = multiply_QCompositeProduct_terms(p1, p2)
    qspace = isempty(p1) ? p2[1].qspace : p1[1].qspace
    results = QComposite[]
    for (c, terms) in states
        append!(results, _QCompositeProduct(qspace, coeff * c, terms, Val(:nosimp)))
    end
    return results 
end
