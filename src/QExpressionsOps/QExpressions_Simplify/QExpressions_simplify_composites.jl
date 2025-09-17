# Helper Functions for QComposite Constructions ==> To keep all coefficients in front  
function separate_coeff_qcomposite(q::T)::Tuple{CFunction, QComposite} where {T <: QComposite}
    return q.coeff_fun, modify_coeff(q, q.statespace.c_one)
end
function separate_coeff_qcomposite(q::QSum)::Tuple{CFunction, QComposite} 
    return q.statespace.c_one, q 
end
function separate_coeff_qcomposites(qs::Vector{QComposite}, statespace::StateSpace)::Tuple{CFunction, Vector{QComposite}} 
    c_one = statespace.c_one 
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

# contains QComposite Vector (for product), changed, go_left, new_coeff_fun
# ======> Pair Sorting & Simplifications <===========================================================================================================
function simplify_pair_composite(a::S, b::T, statespace::StateSpace)::Tuple{Vector{QComposite}, Bool, Bool, Bool} where {S <: QComposite, T <: QComposite}
    if b < a && commutes(a,b)
        return QComposite[b, a], true, true, false
    else
        return QComposite[a, b], false, true, false
    end
end
function simplify_pair_composite(a::QAtomProduct, b::QAtomProduct, statespace::StateSpace)::Tuple{Vector{QComposite}, Bool, Bool, Bool}
    return multiply_QAtomProducts(a, b),  true, true, false
end 
function simplify_pair_composite(a::QExp, b::QExp, statespace::StateSpace)::Tuple{Vector{QComposite}, Bool, Bool, Bool} 
    if commutes(a,b)
        new_expr = a.expr+b.expr 
        if is_numeric(new_expr)
            if length(new_expr) == 0 
                coeff_fun = statespace.c_one 
            elseif length(new_expr) == 1
                coeff_fun = exp(new_expr.terms[1].coeff_fun)
            else
                error("Numeric terms should be of length 0 or 1. ")
            end
            return QComposite[IdentityQAtomProduct(statespace, coeff_fun)], true, false, true
        end
        return QComposite[modify_expr(a, new_expr)], true, false, false
    else
        return QComposite[a, b], false, false, false
    end
end
function simplify_pair_composite(a::QSum, b::QSum, statespace::StateSpace)::Tuple{Vector{QComposite}, Bool, Bool, Bool}
    @inline function hasdupes_sorted(v::AbstractVector)
        for i in 2:length(v)
            if v[i] == v[i-1]
                return true
            end
        end
        return false
    end
    new_indexes = sort!(vcat(a.indexes, b.indexes))
    if hasdupes_sorted(new_indexes)
        error("Two sums cannot have the same summation indexes and be multiplied with one another!")
    end
    return [modify_expr_indexes(a, a.expr*b.expr, new_indexes, Val(:nosort))], true, false, false
end

function simplify_pair_composite(a::QSum, b::QAtomProduct, statespace::StateSpace)::Tuple{Vector{QComposite}, Bool, Bool, Bool}
    return [modify_expr(a, a.expr * b)], true, false, false
end
function simplify_pair_composite(b::QAtomProduct, a::QSum, statespace::StateSpace)::Tuple{Vector{QComposite}, Bool, Bool, Bool}
    return [modify_expr(a, b * a.expr)], true, false, false
end


"""
    add_QComposite_to_QCompositeProduct(terms::Vector{QComposite}, a::QComposite, ss::StateSpace) → Vector{Tuple{ComplexRational, Vector{QComposite}}}  # sum of products

Append `a` and bubble it left:
- if a pair changes (swap / multiply / unify), splice the replacement and step left
- if no change, step left
- branching is preserved
"""
function add_QComposite_to_QCompositeProduct(terms::AbstractVector{<:QComposite}, a::T, ss::StateSpace)::Tuple{CFunction, Vector{QComposite}} where {T<:QComposite}
    seed = QComposite[terms...]
    push!(seed, a)

    start_i = length(seed) > 1 ? length(seed) - 1 : 0
    c = ss.c_one
    t = seed
    i = start_i

    done = false
    while !done 
        if i > 0
            pair, changed, go_left, new_coeff_fun = simplify_pair_composite(t[i], t[i+1], ss)  
            if changed
                nt = Vector{QComposite}()
                append!(nt, t[1:i-1])
                if !new_coeff_fun 
                    append!(nt, pair)          # pair may be length 0, 1, or 2
                else
                    curr_coeff = ss.c_one
                    for p in pair 
                        if !is_numeric(p) 
                            push!(nt, p) 
                        else 
                            curr_coeff *= p.coeff_fun
                        end
                    end
                end
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
                    done = true
                end
                if new_coeff_fun
                    c *= curr_coeff
                end
                t = nt
                i = new_index
            else
                # no local change → move right
                new_i = i+1
                if new_i >= length(t)   
                    done = true
                end
                i = new_i
            end
        else
            done = true
        end
    end
    # Outer vector = sum, inner Vector{QAtom} = product
    return c, t
end

function multiply_QCompositeProduct_terms(p1::AbstractVector{<:QComposite}, p2::AbstractVector{<:QComposite})
    statespace = p1[1].statespace
    c = ComplexRational(1,0,1)
    t = p1
    for a in p2 
        dc, t =  add_QComposite_to_QCompositeProduct(t, a, statespace) 
        c *= dc
    end
    return c, t
end

function multiply_QCompositeProducts(coeff::CFunction, p1::AbstractVector{<:QComposite}, p2::AbstractVector{<:QComposite})
    statespace = p1[1].statespace
    c, t = multiply_QCompositeProduct_terms(p1, p2)
    return [ QCompositeProduct(c * coeff , t)]
end
function multiply_QCompositeProducts(coeff::CFunction, p1::AbstractVector{<:QComposite}, p2::AbstractVector{<:QComposite}, ::Val{:nosimp})
    statespace = p1[1].statespace
    c, t = multiply_QCompositeProduct_terms(p1, p2)
    return [ QCompositeProduct(c * coeff , t, Val(:nosimp)) ]
end