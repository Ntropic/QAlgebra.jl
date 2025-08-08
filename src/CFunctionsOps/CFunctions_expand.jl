export expand
# function expand(r::CRational, order


"""
    expand(f::CFunction, mode::Symbol, target::Type{T}, args...) where T<:CFunction

Walk the symbolic AST `f` and apply a mode-specific expansion whenever a node of type `target` is encountered.

# Arguments
- `f::CFunction`  
    The root of your symbolic expression tree.
- `mode::Symbol`  
    The expansion mode to apply. Supported modes:
    - `:Taylor`     — Perform a Taylor-series expansion on each `CExp` or `CLog` node.  
                      Requires an integer `order` (e.g. `expand(f, :Taylor, CExp, 3)`).
    - `:Rational`   — Distribute `CRational` nodes over any `CSum` in their numerator  
                      (i.e. `(A+B)/D → A/D + B/D`).
    - `:Log`        — Apply algebraic log rules to `CLog` nodes, flattening  
                      `log(a*b/c) → log(a)+log(b)-log(c)` for any mix of `CProd`/`CRational`.
                      The type CLog can be ommited for this mode.
- `target::Type{T}`  
    The specific subtype of `CFunction` to expand (e.g. `CExp`, `CRational`, or `CLog`).
- `args...`  
    Additional parameters required by the mode (for example, the Taylor series `order::Int`).

# Returns
A new `CFunction` tree with the requested expansions applied only at nodes of type `target`.

# Examples
```julia
# 3rd-order Taylor expansion of all exponentials:
g = expand(expr, :Taylor, CExp, 3)

# Distribute rational sums:
h = expand(expr, :Rational, CRational)

# Distribute logs over products and quotients:
l = expand(expr, :Log, CLog)
"""
function expand(f::CFunction, ::Val{:Log})
    newf, _ = _expand(f, Val(:Log), CLog)
    return newf
end
function expand(f::CFunction, mode::Symbol, target::Type{T}, args...) where T<:CFunction
    newf, _ = _expand(f, Val(mode), target, args...)
    return newf
end

_expand(s::CAtom,    ::Val{M}, ::Type{T}, args...) where {M,T} = (s, false)
_expand(f::CFunction, ::Val{M}, ::Type{T}, args...) where {M,T} = (f, false)

function _expand(s::CSum, ::Val{M}, ::Type{T}, args...) where {M,T}
    terms2::Vector{CFunction} = []
    any_exp = false
    for (i, term) in enumerate(s.terms)
        t2, e = _expand(term, Val(M), T, args...)
        if isa(t2, CSum)
            append!(terms2, t2.terms)
        else
            push!(terms2, t2)
        end
        any_exp |= e
    end
    # if nothing changed, just return original
    return any_exp ? (CSum(terms2), true) : (s, false)
end

function _expand(p::CProd, ::Val{M}, ::Type{T}, args...) where {M,T}
    terms2 = Vector{CFunction}(undef, length(p.terms))
    any_exp = false
    for (i, term) in enumerate(p.terms)
        t2, e = _expand(term, Val(M), T, args...)
        terms2[i] = t2
        any_exp |= e
    end
    return any_exp ? (CProd(p.coeff, terms2), true) : (p, false)
end

function _expand(r::CRational, ::Val{M}, ::Type{T}, args...) where {M,T}
    n2, e1 = _expand(r.numer,   Val(M), T, args...)
    d2, e2 = _expand(r.denom,   Val(M), T, args...)
    any_exp = e1 || e2
    return any_exp ? (CRational(n2, d2), true) : (r, false)
end

function _expand(e::CExp, ::Val{M}, ::Type{T}, args...) where {M,T}
    new_x, changed = _expand(e.x, Val(M), T, args...)
    if !changed
        return (e, false)
    end
    return (CExp(e.coeff, new_x), true)
end
function _expand(e::CLog, ::Val{M}, ::Type{T}, args...) where {M,T}
    new_x, changed = _expand(e.x, Val(M), T, args...)
    if !changed
        return (e, false)
    end
    return (CLog(e.coeff, new_x), true)
end

# <========> Specialized Functions <==================================================================>

function _expand(e::CExp, ::Val{:Taylor}, ::Type{CExp}, order::Int)   # coeff * ∑_{n=0}^order x^n / n!
    new_x, changed = _expand(e.x, Val(:Taylor), CExp, order)
    if order <= 0 
        error("Order of Taylor expansion must be positive!")
    end 
    #order 0 
    terms::Vector{CFunction} = [CAtom(ComplexRational(1,0,1), zeros(Int, dims(e)))]
    # order 1 
    coeff = e.coeff
    if isa(new_x, CSum)
        append!(terms, (coeff * new_x).terms)
    else
        push!(terms, coeff * new_x)
    end 
    curr_x = copy(new_x)
    @inbounds for n in 2:order # higher orders
        coeff /= n 
        new_x = new_x * curr_x  # creates the exponent
        if isa(new_x, CSum)
            append!(terms, (new_x * coeff).terms)
        else
            push!(terms, new_x * coeff)
        end
    end
    # create a sum of the terms 
    return (CExp(e.coeff, CSum(terms)), true)
end

function _expand(r::CRational, ::Val{:Rational}, ::Type{CRational}, args...) #Distribute a fraction over a sum in the numerator: (A + B + …) / D  ⇒  A/D  +  B/D  +  …
    new_numer, changed = _expand(r.numer, Val(:Rational), CRational, args...)
    new_denom, changed2 = _expand(r.denom, Val(:Rational), CRational, args...)
    if new_numer isa CSum
        # pre-allocate
        out::Vector{CFunction} = []
        @inbounds for i in 1:length(new_numer.terms)
            new = r.numer.terms[i] / r.denom
            if isa(new, CSum)
                append!(out, new.terms)
            else
                push!(out, new)
            end
        end
        return (CSum(out), true)
    else
        # nothing to expand
        either_changed = changed || changed2
        if either_changed
            return (CRational(new_numer, new_denom), either_changed)
        else
            return (r, false)
        end
    end
end

function _expand(e::CLog, ::Val{:Taylor}, ::Type{CLog}, order::Int)  # log(1+Δ) = ∑_{n=1}^order  (-1)^(n+1) Δ^n / n
    # first descend into the argument
    new_x, _ = _expand(e.x, Val(:Taylor), CLog, order)
    # validate
    if order <= 0
        error("Order of Taylor expansion must be positive!")
    end

    # build Δ = new_x - 1
    one_atom = CAtom(ComplexRational(1,0,1), zeros(Int, dims(new_x)))
    delta = new_x - one_atom
    curr_delta = copy(delta)
    # accumulate terms
    terms::Vector{CFunction} = CFunction[]

    if curr_delta isa CSum
        append!(terms, (curr_delta * e.coeff).terms)
    else
        push!(terms, curr_delta * e.coeff)
    end
    @inbounds for n in 2:order
        # sign = (-1)^(n+1)
        sign = isodd(n) ?  1 : -1
        coef = e.coeff * ComplexRational(sign, 0, n)
        curr_delta = curr_delta * delta
        if curr_delta isa CSum
            append!(terms, (curr_delta * coef).terms)
        else
            push!(terms, curr_delta * coef)
        end
    end

    return (CSum(terms), true)
end



function log_distribute(f::CFunction)
    return [copy(f)] 
end
function log_distribute(s::CSum) 
    return copy.(s.terms) 
end
function log_distribute(f::CRational)
    terms = log_distribute(f.numer)
    append!(terms, .-log_distribute(f.denom)) # negate the log of the denominator
end

function log_distribute(f::CProd) 
    terms = log_distribute.(f.terms) 
    terms[1] *= f.coeff
end
function _expand(l::CLog, ::Val{:Log}, ::Type{CLog}, args...)
    # first recurse so nested CLogs get handled too
    x2, _ = _expand(l.x, Val(:Log), CLog, args...)

    # flatten everything into factors
    facs = log_distribute(x2)

    # wrap each in a CLog with ±coeff
    logs = [ CLog(l.coeff, f) for f in facs ]

    return (CSum(logs), true)
end