module QEqSets

using ..QExpressions: diffQEq, Expectation, QObj, QComposite, QAtomProduct, QExpr

export diffQEqSet



struct diffQEqSet
    equations::Vector{diffQEq}
    loss::Function
    function diffQEqSet(equations::Vector{diffQEq}, loss::Function, ::Val{:raw})
        new(equations, loss)
    end
end

function diffQEqSet(; kwargs...)
    @warn "diffQEqSet construction is not implemented yet; returning an empty set."
    return diffQEqSet(diffQEq[], () -> nothing, Val(:raw))
end

function diffQEqSet(equations::Vector{diffQEq}; loss::Function = () -> nothing)
    return diffQEqSet(equations, loss, Val(:raw))
end

function diffQEqSet(equations::Vector{diffQEq}, loss::Function)
    return diffQEqSet(equations, loss, Val(:raw))
end


end
