module QEqSets

using ..QExpressions: diffQEq, diffQEqOrdered, QExpr, QAtomProduct, QAtomOrdered, OrderedDiffQEq
import ..QExpressions: Order

export diffQEqSet, diffQEqSetOrdered, diff_QEqSet, OrderedDiffQEqSet

struct diffQEqSet
    equations::Vector{diffQEq}
    loss::Function
    function diffQEqSet(equations::Vector{diffQEq}, loss::Function, ::Val{:raw})
        new(copy(equations), loss)
    end
end

struct diffQEqSetOrdered
    equations::Vector{diffQEqOrdered}
    loss::Function
    function diffQEqSetOrdered(equations::Vector{diffQEqOrdered}, loss::Function, ::Val{:raw})
        new(copy(equations), loss)
    end
end

end
