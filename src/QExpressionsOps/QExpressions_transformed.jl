# Types of transformed expressions helpers for explicit indexing. 

struct QAtomOrdered <: QComposite
    qspace::QSpace
    coeff_fun::CFunction
    op_indices::Vector{QAtom}
    permutation::Vector{Int}
    function QAtomOrdered(qspace::QSpace, coeff_fun::CFunction, op_indices::Vector{QAtom}, permutation::Vector{Int})
        length(op_indices) == length(permutation) || error("Permutation length does not match operator count.")
        return new(qspace, coeff_fun, copy(op_indices), copy(permutation))
    end
end
permutation(q::QAtomOrdered) = q.permutation

"""
    OrderedQAtomProduct(q::QAtomProduct; lt=isless)

Return a `QAtomOrdered` where the operators of `q` are sorted according to `lt`.
The permutation field records how the sorted ordering maps back to the original
operator sequence.
"""
function OrderedQAtomProduct(q::QAtomProduct; lt=isless)
    n = length(q.expr)
    if n == 0
        return QAtomOrdered(q.qspace, q.coeff_fun, QAtom[], Int[])
    end
    perm = sortperm(q.expr; lt=lt)
    ordered_atoms = Vector{QAtom}(undef, n)
    @inbounds for i in 1:n
        ordered_atoms[i] = q.expr[perm[i]]
    end
    return QAtomOrdered(q.qspace, q.coeff_fun, ordered_atoms, perm)
end