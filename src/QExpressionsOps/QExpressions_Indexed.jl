import ..ConcreteIndexes
import ..CFunctions: CFunction, CAtomIndexed, to_stringer
import ..StringUtils: indices2str
import ..CFunctions: Indexed

"""
    QAtomIndexed(qspace, coeff_fun, op_indices, ensemble_indices, indices, time_index)

Ordered atom where ensemble operators carry concrete index assignments. The
`indices` argument may be a [`ConcreteIndexes`](@ref) instance or any
vector-of-vectors matching the ensemble layout of `qspace`.
"""
struct QAtomIndexed <: QComposite
    qspace::QSpace
    coeff_fun::CFunction
    op_indices::Vector{Vector{Is}}
    ensemble_indices::Vector{Vector{Int}}
    concrete_indices::ConcreteIndexes
    time_index::Int
    function QAtomIndexed(qspace::QSpace, coeff_fun::CFunction, op_indices::Vector{Vector{Is}},
                          ensemble_indices::Vector{Vector{Int}}, concrete_indices::ConcreteIndexes, time_index::Int)
        param_info = qspace.param_info
        concrete_indices.expected_lengths == param_info.how_many_by_ensemble ||
            error("Concrete indices do not match the ensemble sizes of the provided QSpace.")
        return new(qspace, coeff_fun, op_indices, ensemble_indices, concrete_indices, time_index)
    end
end

function QAtomIndexed(atom::QAtomOrdered, indices)
    param_info = atom.qspace.param_info
    concrete = indices isa ConcreteIndexes ? indices : ConcreteIndexes(param_info, indices)
    return QAtomIndexed(atom.qspace, atom.coeff_fun, atom.op_indices, atom.ensemble_indices, concrete, atom.time_index)
end

OrderbyOperator(q::QAtomIndexed) = q
