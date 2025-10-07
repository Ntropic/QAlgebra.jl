import ..ConcreteIndexes
import ..CFunctions: CFunction, CAtomIndexed, to_stringer
import ..StringUtils: indexes2str
import ..CFunctions: Indexed

"""
    QAtomIndexed(qspace, coeff_fun, op_indices, ensemble_indexes, indexes, time_index)

Ordered atom where ensemble operators carry concrete index assignments. The
`indexes` argument may be a [`ConcreteIndexes`](@ref) instance or any
vector-of-vectors matching the ensemble layout of `qspace`.
"""
struct QAtomIndexed <: QComposite
    qspace::QSpace
    coeff_fun::CFunction
    op_indices::Vector{Vector{Is}}
    ensemble_indexes::Vector{Vector{Int}}
    concrete_indexes::ConcreteIndexes
    time_index::Int
    function QAtomIndexed(qspace::QSpace, coeff_fun::CFunction, op_indices::Vector{Vector{Is}},
                          ensemble_indexes::Vector{Vector{Int}}, concrete_indexes::ConcreteIndexes, time_index::Int)
        param_info = qspace.param_info
        concrete_indexes.expected_lengths == param_info.how_many_by_ensemble ||
            error("Concrete indexes do not match the ensemble sizes of the provided QSpace.")
        return new(qspace, coeff_fun, op_indices, ensemble_indexes, concrete_indexes, time_index)
    end
end

function QAtomIndexed(atom::QAtomOrdered, indexes)
    param_info = atom.qspace.param_info
    concrete = indexes isa ConcreteIndexes ? indexes : ConcreteIndexes(param_info, indexes)
    return QAtomIndexed(atom.qspace, atom.coeff_fun, atom.op_indices, atom.ensemble_indexes, concrete, atom.time_index)
end

OrderbyOperator(q::QAtomIndexed) = q
