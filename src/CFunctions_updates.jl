import .CFunctions: define_cintegral, CFunction
import .QExpressions: ConstrainedIndexBlock
using .QSpaces: QSpace, SubSpaceIndex

function _cintegral_extract_indexes_from_blocks(qspace::QSpace, blocks::Vector{ConstrainedIndexBlock})::Tuple{Vector{SubSpaceIndex}, Vector{Vector{Int}}}
    indexes = SubSpaceIndex[]
    parameter_group_indexes = Vector{Vector{Int}}()
    subspaces = qspace.subspaces

    for block in blocks
        for block_index in block.indexes
            push!(indexes, block_index)
            subspace = subspaces[block_index.outer]
            ensemble = subspace.ensemble
            @assert !isnothing(ensemble) "Summation and integration indexes must be ensemble indexes"
            push!(parameter_group_indexes, copy(ensemble.distribution_group_indices))
        end
    end

    return indexes, parameter_group_indexes
end

function define_cintegral(qspace::QSpace, expr::CFunction, blocks::Vector{ConstrainedIndexBlock})
    indexes, parameter_group_indexes = _cintegral_extract_indexes_from_blocks(qspace, blocks)
    return define_cintegral(qspace.param_info, expr, indexes, parameter_group_indexes)
end
