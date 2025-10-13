using Base: WeakRef
import .CFunctions
import .CFunctions: update_t!
import ..ParameterGroups
using ..ParameterGroups: ParameterGroup, ParameterGroupLike

const _GroupStorageAbstract = Union{ComplexF64, Array{ComplexF64}, Vector{Float64}}

struct AbstractIndexParameters
    qspace::WeakRef
    param_info::CFunctions.ParameterInfo
    group_values::Vector{_GroupStorageAbstract}
    where_which::WhereWhichParamGroup
end
