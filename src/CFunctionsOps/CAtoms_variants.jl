import ..ConcreteIndexes
import Base: LinearIndices, pointer, pointer_from_objref, fieldoffset, unsafe_load
import Base.GC

function _collect_index_values(indices::ConcreteIndexes, ens_indices::AbstractVector)
    values = Int[]
    for ens_idx in ens_indices
        ensemble = ens_idx.outer
        inner = ens_idx.inner
        ensemble <= length(indices.indices) ||
            error("Concrete indices missing ensemble $(ensemble).")
        ensemble_entries = indices.indices[ensemble]
        inner <= length(ensemble_entries) ||
            error("Concrete indices missing entry $(inner) in ensemble $(ensemble).")
        idx = ensemble_entries[inner]
        idx > 0 || error("Concrete index for ensemble $(ensemble) inner $(inner) not set.")
        push!(values, idx)
    end
    return values
end

"""
    CAtomIndexed(param_info, coeff, var_exponents, indices)

Coefficient atom that keeps concrete ensemble indices next to its sparse
exponent vector. `indices` may be a [`ConcreteIndexes`](@ref) instance or any
vector of integer vectors aligned with the ensemble layout of `param_info`.
"""
struct CAtomIndexed <: CAtomic
    param_info::ParameterInfo
    coeff::ComplexRational
    var_exponents::SparseVector{Int,Int}
    indices::ConcreteIndexes
end

function CAtomIndexed(atom::CAtom, indices::ConcreteIndexes)
    return CAtomIndexed(atom.param_info, atom.coeff, atom.var_exponents, indices)
end

function CAtomIndexed(atom::CAtom, indices::AbstractVector{<:AbstractVector{<:Integer}})
    concrete = ConcreteIndexes(atom.param_info, indices)
    return CAtomIndexed(atom, concrete)
end
var_exponents(a::CAtomIndexed) = a.var_exponents
coeff(a::CAtomIndexed)::Vector{ComplexRational} = [a.coeff]
length(::CAtomIndexed) = 1


# --------------------------------------------------------------------------------------------------
# CAtomReferenced — reference variant
# --------------------------------------------------------------------------------------------------

@inline function _factor_value(ptr::Ptr{Float64}, anchor)
    value = 0.0
    GC.@preserve anchor begin
        value = unsafe_load(ptr)
    end
    return value
end

# Creates the pointer
@inline function _factor_anchor_pointer(gv::GroupVals{Float64, OT}, sample_indices::Vector{Int}) where {OT}
    isempty(sample_indices) || error("Scalar parameter cannot carry ensemble indices.")
    anchor = gv
    ptr = Ptr{Float64}(0)
    GC.@preserve anchor begin
        base = Ptr{UInt8}(pointer_from_objref(anchor))
        ptr = Base.unsafe_convert(Ptr{Float64}, base + fieldoffset(typeof(anchor), 1))
    end
    return anchor, ptr
end
@inline function _factor_anchor_pointer(storage::Vector{Float64}, sample_indices::Vector{Int})
    length(sample_indices) == 1 ||
        error("Expected one sample index for vector parameter, got $(length(sample_indices)).")
    idx = sample_indices[1]
    anchor = storage
    ptr = Ptr{Float64}(0)
    GC.@preserve anchor begin
        ptr = pointer(anchor, idx)
    end
    return anchor, ptr
end
@inline function _factor_anchor_pointer(storage::Array{Float64, N}, sample_indices::Vector{Int}) where {N}
    length(sample_indices) == N ||
        error("Expected $(N) sample indices for array parameter, got $(length(sample_indices)).")
    idx_tuple = ntuple(i -> sample_indices[i], N)
    lin_idx = LinearIndices(storage)[idx_tuple...]
    anchor = storage
    ptr = Ptr{Float64}(0)
    GC.@preserve anchor begin
        ptr = pointer(anchor, lin_idx)
    end
    return anchor, ptr
end

struct CAtomReferenced{Mode<:ParameterValuesMode} <: CAtomic
    param_info::ParameterInfo
    coeff::ComplexRational
    ptrs::Vector{Ptr{Float64}}
    anchors::Vector{Any}
    param_indices::Vector{Int}
    exponents::Vector{Int}
    pv::ParameterValues{Mode}
end

function _push_factor!(::Type{Mode}, pv::ParameterValues{Mode}, atom::CAtomIndexed,
                       param_index::Int, exponent::Int,
                       ptrs::Vector{Ptr{Float64}}, anchors::Vector{Any},
                       param_indices::Vector{Int}, exponents::Vector{Int}) where {Mode<:ParameterValuesMode}
    param = atom.param_info.params[param_index]
    group_index = param.group_index
    group = pv.groups[group_index]
    group.of_t && error("CAtomReferenced does not support time-dependent parameter groups.")
    sample_indices = _collect_index_values(atom.indices, param.ensemble_indices)
    gv = pv.group_values[group_index]
    storage = gv.value
    anchor, ptr = _factor_anchor_pointer(gv, sample_indices)
    push!(ptrs, ptr)
    push!(anchors, anchor)
    push!(param_indices, param_index)
    push!(exponents, exponent)
    return nothing
end

function _atom_referenced_components(pv::ParameterValues{Mode}, atom::CAtomIndexed) where {Mode<:ParameterValuesMode}
    ptrs = Ptr{Float64}[]
    anchors = Any[]
    param_indices = Int[]
    exponents = Int[]
    exps = atom.var_exponents
    for idx in exps.nzind
        exponent = exps[idx]
        _push_factor!(Mode, pv, atom, idx, exponent, ptrs, anchors, param_indices, exponents)
    end
    return ptrs, anchors, param_indices, exponents
end

function CAtomReferenced(pv::ParameterValues{Mode}, atom::CAtomIndexed) where {Mode<:ParameterValuesMode}
    ptrs, anchors, param_indices, exponents = _atom_referenced_components(pv, atom)
    return CAtomReferenced{Mode}(atom.param_info, atom.coeff, ptrs, anchors, param_indices, exponents, pv)
end

function CAtomReferenced(pv::ParameterValues{Mode}, atom::CAtom, indices::ConcreteIndexes) where {Mode<:ParameterValuesMode}
    return CAtomReferenced(pv, CAtomIndexed(atom, indices))
end

function CAtomReferenced(pv::ParameterValues{Mode}, atom::CAtom, indices::AbstractVector{<:AbstractVector{<:Integer}}) where {Mode<:ParameterValuesMode}
    return CAtomReferenced(pv, CAtomIndexed(atom, indices))
end

var_exponents(a::CAtomReferenced) = begin
    exps = spzeros(Int, a.param_info.dims)
    for (param_idx, exponent) in zip(a.param_indices, a.exponents)
        exps[param_idx] = exponent
    end
    return exps
end

coeff(a::CAtomReferenced) = [a.coeff]
length(::CAtomReferenced) = 1

# --------------------------------------------------------------------------------------------------
# CEval — fully evaluated constant
# --------------------------------------------------------------------------------------------------

"""
    CEval

Constant evaluation leaf that stores a fully realised complex value for a time-independent
coefficient expression.
"""
struct CEval <: CAtomic
    param_info::ParameterInfo
    value::ComplexF64
end

"""
    CEval(param_info::ParameterInfo, value)

Construct a constant evaluation node from `value`, storing it as `ComplexF64`.
"""
CEval(param_info::ParameterInfo, value) = CEval(param_info, ComplexF64(value))

"""
    CEval(atom::CAtomReferenced)

Collapse a time-independent referenced atom into a constant `CEval` using the cached
parameter slots gathered by `CAtomReferenced`.
"""
function CEval(atom::CAtomReferenced{Mode}) where {Mode<:ParameterValuesMode}
    prod_val = 1.0
    for (ptr, anchor, exponent) in zip(atom.ptrs, atom.anchors, atom.exponents)
        base = _factor_value(ptr, anchor)
        term = exponent == 1 ? base : base ^ exponent
        prod_val *= term
    end
    value = ctimes(atom.coeff, prod_val)
    return CEval(atom.param_info, ComplexF64(value))
end

var_exponents(c::CEval) = spzeros(Int, c.param_info.dims)
coeff(c::CEval) = ComplexF64[c.value]
length(::CEval) = 1
