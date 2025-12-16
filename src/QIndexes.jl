module QIndexes

export AbstractIndex, TimeIndex, SampleIndex, index_label, index_tag

abstract type QIndex end

"""
    AbstractIndex(; label, subspace=nothing, slot=nothing, dag=false)

Symbolic index used to tag operator factors that live in a specific subspace.

- `subspace`: outer-subspace identifier (1-based).
- `ensemble`: ensemble identifier (1-based) 
- `slot`: optional inner position within the subspace (1-based).
"""
struct AbstractIndex <: QIndex
    subspace::Int
    ensemble::Int  # defined by subspace, and hence doesn'T need to be checked in equality and isless checks
    slot::Int
    summation::Bool
    function AbstractIndex(subspace::Int, ensemble::Int, slot::Int, summation::Bool=false)
        return new(subspace, ensemble, slot, summation)
    end
end
function Base.isequal(a::AbstractIndex, b::AbstractIndex)::Bool 
    return a.subspace == b.subspace && a.slot == b.slot && a.summation == b.summation
end
function Base.isless(a::AbstractIndex, b::AbstractIndex)::Bool 
    if a.subspace != b.subspace
        return a.subspace < b.subspace
    end
    if a.slot != b.slot
        return a.slot < b.slot
    end
    return !a.summation && b.summation
end

"""
    TimeIndex(; label=:t, order=0)

Specialised index that tracks explicit time labels.

- `order`: derivative order or discrete time slot. (0 based -> -1for no time index)
"""
struct TimeIndex 
    order::Int 
end
function Base.isequal(a::TimeIndex, b::TimeIndex)::Bool 
    return a.order == b.order 
end
function Base.isless(a::TimeIndex, b::TimeIndex)::Bool 
    return a.order < b.order 
end

"""
    SampleIndex(; label, ensemble, slot)

Index identifying a sample drawn from an ensemble subspace.

- `subspace`: outer-subspace identifier (1-based).
- `ensemble`: ensemble identifier (1-based).
- `slot`: position inside the ensemble (1-based).
"""
struct SampleIndex <: QIndex
    subspace::Int
    ensemble::Int  # defined by subspace, and hence doesn'T need to be checked in equality and isless checks
    slot::Int
    summation::Bool
    function SampleIndex(subspace::Int, ensemble::Int, slot::Int, summation::Bool=false)
        return new(subspace, ensemble, slot, summation)
    end
end
function Base.isequal(a::SampleIndex, b::SampleIndex)::Bool 
    return a.subspace == b.subspace && a.slot == b.slot && a.summation == b.summation
end
function Base.isless(a::SampleIndex, b::SampleIndex)::Bool 
    if a.subspace != b.subspace 
        return a.subspace < b.subspace 
    end
    if a.slot != b.slot
        return a.slot < b.slot
    end
    return !a.summation && b.summation
end

end # module
