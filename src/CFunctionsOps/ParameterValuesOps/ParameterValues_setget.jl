# storage types you support
const _GroupStorage = Union{Float64, Vector{Float64}, Array{Float64}}

# Wrapper that carries storage type S and of_t flag OT in the type.
# Mutable so Float64 payloads can be reassigned.
mutable struct GroupVals{S<:_GroupStorage, OT}
    value::S
end

# ergonomic constructor from Bool
GroupVals(value::S, of_t::Bool) where {S<:_GroupStorage} = GroupVals{S, of_t}(value)

const GroupValsAny = GroupVals{S,OT} where {S<:_GroupStorage, OT}

# (optional) basic array interface
Base.eltype(::Type{GroupVals{S,OT}}) where {S,OT} = Float64

function Base.size(G::GroupVals{Float64, OT}) where {OT}
    return (1,)  # expose scalar as length-1 along "time"
end
function Base.size(G::GroupVals{Vector{Float64}, true})
    return (length(G.value),)
end
function Base.size(G::GroupVals{Vector{Float64}, false})
    return (1, length(G.value))
end
function Base.size(G::GroupVals{Array{Float64}, true})
    return size(G.value)
end
function Base.size(G::GroupVals{Array{Float64}, false})
    return (1, size(G.value)...)
end
Base.axes(G::GroupVals) = map(Base.OneTo, size(G))

using Base: @propagate_inbounds, checkbounds

# -------------------- getindex --------------------

# allow explicit empty trailing indices call
@propagate_inbounds function Base.getindex(G::GroupVals{Float64, OT}, t::Int, I...) where {OT}
    @boundscheck @assert isempty(I) "Scalar takes no sample indices"
    return getindex(G, t)
end

@inline _shift_time_index(t::Integer) = t + 1
@inline _shift_time_index(::Colon) = error("Colon slicing not supported for time axis; access individual time indices instead.")
@inline _shift_time_index(t::AbstractRange{<:Integer}) = error("Range slicing not supported for time axis; access individual time indices instead.")
@inline _shift_time_index(t::AbstractVector{<:Integer}) = error("Vector slicing not supported for time axis; access individual time indices instead.")
@inline _shift_time_index(t::Tuple) = map(_shift_time_index, t)

@propagate_inbounds function Base.getindex(G::GroupVals{Vector{Float64}, true}, t, I...)
    isempty(I) || @assert false "Unexpected sample indices for time-dependent vector"
    ti = _shift_time_index(t)
    @inbounds return G.value[ti]
end

# vector, no time: ignore dummy t (allow 0 or colon)
@propagate_inbounds function Base.getindex(G::GroupVals{Vector{Float64}, false}, t, I...)
    t == 0 || t === nothing || error("Time index must be 0 for time-independent vectors")
    @inbounds return G.value[I...]
end

# array, time-first: support slices by shifting time index appropriately
@propagate_inbounds function Base.getindex(G::GroupVals{Array{Float64, N}, true}, t, I...) where {N}
    ti = _shift_time_index(t)
    @inbounds return G.value[ti, I...]
end

# array, no time: ignore dummy t (allow 0 or colon)
@propagate_inbounds function Base.getindex(G::GroupVals{Array{Float64, N}, false}, t, I...) where {N}
    t == 0 || t === nothing || error("Time index must be 0 for time-independent arrays")
    @inbounds return G.value[I...]
end

# -------------------- setindex! --------------------

@propagate_inbounds function Base.setindex!(G::GroupVals{Float64, OT}, v, t::Int, I...) where {OT}
    @boundscheck @assert isempty(I)
    return G.value = v
end

# vector, time-first (translate t + 1)
@propagate_inbounds function Base.setindex!(G::GroupVals{Vector{Float64}, true}, v, t::Int, idx)
    t1 = t + 1
    @boundscheck checkbounds(G.value, t1)
    @inbounds G.value[t1] = v
    return v
end

# vector, no time (t must be 0)
@propagate_inbounds function Base.setindex!(G::GroupVals{Vector{Float64}, false}, v, t::Int, idx)
    @boundscheck checkbounds(G.value, idx)
    @inbounds G.value[idx] = v
    return v
end

# array, time-first (translate t + 1)
@propagate_inbounds function Base.setindex!(G::GroupVals{Array{Float64}, true}, v, t::Int, I...)
    t1 = t + 1
    @boundscheck checkbounds(G.value, t1, I...)
    @inbounds G.value[t1, I...] = v
    return v
end

# array, no time (t must be 0)
@propagate_inbounds function Base.setindex!(G::GroupVals{Array{Float64}, false}, v, t::Int, I...)
    @boundscheck checkbounds(G.value, I...)
    @inbounds G.value[I...] = v
    return v
end
