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

# vector, no time: one sample index; t is dummy and must be 0
@propagate_inbounds function Base.getindex(G::GroupVals{Vector{Float64}, false}, t::Int, idx)
    @boundscheck checkbounds(G.value, idx)
    @inbounds return G.value[idx]
end

# array, time-first: (t, I...) with time translated by +1
@propagate_inbounds function Base.getindex(G::GroupVals{Array{Float64}, true}, t::Int, I...)
    t1 = t + 1
    @boundscheck checkbounds(G.value, t1, I...)
    @inbounds return G.value[t1, I...]
end

# array, no time: (t, I...) with t==0, forward I...
@propagate_inbounds function Base.getindex(G::GroupVals{Array{Float64}, false}, t::Int, I...)
    @boundscheck checkbounds(G.value, I...)
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
