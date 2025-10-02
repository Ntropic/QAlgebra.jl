module OffsetArrays

using Base: OneTo


"""
    OffsetArray{T,N}

Dense, offset-aware storage for children indexed by integer tuples.

- `offset` stores the (hard) minimum valid index per axis.
- `max_ints` only sets the *initial* maximum; writing beyond it expands the array.
- Elements default to `nothing` and storage grows on demand (upper side only).
"""
mutable struct OffsetArray{T,N} 
    data::Vector{Union{Nothing,T}}
    dims::NTuple{N,Int}          # lengths per axis
    offset::NTuple{N,Int}        # hard lower bound per axis
    count::Int                   # number of non-nothing entries
end

# ── Constructors ───────────────────────────────────────────────────────────────

"""
    OffsetArray{T}(min_ints::AbstractVector{Int}, max_ints::AbstractVector{Int})

`max_ints` defines the initial upper bound; it is *not* a cap.
Use `-1` in `max_ints` to start with length 0 along that axis.
"""
function OffsetArray{T}(min_ints::AbstractVector{Int}, max_ints::AbstractVector{Int}) where {T}
    length(min_ints) == length(max_ints) || error("min_ints and max_ints must have the same length")
    N = length(min_ints)

    mins_tuple = ntuple(i -> min_ints[i], N)

    dims = ntuple(i -> begin
        maxv = max_ints[i]
        if maxv == -1
            0
        else
            maxv >= min_ints[i] || error("max_ints[$i] must be ≥ min_ints[$i] (or -1)")
            maxv - min_ints[i] + 1
        end
    end, N)

    total = N == 0 ? 0 : prod(dims)
    data  = Vector{Union{Nothing,T}}(undef, total)
    fill!(data, nothing)
    return OffsetArray{T,N}(data, dims, mins_tuple, 0)
end

# 1-D convenience
OffsetArray{T}(min::Int, max::Int) where {T} = OffsetArray{T}(Int[min], Int[max])

# ── Basic interface ────────────────────────────────────────────────────────────

Base.ndims(::OffsetArray{T,N}) where {T,N} = N

@inline function Base.axes(children::OffsetArray{T,N}) where {T,N}
    N == 0 && return ()
    return ntuple(i -> begin
        len = children.dims[i]
        len == 0 ? (children.offset[i]:(children.offset[i]-1)) :
                   (children.offset[i]:(children.offset[i] + len - 1))
    end, N)
end

Base.size(children::OffsetArray{T,N}) where {T,N} = children.dims
Base.IndexStyle(::Type{<:OffsetArray}) = Base.IndexCartesian()

# ── Internals ──────────────────────────────────────────────────────────────────

@inline function _linear_index(dims::NTuple{N,Int}, subs::NTuple{N,Int}) where {N}
    ranges = ntuple(i -> OneTo(dims[i]), N)
    return LinearIndices(ranges)[subs...]
end

# Resize (upper-side only) and preserve data. Offsets remain fixed.
function _resize_upper!(children::OffsetArray{T,N}, new_dims::NTuple{N,Int}) where {T,N}
    old_dims   = children.dims
    old_total  = N == 0 ? 0 : prod(old_dims)
    new_total  = N == 0 ? 0 : prod(new_dims)

    # no-op if nothing changes
    new_dims == old_dims && return

    new_data = Vector{Union{Nothing,T}}(undef, new_total)
    fill!(new_data, nothing)

    if old_total > 0
        old_ranges = ntuple(i -> OneTo(old_dims[i]), N)
        new_ranges = ntuple(i -> OneTo(new_dims[i]), N)
        li_old = LinearIndices(old_ranges)
        li_new = LinearIndices(new_ranges)

        # copy block as-is (offset didn't change)
        for rel in CartesianIndices(old_ranges)
            old_idx = li_old[Tuple(rel)...]
            v = children.data[old_idx]
            v === nothing && continue
            new_idx = li_new[Tuple(rel)...]
            new_data[new_idx] = v
        end
    end

    children.data = new_data
    children.dims = new_dims
end

# Ensure all given indices are within the current upper bounds (expanding as needed).
function _ensure_upper!(children::OffsetArray{T,N}, idxs::Vararg{Int,N}) where {T,N}
    requested = ntuple(i -> idxs[i], N)
    # enforce hard lower bound
    @inbounds for i in 1:N
        requested[i] < children.offset[i] &&
            error("Index $(requested[i]) below minimum offset $(children.offset[i]) on axis $i")
    end

    new_dims = collect(children.dims)
    need = false
    @inbounds for i in 1:N
        lo  = children.offset[i]
        len = children.dims[i]
        hi  = lo + len - 1
        if len == 0 || requested[i] > hi
            new_len = requested[i] - lo + 1
            new_dims[i] = new_len
            need = true
        end
    end
    need && _resize_upper!(children, Tuple(new_dims))
    return children
end

# ── Indexing ───────────────────────────────────────────────────────────────────

function Base.getindex(children::OffsetArray{T,N}, idxs::Vararg{Int,N}) where {T,N}
    N == 0 && return nothing
    # out-of-bounds → nothing
    @inbounds for i in 1:N
        len = children.dims[i]
        len == 0 && return nothing
        lo  = children.offset[i]
        hi  = lo + len - 1
        (idxs[i] < lo || idxs[i] > hi) && return nothing
    end
    rel = ntuple(i -> idxs[i] - children.offset[i] + 1, N)
    return children.data[_linear_index(children.dims, rel)]
end

function Base.setindex!(children::OffsetArray{T,N}, value, idxs::Vararg{Int,N}) where {T,N}
    _ensure_upper!(children, idxs...)             # auto-expand upward as needed
    rel = ntuple(i -> idxs[i] - children.offset[i] + 1, N)
    li  = _linear_index(children.dims, rel)
    old = children.data[li]
    if old === nothing && value !== nothing
        children.count += 1
    elseif old !== nothing && value === nothing
        children.count -= 1
    end
    children.data[li] = value
    return value
end

Base.isempty(children::OffsetArray) = children.count == 0

# ── Simple tests (comment out in production) ───────────────────────────────────
if abspath(PROGRAM_FILE) == @__FILE__
    # 2D, lower bound fixed at (0,0). Start with initial max (1,1).
    A = OffsetArray{Int}([0,0], [1,1])
    @assert axes(A) == (0:1, 0:1)

    # write inside → no resize
    A[0,0] = 10
    @assert A[0,0] == 10

    # write beyond initial max → expands (upper only)
    A[3,1] = 1
    @assert axes(A) == (0:3, 0:1)
    @assert A[3,1] == 1

    # lower-than-min should throw
    thrown = false
    try
        A[-1,0] = 9
    catch
        thrown = true
    end
    @assert thrown
end

end # module
