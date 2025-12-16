struct SubSpaceIndex
    outer::Int
    inner::Int
    expanded::Int
end
Base.copy(x::SubSpaceIndex) = SubSpaceIndex(x.outer, x.inner, x.expanded)

struct EnsembleIndex
    outer::Int
    inner::Int
end
Base.copy(x::EnsembleIndex) = EnsembleIndex(x.outer, x.inner)

"""
    ConcreteIndexes(expected_lengths, indices)

Container that stores concrete ensemble indices for ordered atoms.
`expected_lengths` must match the number of tracked ensemble subspaces, and each
entry in `indices` must have the corresponding length.
"""
struct ConcreteIndexes
    expected_lengths::Vector{Int}
    indices::Vector{Vector{Int}}
    function ConcreteIndexes(expected_lengths::Vector{Int}, indices::Vector{Vector{Int}})
        length(expected_lengths) == length(indices) ||
            error("ConcreteIndexes: expected $(length(expected_lengths)) ensemble entries, got $(length(indices)).")
        copied_expected = copy(expected_lengths)
        copied_indices = Vector{Vector{Int}}(undef, length(indices))
        @inbounds for i in eachindex(indices)
            curr = copy(indices[i])
            length(curr) == copied_expected[i] ||
                error("ConcreteIndexes: ensemble $(i) expects $(copied_expected[i]) entries, got $(length(curr)).")
            copied_indices[i] = curr
        end
        return new(copied_expected, copied_indices)
    end
end
ConcreteIndexes(expected_lengths::AbstractVector{<:Integer}) = begin
    lengths = Vector{Int}(expected_lengths)
    return ConcreteIndexes(lengths, [fill(0, lengths[i]) for i in eachindex(lengths)])
end

"""
    pushindex!(ci::ConcreteIndexes, idx::EnsembleIndex, value::Int)
    pushindex!(ci::ConcreteIndexes, updates...)
    pushindex!(ci::ConcreteIndexes, collection)

Store integer values for specific ensemble slots within `ci`. The primary form
accepts an `EnsembleIndex` alongside the value to record. Additional overloads
allow applying multiple updates either via tuples/pairs or any iterable of such
pairs. Each update overwrites the targeted `(outer, inner)` entry.
"""
function pushindex!(ci::ConcreteIndexes, idx::EnsembleIndex, value::Integer)
    1 ≤ idx.outer ≤ length(ci.indices) || throw(ArgumentError("Ensemble index $(idx.outer) out of bounds (expected 1:$(length(ci.indices)))."))
    entries = ci.indices[idx.outer]
    1 ≤ idx.inner ≤ length(entries) || throw(ArgumentError("Inner index $(idx.inner) out of bounds for ensemble $(idx.outer) (expected 1:$(length(entries)))."))
    entries[idx.inner] = Int(value)
    return ci
end
@inline pushindex!(ci::ConcreteIndexes, item::Tuple{EnsembleIndex,<:Integer}) = pushindex!(ci, item[1], item[2])
@inline pushindex!(ci::ConcreteIndexes, item::Pair{EnsembleIndex,<:Integer}) = pushindex!(ci, first(item), last(item))
function pushindex!(ci::ConcreteIndexes, items::AbstractVector{T}) where {T}
    for item in items
        pushindex!(ci, item)
    end
    return ci
end
function pushindex!(ci::ConcreteIndexes, items::Vararg{Union{Tuple{EnsembleIndex,<:Integer},Pair{EnsembleIndex,<:Integer}}})
    for item in items
        pushindex!(ci, item)
    end
    return ci
end

Base.getindex(ci::ConcreteIndexes, i::Int) = ci.indices[i]
Base.length(ci::ConcreteIndexes) = length(ci.indices)
Base.iterate(ci::ConcreteIndexes) = iterate(ci.indices)
Base.iterate(ci::ConcreteIndexes, state) = iterate(ci.indices, state)


function vecvec_or(A::AbstractVector{<:AbstractVector{Bool}}, B::AbstractVector{<:AbstractVector{Bool}})
    out = Vector{BitVector}(undef, length(A))
    @inbounds for i in eachindex(B)
        ai = A[i]; bi = B[i]
        n = length(bi)  # result has the same length as B[i]
        outi = BitVector(undef, n)
        @inbounds @simd for j in 1:n
            outi[j] = ai[j] | bi[j]
        end
        out[i] = outi
    end
    @inbounds for i in length(B)+1:length(A)
        out[i] = A[i]   
    end
    return out
end
function vecvec_or!(A::AbstractVector{<:AbstractVector{Bool}}, B::AbstractVector{<:AbstractVector{Bool}})
    # length(A) >= length(B) || throw(DimensionMismatch("B is longer than A ($(length(B)) > $(length(A)))"))
    @inbounds for i in eachindex(B)        # only iterate existing B[i]
        ai = A[i]; bi = B[i]
        @inbounds @simd for j in eachindex(ai, bi)  # up to length(bi)
            ai[j] |= bi[j]
        end
    end
    return B
end


function unique_sorted!(v::Vector{T}) where {T}
    n = length(v)
    n ≤ 1 && return v
    w = 1
    @inbounds for i in 2:n
        if v[i] != v[w]
            w += 1
            v[w] = v[i]
        end
    end
    return resize!(v, w)
end

function sort_unique!(v::Vector{T})::Vector{T} where {T}
    return unique_sorted!(sort!(v))
end

function variants_C(name::String)
    clean = strip(name)

    # Determine core
    if startswith(clean, "C")
        core = clean[2:end]   # drop the first C
    else
        core = clean
    end

    core = lowercase(core)

    cname = "C" * uppercasefirst(core)
    name_cap = uppercasefirst(core)
    name_low = lowercase(core)

    return (cname, name_cap, name_low)
end

# find first variant, that first searches to the right of a vector, and only if nothing is found searches to the left 
function findfirstfreeafterbefore(x::BitVector, start_index::Int)::Union{Int, Nothing}
    n = length(x) 
    @assert 1 <= start_index <= n "Starting index must be in bounds."
    # Search to the right first, then left
    @inbounds for i in Base.Iterators.flatten((start_index+1:n, 1:start_index-1)) 
        if !x[i]
            return i 
        end
    end
    return nothing
end

module SparsePermutationTools
using SparseArrays
export SparsePermutation, sparseperm, identityperm, denseperm, perm_image, applyperm, applyperm!, composeperm, as_repartition_moves

struct SparsePermutation <: AbstractVector{Int}
    len::Int
    delta::SparseVector{Int,Int}
    function SparsePermutation(len::Int, delta::SparseVector{Int,Int})
        len >= 0 || throw(ArgumentError("length must be non-negative"))
        length(delta) == len || throw(DimensionMismatch("delta length $(length(delta)) does not match permutation length $len"))
        for (idx, shift) in zip(delta.nzind, delta.nzval)
            1 <= idx <= len || throw(ArgumentError("source index out of bounds"))
            dest = idx + shift
            1 <= dest <= len || throw(ArgumentError("target index out of bounds"))
        end
        return new(len, delta)
    end
end

SparsePermutation(len::Int) = SparsePermutation(len, spzeros(Int, len))

function sparseperm(perm::AbstractVector{<:Int})
    n = length(perm)
    idxs = Int[]
    vals = Int[]
    for (i, val) in enumerate(perm)
        1 <= val <= n || throw(ArgumentError("value out of bounds"))
        diff = val - i
        diff == 0 && continue
        push!(idxs, i)
        push!(vals, diff)
    end
    return SparsePermutation(n, sparsevec(idxs, vals, n))
end

identityperm(n::Int) = SparsePermutation(n)

function denseperm(sp::SparsePermutation)
    perm = collect(1:sp.len)
    perm .+= sp.delta
    return perm
end

function perm_image(sp::SparsePermutation, idx::Int)
    1 <= idx <= sp.len || throw(BoundsError(sp, idx))
    return idx + sp.delta[idx]
end

function applyperm(sp::SparsePermutation, data::AbstractVector)
    length(data) == sp.len || throw(DimensionMismatch("expected length $(sp.len), got $(length(data))"))
    result = similar(data)
    copyto!(result, data)
    for i in 1:sp.len
        dest = perm_image(sp, i)
        result[dest] = data[i]
    end
    return result
end

function applyperm!(data::AbstractVector, sp::SparsePermutation)
    tmp = applyperm(sp, data)
    copyto!(data, tmp)
    return data
end

function composeperm(p::SparsePermutation, q::SparsePermutation)
    p.len == q.len || throw(DimensionMismatch("permutations act on different lengths"))
    affected = Set{Int}()
    union!(affected, q.delta.nzind)
    union!(affected, p.delta.nzind)
    idxs = Int[]
    vals = Int[]
    for idx in affected
        img = perm_image(p, perm_image(q, idx))
        diff = img - idx
        diff == 0 && continue
        push!(idxs, idx)
        push!(vals, diff)
    end
    return SparsePermutation(p.len, sparsevec(idxs, vals, p.len))
end

function composeperm(first::SparsePermutation, rest::SparsePermutation...)
    acc = first
    for perm in rest
        acc = composeperm(acc, perm)
    end
    return acc
end

function as_repartition_moves(sp::SparsePermutation)::Vector{Tuple{Int,Int}}
    moves = Tuple{Int,Int}[]
    for idx in sp.delta.nzind
        push!(moves, (idx, idx + sp.delta[idx]))
    end
    return moves
end

import Base: length, copy, size, axes, IndexStyle, iterate

length(sp::SparsePermutation) = sp.len
size(sp::SparsePermutation) = (sp.len,)
axes(sp::SparsePermutation) = (Base.OneTo(sp.len),)
IndexStyle(::Type{SparsePermutation}) = IndexLinear()

function iterate(sp::SparsePermutation, state::Int=1)
    state > sp.len && return nothing
    return (sp[state], state + 1)
end

function copy(sp::SparsePermutation)
    return SparsePermutation(sp.len, copy(sp.delta))
end

function Base.getindex(sp::SparsePermutation, idx::Int)
    1 <= idx <= sp.len || throw(BoundsError(sp, idx))
    return idx + sp.delta[idx]
end

end
