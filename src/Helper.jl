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
            bi[j] |= ai[j]
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

function sorted_unique_push!(arr::Vector{T}, x::T) where T
    # Find insertion index with binary search
    i = searchsortedfirst(arr, x)
    # Only insert if element is not already there
    if i > length(arr) || arr[i] != x
        insert!(arr, i, x)
    end
    return arr
end