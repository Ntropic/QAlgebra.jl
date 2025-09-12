function vecvec_or(A::Vector{Vector{Bool}}, B::Vector{Vector{Bool}})
    # Assume they are equally shaped. 
    return broadcast.(|, A, B)
end
function vecvec_or!(A::Vector{<:AbstractVector{Bool}}, B::Vector{<:AbstractVector{Bool}})
    # length(A) == length(B) || throw(DimensionMismatch("outer lengths differ ($(length(A)) vs $(length(B)))"))
    @inbounds for i in eachindex(A, B)
        ai, bi = A[i], B[i]
        #length(ai) == length(bi) || throw(DimensionMismatch("inner lengths differ at i=$i ($(length(ai)) vs $(length(bi)))"))
        ai .|= bi                # elementwise OR, in place
    end
    return A
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
