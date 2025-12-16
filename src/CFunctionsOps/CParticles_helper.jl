Base.copy(p::CParticle{T}) where {T<:QIndex} = CParticle{T}(p.group_index, p.exponent, copy(p.abstract_indices), p.time_index)

@inline _with_exponent(p::CParticle{T}, exponent::Int) where {T<:QIndex} = CParticle{T}(p.group_index, exponent, p.abstract_indices, p.time_index)

@inline function Base.isequal(a::CParticle{T}, b::CParticle{T}) where {T<:QIndex}
    a.group_index == b.group_index &&
    isequal(a.time_index, b.time_index) &&
    isequal(a.abstract_indices, b.abstract_indices)
end

# strict lexicographic order on vectors of QIndex
@inline function _lexcmp_idxs(a::Vector{<:QIndex}, b::Vector{<:QIndex})
    na = length(a); nb = length(b)
    n = ifelse(na < nb, na, nb)
    @inbounds for i in 1:n
        ai = a[i]; bi = b[i]
        if !isequal(ai, bi)
            return isless(ai, bi) ? -1 : 1
        end
    end
    return na == nb ? 0 : (na < nb ? -1 : 1)
end

@inline function Base.isless(a::CParticle{T}, b::CParticle{T}) where {T<:QIndex}
    # 1) primary: group_index
    a.group_index != b.group_index && return a.group_index < b.group_index

    # 2) secondary: abstract_indices (lexicographic)
    c = _lexcmp_idxs(a.abstract_indices, b.abstract_indices)
    c != 0 && return c < 0

    # 3) tertiary: time_index
    return isless(a.time_index, b.time_index)
end


function append(p::Vector{CParticle{T}}, q::Vector{CParticle{T}}) where {T<:QIndex}
    nd, ns = length(p), length(q)
    nd == 0 && return [copy(b) for b in q if b.exponent != 0]
    ns == 0 && return copy(p)

    out = Vector{CParticle{T}}(undef, nd + ns)
    i = 1; j = 1; k = 0
    @inbounds while i <= nd && j <= ns
        a = p[i]; b = q[j]
        if isequal(a, b)
            s = a.exponent + b.exponent
            if s != 0
                k += 1
                out[k] = _with_exponent(a, s)
            end
            i += 1; j += 1
        elseif isless(a, b)
            k += 1
            out[k] = a
            i += 1
        else
            if b.exponent != 0
                k += 1
                out[k] = copy(b)
            end
            j += 1
        end
    end
    @inbounds while i <= nd
        k += 1
        out[k] = p[i]
        i += 1
    end
    @inbounds while j <= ns
        b = q[j]
        if b.exponent != 0
            k += 1
            out[k] = copy(b)
        end
        j += 1
    end
    resize!(out, k)
    return out
end

function normalise_particles!(particles::Vector{CParticle{T}}) where {T<:QIndex}
    filter!(p -> p.exponent != 0, particles)
    return particles
end
