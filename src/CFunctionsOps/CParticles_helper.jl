Base.copy(p::CParticle{T}) where {T<:QIndex} = CParticle{T}(p.group_index, p.exponent, copy(p.abstract_indices), p.time_index)

@inline _with_exponent(p::CParticle{T}, exponent::Int) where {T<:QIndex} = CParticle{T}(p.group_index, exponent, p.abstract_indices, p.time_index)

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
