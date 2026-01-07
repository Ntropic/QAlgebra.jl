Base.copy(p::CParticle{T}) where {T<:QIndex} = CParticle{T}(p.group_index, p.exponent, copy(p.abstract_indices), p.time_index)

@inline _with_exponent(p::CParticle{T}, exponent::Int) where {T<:QIndex} = CParticle{T}(p.group_index, exponent, p.abstract_indices, p.time_index)

function normalise_particles!(particles::Vector{CParticle{T}}) where {T<:QIndex}
    filter!(p -> p.exponent != 0, particles)
    return particles
end
