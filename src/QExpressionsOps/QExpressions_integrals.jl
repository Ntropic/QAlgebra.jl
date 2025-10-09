using ..QSpaces: QSpace, SubSpaceIndex
using ..Sampler: QInterpolator, ContinuousSamples
using ..CFunctions: ParameterInfo, ParameterValues, CFunction, evaluate, set_param!
using ..QAlgebra: ConcreteIndexes

struct IntegralDimensionDescriptor
    ensemble_outer::Int
    inner_index::Int
    group_index::Int
end

function _build_integral_interpolator(qspace::QSpace,
                                      indexes::Vector{Vector{SubSpaceIndex}})::Tuple{QInterpolator,Vector{Int},Vector{Function},Vector{IntegralDimensionDescriptor}}
    param_info = qspace.param_info
    subspace_info = param_info.subspace_info
    subspace_info === nothing &&
        error("ParameterInfo is missing subspace information required for integrals.")

    node_vectors = Vector{Vector{Float64}}()
    bounds = Vector{Tuple{Float64,Float64}}()
    axis_lengths = Int[]
    pdfs = Vector{Function}()
    dim_info = IntegralDimensionDescriptor[]
    endpoints_flag = true

    @inbounds for ensemble_indexes in indexes
        for sub_idx in ensemble_indexes
            subspace = qspace.subspaces[sub_idx.outer]
            ens = subspace.ensemble
            ens === nothing && error("Subspace $(subspace.key) is not backed by an ensemble; cannot define integrals over it.")
            sample = ens.sampler
            sample === nothing && error("Ensemble $(subspace.key) has no registered sampler; ensure sampling was attached before defining integrals.")
            sample isa ContinuousSamples ||
                error("Ensemble $(subspace.key) provides discrete samples; coefficient integrals currently require continuous samplers.")
            inter = sample.interpolator
            endpoints_flag &= inter.endpoints
            @inbounds for dim_idx in 1:inter.dims
                dist = sample.distributions[dim_idx]
                push!(node_vectors, copy(inter.nodes[dim_idx]))
                push!(bounds, inter.bounds[dim_idx])
                push!(axis_lengths, length(inter.nodes[dim_idx]))
                let d = dist
                    push!(pdfs, x -> d.pdf(x))
                end
                push!(dim_info, IntegralDimensionDescriptor(sub_idx.outer, sub_idx.inner, sample.group_indices[dim_idx]))
            end
        end
    end

    isempty(node_vectors) && error("Integral definition requires at least one continuous ensemble dimension.")
    joint = QInterpolator(node_vectors, bounds; method=:custom, endpoints=endpoints_flag)
    return joint, axis_lengths, pdfs, dim_info
end

function _build_integral_assignments(param_info::ParameterInfo,
                                     dim_info::Vector{IntegralDimensionDescriptor})::Vector{Vector{Int}}
    assignments = Vector{Vector{Int}}(undef, length(dim_info))
    for (dim_idx, info) in enumerate(dim_info)
        group_params = param_info.param_groups[info.group_index].parameter_indices
        selected = Int[]
        for param_idx in group_params
            tuples = param_info.param_index_tuples[param_idx]
            any(t -> t[1] == info.ensemble_outer && t[2] == info.inner_index, tuples) || continue
            push!(selected, param_idx)
        end
        assignments[dim_idx] = unique(selected)
    end
    return assignments
end

function _make_integrand(expr::CFunction,
                         pv::ParameterValues,
                         assignments::Vector{Vector{Int}})::Function
    default_indexes = ConcreteIndexes(expr.param_info)
    function integrand(xpt::AbstractVector{<:Real})
        @inbounds for (dim_idx, val) in enumerate(xpt)
            for param_idx in assignments[dim_idx]
                set_param!(pv, param_idx, Float64(val))
            end
        end
        return ComplexF64(Complex(evaluate(expr, pv, default_indexes)))
    end
    return integrand
end
