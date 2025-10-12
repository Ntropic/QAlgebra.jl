"""
    assign_ensemble_samples!(subspaces, param_values)

Inspect every ensemble in `subspaces` and (re)build sampling data for those
whose distribution-backed parameter groups possess payloads. Samplers are
instantiated only when *every* group registered with the ensemble has an
assigned `QDistribution`; otherwise the sampler is left unset so dependent code
does not observe stale samples.

The function is invoked during `QSpace` construction and again whenever
`resolve_param!` updates a distribution payload.
"""
function assign_ensemble_samples!(groups::AbstractVector{ParameterGroupLike}, subspace_info::SubSpaceInfo, param_values::ParameterValues)
    info = param_values.param_info
    outer_symbols = info.outer_labels_symbols
    outer_names = info.outer_labels
    for ss in subspaces
        ens = ss.ensemble
        ens === nothing && continue
        dist_idxs = ens.distribution_group_indices
        isempty(dist_idxs) && continue
        all_assigned = true
        dists = Vector{QDistribution}(undef, length(dist_idxs))
        for (pos, idx) in enumerate(dist_idxs)
            payload = groups[idx].payload
            if payload isa QDistribution
                dists[pos] = payload
            else
                all_assigned = false
                break
            end
        end
        if all_assigned
            group_symbols = outer_symbols[dist_idxs]
            group_names = outer_names[dist_idxs]
            method = ens.sample_method === :default ?
                (ens.as_continuum ? :chebychev : :random) : ens.sample_method
            sample = if ens.as_continuum
                build_continuous_samples(ens, dist_idxs, group_symbols, group_names, dists;
                    method=method)
            else
                build_discrete_samples(ens, dist_idxs, group_symbols, group_names, dists;
                    method=method,
                    num_nodes=ens.sample_num_nodes,
                    atol=ens.sample_atol,
                    rtol=ens.sample_rtol,
                    max_iter=ens.sample_max_iter)
            end
            ens.sampler = sample::AbstractEnsembleSample
            attach_samples!(param_values, sample)
            sample_count = size(sample.samples, 1)
            register_ensemble_sample_size!(param_values, ss.ss_outer_ind, sample_count)
        else
            ens.sampler = nothing
        end
    end
    ensure_functions!(param_values)
    return nothing
end

function update_t!(qspace::QSpace, value::Float64; slot::Int=0)
    update_t!(qspace.sample_index_param_values, value; slot=slot)
    return qspace
end

"""
    resolve_param!(qspace, name, payload)

Attach or update the payload for a parameter group on an existing `qspace`. This could be a a scalar value or a Function. 
    The function checks if the argument is in line with the requirements of the parameter group. 
"""
#function resolve_param!(qspace::QSpace, name::Union{Symbol,String}, payload::ParameterGroupStorageUnion)
#    group_idx = get_parameter_group(qspace, name)
#    return _resolve_param!(qspace, group_idx, payload)
#end
using ..ParameterGroups: ParameterGroupLike
