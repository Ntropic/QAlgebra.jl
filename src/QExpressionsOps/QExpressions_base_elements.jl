using ..QIndexes: TimeIndex
using ..StringUtils: parse_time_brace, stringparse4base_operators
using ..QSpaces: CR_ONE, CR_ZERO
include("QExpressions_base_elements/QExpressions_base_elements_accessors.jl")
# Helpers 

# Build the vector of AbstractIndex values described by the given string pieces so parameters can be fully specified.
function param_comps2AbstractIndexes(comps::Vector{String}, param_group::ParameterGroup{T}) where {T}
    abstract_indices = Vector{AbstractIndex}(undef, length(comps)) 
    @inbounds for (i, (comps, ss, es, index_pair)) in enumerate(zip(comps, param_group.subspace_indices, param_group.ensemble_indices, param_group.index_string_pairs))
        letter, number = split_index(comps)
        # find out if letter is in index_pair 
        if letter == index_pair[1]
            w = 1
        elseif letter == index_pair[2]
            w = 2
        else
            error("index letter $letter not found in index pair $index_pair.")
        end
        abstract_indices[i] = AbstractIndex(ss, es, number, w == 2)
    end
    return abstract_indices 
end
# returns the index (-1 if not found) of the subspace with name as an index name, and a Boolean to specify the index is a summation index
function find_indexname_in_subspaces(qspace::QSpace, name::String)::Tuple{Int, Bool}
    for (idx, subspace) in enumerate(qspace.subspaces)
        if subspace.is_ensemble_ss
            ensemble::Ensemble = subspace.ensemble 
            if name == ensemble.non_sum_string
                return idx, false
            elseif name == ensemble.sum_string
                return idx, true
            end
        else
            if name == subspace.key 
                return idx, false 
            end
        end
    end
    return -1, true #nothing found 
end
function find_indexname_in_ensembles(qspace::QSpace, name::String)::Tuple{Int, Bool}
        for (idx, subspace) in enumerate(qspace.subspaces)
        if subspace.is_ensemble_ss
            ensemble::Ensemble = subspace.ensemble 
            if name == ensemble.non_sum_string
                return idx, false
            elseif name == ensemble.sum_string
                return idx, true
            end
        end
    end
    return -1, true #nothing found 
end

function find_name_in_abstracts(qspace::QSpace, name::String)::Int
    for (idx, abstract) in enumerate(qspace.operatortypes)
        if abstract.name == name 
            return idx 
        end
    end
    return -1 
end


# Get Parameter
function get_parameter(qspace::QSpace, name::String)::Union{ParameterGroupLike, ParameterGroupAccessor, QExpr}
    param_info = qspace.param_info
    prefix, comps, t_spec = stringparse4base_operators(name)
    # match prefix 
    prefix_sym = Symbol(prefix) 
    # find index of match prefix_sym in param_info.params_symbols
    idx = findfirst(==(prefix_sym), param_info.params_symbols) 
    if idx === nothing
        error("prefix $prefix not found in param_info. Possibilities are $(param_info.params_symbols).") 
    end
    param_group = param_info.param_groups[idx] 
    expected_len_comps = length(param_group.indices)
    of_t = param_group.of_t
    len_curr_comps = length(comps)
    curr_of_t = t_spec != -1
    # either all indexes specified or none 
    if len_curr_comps == expected_len_comps
        if expected_len_comps == 0 && of_t && t_spec == -1
            return ParameterGroupAccessor(qspace, param_group, idx, param_group.subspace_indices,
                                          param_group.ensemble_indices, param_group.index_string_pairs, nothing)
        end
        if of_t && t_spec == -1
            t_spec = 0
        end
        # build the parameter that was specified by first building AbstractIndex Vector 
        abstract_indices = param_comps2AbstractIndexes(comps, param_group)
        curr_cparticle::CParticle{AbstractIndex} = CParticle(idx, 1,  abstract_indices, TimeIndex(t_spec)) 
        particle_vec = CParticle{AbstractIndex}[curr_cparticle]
        return QExpr(qspace, QComposite[QAtomProduct(qspace, CAtom(param_info, CR_ONE, particle_vec), QAtom[])])
    elseif length(comps) == 0 && expected_len_comps > 0 #fully unspecified indices
        if curr_of_t && !of_t
            error("Cannot specify a time index for parameter $(param_group.param_symbol) that isn't time dependent.")
        end
        fixed_time = of_t && curr_of_t ? t_spec : nothing
        return ParameterGroupAccessor(qspace, param_group, idx, param_group.subspace_indices,
                                      param_group.ensemble_indices, param_group.index_string_pairs, fixed_time)
    else
        error("Invalid parameter specification for $name: expected indices=$(expected_len_comps), of_t=$(of_t), got indices=$(len_curr_comps), of_t=$(curr_of_t), t_spec=$(t_spec).")
    end
end
function (param_group::ParameterGroup{G})(params::Vararg{T})::QExpr where {G, T<:Union{String,Symbol}}
    @assert length(params) == length(param_group.indices) +  param_group.of_t 
    abstract_indexes = param_comps2AbstractIndexes(params[1:length(param_group.indices)], param_group) 
    if param_group.of_t 
        t, t_spec = split_index(params[end])
        @assert t == "t"
    else
        t_spec = -1 
    end
    curr_cparticle::CParticle{AbstractIndex} = CParticle(idx, 1,  abstract_indices, TimeIndex(t_spec)) 
    particle_vec = CParticle{AbstractIndex}[curr_cparticle]
    return QExpr(qspace, QComposite[QAtomProduct(qspace, CAtom(param_info, CR_ONE, particle_vec), QAtom[])])
end
function Base.getindex(param_group::ParameterGroup{G}, params::Vararg{T})::QExpr where {G, T<:Union{String,Symbol}}
    return param_group(params...)   # just forward to the call
end

# Get Operator
function get_operator(qspace::QSpace, name::String)::Union{SubSpaceAccessor, AbstractOperatorAccessor, QExpr}
    # check if any of the subspaces use the name as index 
    base_name, args = brace_separate(name)
    idx, is_summation = find_indexname_in_subspaces(qspace, base_name)
    if idx > 0 && !is_summation
        @assert length(args) == 0 "Cannot specify a subspace by time (got $(args[1]))."
        subspace = qspace.subspaces[idx]
        op_set = subspace.op_set
        if length(op_set.ops) == 1 && op_set.ops[1] == ""
            ensemble_idx = qspace.subspace_info.ensemble_index_by_subspace_index[idx]
            abstract_index = AbstractIndex(idx, ensemble_idx, 0, false)
            operator = QParticle(op_set.base_ops[1], abstract_index)
            return QExpr(qspace, QComposite[QAtomProduct(qspace, qspace.c_one, QAtom[QTerm(QParticle[operator], default_time_index(qspace))])])
        end
        ensemble_idx = qspace.subspace_info.ensemble_index_by_subspace_index[idx]
        sum_label = ""
        if subspace.is_ensemble_ss && subspace.ensemble !== nothing
            sum_label = subspace.ensemble.sum_string
        end
        return SubSpaceAccessor(qspace, subspace, idx, ensemble_idx, (subspace.key, sum_label))
    else   # detect the concrete operator from string
        t_spec = parse_time_brace(args)
        op_name, index_names = normalize_underscore_indices(base_name)
        abstract_idx = find_name_in_abstracts(qspace, op_name)
        if abstract_idx > 0    # Is abstract operator! 
            curr_abstract = qspace.operatortypes[abstract_idx]
            if isempty(index_names)
                if !curr_abstract.of_time && t_spec != -1
                    error("Cannot specify the time index ($(args[1])) for abstract operator $(curr_abstract.name) that isn't time dependent.")
                end
                fixed_time = curr_abstract.of_time && t_spec != -1 ? t_spec : nothing
                return AbstractOperatorAccessor(qspace, curr_abstract, abstract_idx, fixed_time)
            end

            @assert length(index_names) == 1 "Operators must have a single index to specify the subspace" 
            index_name = index_names[1]
            main_index, sub_index = split_index(index_name)
            @assert length(main_index) == 0 "Abstract operator index must be numeric, got: $index_name instead."
            if curr_abstract.of_time 
                if t_spec == -1 
                    t_spec = 0
                end
            else
                @assert t_spec == -1 "Cannot specify the time index ($(args[1])) for abstract operator $(curr_abstract.name) that isn't time dependent."
            end 
            abstrac_op = QAbstract(curr_abstract, abstract_idx, sub_index, 1, false, t_spec) 
            return QExpr(qspace, QComposite[QAtomProduct(qspace, qspace.c_one, QAtom[abstrac_op])])
        end
        @assert length(index_names) == 1 "Operators must have a single index to specify the subspace" 
        index_name = index_names[1]
        # check if the op_name is an abstract operator 
        main_index, sub_index = split_index(index_name)
        # Or is it a concrete subspace operator
        time_index = resolve_time_index(qspace, t_spec == -1 ? nothing : t_spec)
        subspace_idx, is_summation = find_indexname_in_subspaces(qspace, main_index)
        # is op_name defined for that subspace? 
        op_set = qspace.subspaces[subspace_idx].op_set
        op_name_idx = findfirst(==(op_name), op_set.ops)
        @assert !isnothing(op_name_idx) "Unrecognized operator $op_name, subspace supports: $(op_set.ops)."
        ensemble_idx = qspace.subspace_info.ensemble_index_by_subspace_index[subspace_idx]
        abstract_index = AbstractIndex(subspace_idx, ensemble_idx, sub_index, is_summation)
        operator = QParticle(op_set.base_ops[op_name_idx], abstract_index)
        return QExpr(qspace, QComposite[QAtomProduct(qspace, qspace.c_one, QAtom[QTerm(QParticle[operator], time_index)])], Val(:nosimp))
    end
end

# Get Operators (all operators of a subspace)
function get_operators(qspace::QSpace, name::String)::Union{Vector{QExpr}, SubSpaceOperatorsAccessor}
    base_name, args = brace_separate(name)
    t_spec = parse_time_brace(args)
    prefix, comps = normalize_underscore_indices(base_name)
    subspace_idx, is_summation = find_indexname_in_subspaces(qspace, prefix)
    subspace_idx > 0 || error("No subspace named $(prefix) found in QSpace.")
    subspace = qspace.subspaces[subspace_idx]
    time_index = resolve_time_index(qspace, t_spec == -1 ? nothing : t_spec)
    ensemble_idx = qspace.subspace_info.ensemble_index_by_subspace_index[subspace_idx]
    sum_label = ""
    if subspace.is_ensemble_ss && subspace.ensemble !== nothing
        sum_label = subspace.ensemble.sum_string
    end
    index_pair = (subspace.key, sum_label)

    if subspace.is_ensemble_ss
        if length(comps) == 0
            return SubSpaceOperatorsAccessor(qspace, subspace, subspace_idx, ensemble_idx, index_pair)
        elseif length(comps) == 1
            abstract_index = AbstractIndex(subspace_idx, ensemble_idx, parse(Int, comps[1]), is_summation)
            return _subspace_ops_exprs(qspace, subspace, abstract_index, time_index)
        else
            error("Invalid subspace index specification for $name.")
        end
    else
        length(comps) == 0 || error("Non-ensemble subspace $(subspace.key) does not take an index.")
        abstract_index = AbstractIndex(subspace_idx, ensemble_idx, 0, false)
        return _subspace_ops_exprs(qspace, subspace, abstract_index, time_index)
    end
end
"""
    base_operators(qspace::QSpace, name::String; do_fun::Bool=false, by_ensemble::Bool=false)

Returns variables and/or operators in the state space `qspace`.

Time handling:
- Subspace operators use `qspace.of_time` to decide time dependence. If `false`, time indices are rejected and default to `-1`. If `true`, the default is `0` unless `name` ends with `_tN` (e.g. `"A_t3"`), in which case `time_index = N`.
- Abstract operators use `operator_type.of_time` with the same defaults.
"""
function base_operators(qspace::QSpace, name::String; do_fun::Bool=false, by_ensemble::Bool=false)
    name_base, t_spec = _parse_time_suffix(name)
    if name_base == "I"
        return _identity_expr(qspace)
    end

    result = get_parameter(qspace, name; do_fun=do_fun, name_base=name_base, t_spec=t_spec, error_if_missing=false)
    result !== nothing && return result

    result = get_operator(qspace, name; do_fun=do_fun, by_ensemble=by_ensemble, name_base=name_base, t_spec=t_spec, error_if_missing=false)
    result !== nothing && return result

    result = get_abstract(qspace, name; do_fun=do_fun, name_base=name_base, t_spec=t_spec, error_if_missing=false)
    result !== nothing && return result

    error("No variable, subspace component or abstract operator found for key='$name'.")
end

function base_operators(qspace::QSpace, name::Symbol; do_fun::Bool=false, by_ensemble::Bool=false)
    return base_operators(qspace, String(name); do_fun=do_fun, by_ensemble=by_ensemble)
end

function base_operators(qspace::QSpace, names::Vector{String}; do_fun::Bool=false, by_ensemble::Bool=false)
    return [base_operators(qspace, name; do_fun=do_fun, by_ensemble=by_ensemble) for name in names]
end

function base_operators(qspace::QSpace, names::Vector{Symbol}; do_fun::Bool=false, by_ensemble::Bool=false)
    return [base_operators(qspace, name; do_fun=do_fun, by_ensemble=by_ensemble) for name in names]
end

function Base.getindex(qspace::QSpace, key::Union{String,Symbol})
    name = key isa Symbol ? String(key) : key
    name_base, t_spec = _parse_time_suffix(name)
    if name_base == "I"
        return _identity_expr(qspace)
    end

    result = get_parameter(qspace, name; name_base=name_base, t_spec=t_spec, error_if_missing=false)
    result !== nothing && return result

    result = get_operator(qspace, name; name_base=name_base, t_spec=t_spec, error_if_missing=false)
    result !== nothing && return result

    result = get_abstract(qspace, name; name_base=name_base, t_spec=t_spec, error_if_missing=false)
    result !== nothing && return result

    error("No variable, subspace component or abstract operator found for key='$name'.")
end

function Base.getindex(qspace::QSpace, keys::AbstractVector{<:Union{String,Symbol}})
    return [qspace[key] for key in keys]
end
