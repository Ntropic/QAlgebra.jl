module QSpaces

using ComplexRationals
using ..CFunctions
import ..CFunctions: update_t!, resolve_param!
using ..StringUtils
using ..Cumulants: ReducedCumulantList
using ..Sampler
using ..Sampler: AbstractEnsembleSample, DiscreteSamples, ContinuousSamples
using ..ParameterGroups: ParameterGroup, ParameterGroupKind, ParameterGroupDistribution, ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction, WhereWhichParamGroup
using Base: WeakRef, GC
using SparseArrays

export OperatorSet, operator_magnitude, max_operator_magnitude, SubSpaceDicts, AbstractOperatorDicts
export Ensemble, SubSpace, SubSpaceDefinitions, SubSpaceInfo, SubSpaceIndex, outer, inner, expanded, Index2Symbol, Index2String, Index2Ensemble, Index2Ensemble_and_Summation, SummationIndex2SubSpaceIndex
export AbstractEnsembleSample, DiscreteSamples, ContinuousSamples
export OperatorType, OperatorTypeInfo, OperatorDefinitions
export Parameter, ParameterDefinitions, set_parameter_group_definition!, map_by_subspace, map_by_tindex
export QSpace

Is = Vector{Int}

"""
    OperatorSet(name, particle_type, len, neutral_element, base_ops, ops, op_product, op_dag, op2str, op2latex; kwargs...)

Describe the algebra associated with a family of operators. An `OperatorSet`
encodes multiplication, adjoint, and formatting behaviour used by subspaces and
ensembles.

Arguments:
- `name::String`: Human-readable identifier.
- `particle_type::String`: `fermion`, `boson`, `anyon`, … used for metadata.
- `len::Int`: Number of indices describing each operator.
- `neutral_element::Vector{Int}`: Index tuple representing the identity element.
- `base_ops::Vector{Vector{Int}}`: Basis operators (one per generator).
- `ops::Vector{String}`: Printable symbols corresponding to `base_ops`.
- `op_product::Function`: Binary product rule returning coefficient/index tuples.
- `op_dag::Function`: Adjoint involution on index tuples.
- `op2str` / `op2latex`: Formatting callbacks for console and LaTeX output.

Keyword arguments:
- `commutes::Union{Nothing,Function} = nothing`: Custom commutativity test.
- `operator_magnitude::Union{Nothing,Function} = nothing`: Ranking function for operator ordering.
- `min_ints::Union{Nothing,Vector{Int}} = nothing`: Component-wise minima for valid indices.
- `max_ints::Union{Nothing,Vector{Int}} = nothing`: Component-wise maxima (`-1` denotes unbounded).
- `max_magnitude::Int = -1`: Cached upper bound from `operator_magnitude` (`-1` lets it be inferred).
"""
struct OperatorSet
    name::String
    particle_type::String   # fermion, boson, anyon...
    len::Int                # length of indexes describing operator
    neutral_element::Vector{Int}   # neutral element of the operator set
    base_ops::Vector{Vector{Int}}
    ops::Vector{String}     # operator symbols
    op_product::Function    # takes operator indexes of two operators of this set and outputs a vector of tuples of coefficients and associated indexes for the resulting operators in this set
    op_dag::Function        # Create Complex Transpoose Conjugate
    op2str::Function        # transforms an operator index into a string for console printing
    op2latex::Function      # transforms an operator index into a LaTeX string for formatted LaTeXStrings
    commutes::Function
    operator_magnitude::Function
    min_ints::Vector{Int}   # component-wise minimum integer index
    max_ints::Vector{Int}   # component-wise maximum index (-1 marks unbounded)
    max_magnitude::Int       # -1 denotes unknown cap
    function OperatorSet(name::String, particle_type::String, len::Int, neutral_element::Vector{Int}, base_ops::Vector{Vector{Int}}, ops::Vector{String}, op_product::Function, op_dag::Function, op2str::Function, op2latex::Function; commutes::Union{Nothing,Function}=nothing, operator_magnitude::Union{Nothing,Function}=nothing, min_ints::Union{Nothing,Vector{Int}}=nothing, max_ints::Union{Nothing,Vector{Int}}=nothing, max_magnitude::Int=-1)
        length(neutral_element) == len || error("neutral_element length must match len")
        for op in base_ops
            length(op) == len || error("base_ops entries must match len")
        end
        extrema_vectors = Vector{Vector{Int}}()
        push!(extrema_vectors, neutral_element)
        append!(extrema_vectors, base_ops)
        comp_min = fill(typemax(Int), len)
        comp_max = fill(typemin(Int), len)
        for vec in extrema_vectors
            for i in 1:len
                vi = vec[i]
                if vi < comp_min[i]
                    comp_min[i] = vi
                end
                if vi > comp_max[i]
                    comp_max[i] = vi
                end
            end
        end
        inferred_min = isnothing(min_ints) ? comp_min : copy(min_ints)
        inferred_max = isnothing(max_ints) ? comp_max : copy(max_ints)
        length(inferred_min) == len || error("min_ints length must match len")
        length(inferred_max) == len || error("max_ints length must match len")
        for i in 1:len
            inferred_max[i] == -1 && continue
            inferred_max[i] >= inferred_min[i] || error("max_ints must be >= min_ints (or -1 for unbounded)")
        end
        default_commutes = let neutral = neutral_element, op_product = op_product
            function commutes_default(op1::Vector{Int}, op2::Vector{Int})::Bool
                if op1 == op2 || op1 == neutral || op2 == neutral
                    return true
                end
                prod_1 = op_product(op1, op2)
                prod_2 = op_product(op2, op1)
                if length(prod_1) != length(prod_2)
                    return false
                end
                sort!(prod_1, by=x -> x[2])
                sort!(prod_2, by=x -> x[2])
                for k in eachindex(prod_1)
                    if prod_1[k][2] != prod_2[k][2] || prod_1[k][1] != -prod_2[k][1]
                        return false
                    end
                end
                return true
            end
            commutes_default
        end
        default_magnitude = let neutral = neutral_element
            function magnitude_default(op::Is)::Int
                return op == neutral ? 0 : 1
            end
            magnitude_default
        end
        commutes_fun = isnothing(commutes) ? default_commutes : commutes
        magnitude_fun = isnothing(operator_magnitude) ? default_magnitude : operator_magnitude
        max_mag_int = Int(max_magnitude)
        if max_mag_int == -1 && all(max_val != -1 for max_val in inferred_max)
            max_mag_int = magnitude_fun(copy(inferred_max))
        end
        return new(name, particle_type, len, neutral_element, base_ops, ops, op_product, op_dag, op2str, op2latex, commutes_fun, magnitude_fun, inferred_min, inferred_max, max_mag_int)
    end
end
#function OperatorSet
function Base.show(io::IO, os::OperatorSet)
    op_strs = String[]
    for curr_ind in os.base_ops
        push!(op_strs, os.op2str(curr_ind, "p"))
    end
    print(io, os.name, " (", os.particle_type, "):  " * join(op_strs, ","))
end

operator_magnitude(os::OperatorSet, op::Is)::Int = os.operator_magnitude(op)

function max_operator_magnitude(os::OperatorSet)::Int
    os.max_magnitude != -1 && return os.max_magnitude
    any(==( -1), os.max_ints) && return -1
    return operator_magnitude(os, copy(os.max_ints))
end
include("OperatorSets/Qubit_Pauli.jl")
include("OperatorSets/Qubit_PM.jl")
include("OperatorSets/Ladder.jl")

include("QSpaceOps/QSpace_subspaces.jl")
include("QSpaceOps/QSpace_abstract.jl")
include("QSpaceOps/QSpace_parameters.jl")
using ..Sampler: build_discrete_samples, build_continuous_samples

struct SubSpaceDicts
    by_outer::Dict{Symbol,Int}
    by_inner::Dict{Symbol,Tuple{Int,Int}}
end
struct AbstractOperatorDicts
    by_name::Dict{Symbol,Int}
end
function build_subspace_dicts(subspaces::Vector{SubSpace})::SubSpaceDicts
    outer_map = Dict{Symbol,Int}()
    inner_map = Dict{Symbol,Tuple{Int,Int}}()
    for (idx, ss) in enumerate(subspaces)
        if haskey(outer_map, ss.key_symbol)
            error("Duplicate outer subspace key $(ss.key_symbol) detected while building QSpace.")
        end
        outer_map[ss.key_symbol] = idx
        for (inner_idx, sym) in enumerate(ss.keys_symbols)
            if haskey(inner_map, sym)
                error("Duplicate inner subspace key $(sym) detected while building QSpace.")
            end
            inner_map[sym] = (idx, inner_idx)
        end
    end
    return SubSpaceDicts(outer_map, inner_map)
end
function build_operator_dicts(operatortypes::Vector{OperatorType})::AbstractOperatorDicts
    map = Dict{Symbol,Int}()
    for (idx, optype) in enumerate(operatortypes)
        if haskey(map, optype.name_sym)
            error("Duplicate operator type symbol $(optype.name_sym) detected while building QSpace.")
        end
        map[optype.name_sym] = idx
    end
    return AbstractOperatorDicts(map)
end

"""
    QSpace(subspace_def, op_def, param_def; max_t_ind=0)

Create the full working space that ties together subspaces, operator definitions, and
parameter families. A `QSpace` keeps:
- subspace topology (`SubSpaceDefinitions`) plus the induced ensembles,
- operator information (`OperatorDefinitions`) for constructing concrete atoms, and
- parameter collections (`ParameterDefinitions`) together with their current values.

Keyword arguments:
- `max_t_ind::Int = 0`: Highest time index admitted when constructing time-dependent expressions.
"""
mutable struct QSpace
    # Subspace definitions:
    subspaces::Vector{SubSpace}
    subspace_info::SubSpaceInfo    # Info object containing references to all the indexing of outer and inner subspaces
    ensembles::Vector{Ensemble}

    # Abstract operators
    operatortypes::Vector{OperatorType}
    operatortype_info::OperatorTypeInfo

    # Parameter fields:
    params::Vector{Parameter}
    param_info::ParameterInfo
    where_which_param_groups::WhereWhichParamGroup
    sample_index_param_values::ParameterValues
    parameter_dicts::ParameterDicts
    subspace_dicts::SubSpaceDicts
    operator_dicts::AbstractOperatorDicts

    I_op::Vector{Is}               # Neutral Vector of all expanded subspaces
    I_ensemble_op::Vector{Vector{Is}}      # Neutral Vector of all expanded ensemble subspaces
    c_one::CAtom                            # onelike function in CFunctions 
    c_zero::CAtom                           # zerolike function in CFunctions 
    cumulant_cache::ReducedCumulantList
    max_t_ind::Int

    function QSpace(subspace_def::SubSpaceDefinitions, op_def::OperatorDefinitions, param_def::ParameterDefinitions; max_t_ind::Int=0)
        # ==========> 1st Subspaces <==========
        subspaces = subspace_def.subspaces
        subspace_info = SubSpaceInfo(subspaces)
        used_symbols = subspace_def.used_symbols
        I_op, I_ensemble_op = subspace_def.I_op, subspace_def.I_ensemble_op  # These are

        # ==========> 2nd Abstract Operators <==========
        operatortypes = OperatorDefinitions2OperatorType(op_def, subspace_def)
        operatortype_info = OperatorTypeInfo(operatortypes, commute_fun=op_def.commute_fun, check_n=op_def.check_n) 

        # ==========> 3rd Parameters <==========
        params, param_info, sample_index_param_values, parameter_dicts = ParameterDefinitions2Parameters(param_def, subspace_info, subspaces, used_symbols, max_t_ind)
        where_which = WhereWhichParamGroup(param_info.param_groups)

        subspace_dicts = build_subspace_dicts(subspaces)
        operator_dicts = build_operator_dicts(operatortypes)

        # Generate the string representations
        c_one = CAtom(param_info, spzeros(Int, length(params)))
        c_zero = CAtom(param_info, ComplexRational(0,0, 1), spzeros(Int, length(params)))
        cumulant_cache = ReducedCumulantList(1)
        ensembles = Ensemble[ss.ensemble for ss in subspaces if ss.ensemble !== nothing]

        qss = new(subspaces, subspace_info, ensembles,                           # Subspaces 
                operatortypes, operatortype_info,                                 # Abstract Operators 
                params, param_info, where_which, sample_index_param_values, parameter_dicts, subspace_dicts, operator_dicts,
                I_op, I_ensemble_op, c_one, c_zero, cumulant_cache, max_t_ind)    # Precomputed operator blueprints 

        
        GC.@preserve qss begin
            for ens in ensembles
                ens.qspace_ref = WeakRef(qss)
            end
        end

        return qss
    end
end
# Define the custom show for QSpace.
function Base.show(io::IO, qspace::QSpace)
    if get(io, :compact, false)
        param_str = if length(qspace.params) < 12
            join((p.param_str for p in qspace.params), ",")
        else
            group_labels = String[]
            for g in 1:length(qspace.param_info.param_groups)
                idx = findfirst(p -> p.group_index == g, qspace.params)
                idx === nothing && continue
                push!(group_labels, qspace.params[idx].param_str)
            end
            join(group_labels, ",")
        end
        subs = [join(ss.keys[1:ss.num_operator_indexes], ",") for ss in qspace.subspaces]
        ops  = string.(qspace.operatortypes)
        print(io, "QSpace([", param_str, "], sub=", subs, ", ops=", ops, ")")
        return
    end
    # Header line
    if length(qspace.params) < 12
        param_str = join([p.param_str for p in qspace.params], ",")
    else
        group_labels = String[]
        for g in 1:length(qspace.param_info.param_groups)
            idx = findfirst(p -> p.group_index == g, qspace.params)
            idx === nothing && continue
            push!(group_labels, qspace.params[idx].param_str)
        end
        param_str = join(group_labels, ",")
    end
    println(io, "QSpace: [" * param_str * "]")

    # Build LHS and RHS strings for each subspace
    lhs_list = String[]
    rhs_list = Any[]
    for ss in qspace.subspaces
        prefix = ss.is_ensemble_ss ? "Ensemble: " : "Subspace: "
        op_keys = ss.keys[1:ss.num_operator_indexes]
        lhs = prefix * join(op_keys, ",")
        if ss.num_sum_indexes > 0
            sum_keys = ss.keys[ss.num_operator_indexes+1:end]
            lhs *= ", ∑ " * join(sum_keys, ",")
        end
        push!(lhs_list, lhs)
        push!(rhs_list, ss.op_set)
    end

    # Find maximum lhs length
    maxlen = maximum(length, lhs_list)

    # Print aligned
    for (lhs, rhs) in zip(lhs_list, rhs_list)
        padded_lhs = rpad(lhs, maxlen)
        print(io, "   - ", padded_lhs, " → ")
        show(io, rhs)
        println(io)
    end
    # Operator types
    op_strs = String[operator_type2string(op) for op in qspace.operatortypes]
    if length(op_strs) > 0 
        println(io, "   - Abstract Ops:", join(op_strs, ", "))
    end
end

include("QSpaceOps/QSpace_get_types.jl")
using ..ParameterGroups: ParameterGroupStorageUnion
""" 
    update_t!(qspace::QSpace, value::Float64; slot::Int=0)
Change the time for one of the time parameters, selected via slot. 
""" 
function update_t!(qspace::QSpace, value::Float64; slot::Int=0)
    update_t!(qspace.sample_index_param_values, value; slot=slot)
    return qspace
end

"""
    resolve_param!(qspace, name, payload)

Attach or update the payload for a parameter group on an existing `qspace`. This could be a a scalar value or a Function. 
    The function checks if the argument is in line with the requirements of the parameter group. 
"""
function resolve_param!(qspace::QSpace, name::Union{Symbol,String}, payload::ParameterGroupStorageUnion)
    group_idx = get_parameter_group(qspace, name)
    return resolve_param!(qspace.sample_index_param_values, group_idx, payload)
end

end # module QSpaces
