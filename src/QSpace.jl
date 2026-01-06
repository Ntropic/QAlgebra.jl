module QSpaces

using ComplexRationals
using ..CFunctions
import ..CFunctions: CParticle, CAtom, ParameterValues, AbstractIndexMode, ParameterInfo, update_t!, resolve_param!
import ..ConcreteIndexes
using ..StringUtils
using ..Cumulants: ReducedCumulantList
using ..Sampler
using ..Sampler: AbstractEnsembleSample, DiscreteSamples, ContinuousSamples
using ..ParameterGroups: ParameterGroup, ParameterGroupKind, ParameterGroupDistribution, ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction, WhereWhichParamGroup
using Base: WeakRef, GC
using SparseArrays
using ..QAlgebra: PRINT_NON_ENSEMBLE_INDEXES

export OperatorSet, operator_magnitude, max_operator_magnitude, SubSpaceDicts, AbstractOperatorDicts
export CR_ZERO, CR_ONE
export Ensemble, SubSpace, SubSpaceDefinitions, SubSpaceInfo, SubSpaceIndex, EnsembleIndex, outer, inner, expanded, Index2Symbol, Index2String, Index2Ensemble, Index2Ensemble_and_Summation, SummationIndex2SubSpaceIndex, SubSpaceIndex2EnsembleIndex, pushindex!
export AbstractEnsembleSample, DiscreteSamples, ContinuousSamples
export OperatorType, OperatorTypeInfo, OperatorDefinitions
export Parameter, ParameterDefinitions, set_parameter_group_definition!, param2string, default_group_signature
export QSpace
export AbstractIndexParameters

Is = Vector{Int}
const CR_ZERO = ComplexRational(0,0,1)
const CR_ONE  = ComplexRational(1,0,1)

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
    len::Int                # length of indices describing operator
    neutral_element::Vector{Int}   # neutral element of the operator set
    base_ops::Vector{Vector{Int}}
    ops::Vector{String}     # operator symbols
    op_product::Function    # takes operator indices of two operators of this set and outputs a vector of tuples of coefficients and associated indices for the resulting operators in this set
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
    add_index = PRINT_NON_ENSEMBLE_INDEXES
    for curr_ind in os.base_ops
        push!(op_strs, os.op2str(curr_ind, "p"; add_index=add_index))
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
include("QSpaceOps/QSpaceDicts.jl")
using ..Sampler: build_discrete_samples, build_continuous_samples

"""
    QSpace(subspace_def, op_def, param_def)

Create the full working space that ties together subspaces, operator definitions, and
parameter families. A `QSpace` keeps:
- subspace topology (`SubSpaceDefinitions`) plus the induced ensembles,
- operator information (`OperatorDefinitions`) for constructing concrete atoms, and
- parameter collections (`ParameterDefinitions`) together with their current values.
"""
mutable struct QSpace
    # Subspace definitions:
    subspaces::Vector{SubSpace}
    subspace_info::SubSpaceInfo    # Info object containing references to all the indexing of outer and inner subspaces
    ensembles::Vector{Ensemble}
    of_time::Bool

    # Abstract operators
    operatortypes::Vector{OperatorType}
    operatortype_info::OperatorTypeInfo

    # Parameter fields:
    param_info::ParameterInfo
    where_which_param_groups::WhereWhichParamGroup
    sample_index_param_values::ParameterValues
    parameter_dicts::ParameterDicts
    subspace_dicts::SubSpaceDicts
    operator_dicts::AbstractOperatorDicts

    c_one::CAtom                            # onelike function in CFunctions 
    c_zero::CAtom                           # zerolike function in CFunctions 
    cumulant_cache::ReducedCumulantList

    function QSpace(subspace_def::SubSpaceDefinitions, op_def::OperatorDefinitions, param_def::ParameterDefinitions; of_time::Bool=false)
        subspace_def = deepcopy(subspace_def)
        op_def = deepcopy(op_def)
        param_def = deepcopy(param_def)
        # ==========> 1st Subspaces <==========
        subspaces = subspace_def.subspaces
        subspace_info = SubSpaceInfo(subspaces; of_time=of_time)
        used_symbols = subspace_def.used_symbols


        # ==========> 2nd Abstract Operators <==========
        operatortypes = OperatorDefinitions2OperatorType(op_def, subspace_def)
        operatortype_info = OperatorTypeInfo(operatortypes, commute_fun=op_def.commute_fun, check_n=op_def.check_n) 

        # ==========> 3rd Parameters <==========
        param_info, sample_index_param_values, parameter_dicts = ParameterDefinitions2Parameters(param_def, subspace_info, subspaces, used_symbols)
        where_which = WhereWhichParamGroup(param_info.param_groups)

        subspace_dicts = build_subspace_dicts(subspaces)
        operator_dicts = build_operator_dicts(operatortypes)

        # Generate the string representations
        c_zero = CAtom(param_info, CR_ZERO, CParticle{AbstractIndex}[])
        c_one = CAtom(param_info, CR_ONE, CParticle{AbstractIndex}[])
        cumulant_cache = ReducedCumulantList(1)
        ensembles = Ensemble[ss.ensemble for ss in subspaces if ss.ensemble !== nothing]

        qss = new(subspaces, subspace_info, ensembles, of_time,                  # Subspaces 
                operatortypes, operatortype_info,                                 # Abstract Operators 
                param_info, where_which, sample_index_param_values, parameter_dicts, subspace_dicts, operator_dicts,
                c_one, c_zero, cumulant_cache)    # Precomputed operator blueprints 

        
        GC.@preserve qss begin
            for ens in ensembles
                ens.qspace_ref = WeakRef(qss)
            end
        end

        return qss
    end

# close mutable struct QSpace
end

function Base.getproperty(qspace::QSpace, s::Symbol)
    if s === :max_t_ind
        return qspace.sample_index_param_values.max_t_ind
    else
        return getfield(qspace, s)
    end
end

# Define the custom show for QSpace.
function Base.show(io::IO, qspace::QSpace)
    groups = qspace.param_info.param_groups
    group_labels = [_format_group_signature(group) for group in groups]
    if get(io, :compact, false)
        param_str = join(group_labels, ",")
        subs = [join(subspace_labels(ss), ",") for ss in qspace.subspaces]
        ops  = string.(qspace.operatortypes)
        print(io, "QSpace([", param_str, "], sub=", subs, ", ops=", ops, ")")
        return
    end
    # Header line
    param_str = join(group_labels, ",")
    println(io, "QSpace: [" * param_str * "]")

    # Build LHS and RHS strings for each subspace
    lhs_list = String[]
    rhs_list = Any[]
    for ss in qspace.subspaces
        prefix = ss.is_ensemble_ss ? "Ensemble: " : "Subspace: "
        labels = subspace_labels(ss)
        lhs = prefix * labels[1]
        if has_summation(ss)
            lhs *= ", ∑ " * labels[2]
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

# ------------------------------------------------------------------
# ConcreteIndexes convenience constructors
@inline function _ensemble_expected_lengths(qspace::QSpace)::Vector{Int}
    lengths = Int[]
    for ss in qspace.subspaces
        ss.is_ensemble_ss || continue
        push!(lengths, has_summation(ss) ? 2 : 1)
    end
    return lengths
end
@inline ConcreteIndexes(qspace::QSpace) = ConcreteIndexes(_ensemble_expected_lengths(qspace))
@inline function ConcreteIndexes(qspace::QSpace, indices::AbstractVector{<:AbstractVector{<:Integer}})
    expected = _ensemble_expected_lengths(qspace)
    vectors = [Vector{Int}(idx) for idx in indices]
    return ConcreteIndexes(expected, vectors)
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

"""
    AbstractIndexParameters(qspace::QSpace)

Construct `ParameterValues{AbstractIndexMode}` for the given `qspace`. 
Allows evaluating CFunctions using specific parameter values mapped to the abstract ensemble subspace indices.
"""
function AbstractIndexParameters(qspace::QSpace)::ParameterValues
    subspaces = qspace.subspaces
    info = qspace.param_info
    maps = Vector{Vector{Int}}(undef, length(subspaces))
    @inbounds for (idx, ss) in enumerate(subspaces)
        if ss.ensemble === nothing
            maps[idx] = Int[]
        else
            maps[idx] = copy(ss.ensemble.distribution_group_indices)
        end
    end
    return ParameterValues(info; mode=AbstractIndexMode(), max_t_ind=qspace.sample_index_param_values.max_t_ind, ensemble_distribution_groups=maps)
end

end # module QSpaces
