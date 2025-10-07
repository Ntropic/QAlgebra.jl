module QSpaces

using ComplexRationals
using ..CFunctions
using ..StringUtils
using ..Cumulants: ReducedCumulantList
using Base: WeakRef, GC
using SparseArrays

export OperatorSet, operator_magnitude, max_operator_magnitude
export Ensemble, SubSpace, SubSpaceDefinitions, SubSpaceInfo, SubSpaceIndex, outer, inner, expanded, Index2Symbol, Index2String, Index2Ensemble, Index2Ensemble_and_Summation, SummationIndex2SubSpaceIndex
export OperatorType, OperatorTypeInfo, OperatorDefinitions
export Parameter, ParameterDefinitions, map_by_subspace, map_by_tindex
export QSpace

Is = Vector{Int}
"""
    OperatorSet(name::String, particle_type::String, len::Int, neutral_element::Vector{Int}, base_ops::Vector{Vector{Int}}, ops::Vector{String}, op_product::Function, op_dag::Function, op2str::Function, op2latex::Function; commutes::Union{Nothing,Function}=nothing, operator_magnitude::Union{Nothing,Function}=nothing, min_ints::Union{Nothing,Vector{Int}}=nothing, max_ints::Union{Nothing,Vector{Int}}=nothing, max_magnitude::Int=-1)

OperatorSets define the algebraic structure of a quantum system, defining ways to multiply and conjugate operators within the space, how to print them (both for plain and latex formatting), how to extract operators from strings.
We provide a few standard operator sets, such as QubitPauli, QubitPM and Ladder.
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

"""
    QSpace(subspace_def::SubSpaceDefinitions, op_def::OperatorDefinitions, param_def::ParameterDefinitions; max_t_ind::Int=0) -> QSpace

Constructs a combined Hilbert and Parameter space. The Hilbert space consists of different subspaces, themselves composed of different operator sets.
The Parameter space also defines the variables, that are needed to describe equations on the Hilbert space and abstract operators, that are not yet specified. 
Optionally you can also allow for multiple time dimensions, which can be useful for solving nested integrals over different time parameters.
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
        params, param_info = ParameterDefinitions2Parameters(param_def, subspace_info, subspaces, used_symbols, max_t_ind)
        final_group_count = length(param_info.outer_labels_symbols)
        for ss in subspaces
            resize!(ss.parameter_group_acting, final_group_count)
            resize!(ss.parameter_group_distribution, final_group_count)
        end
    
        # Generate the string representations
        c_one = CAtom(param_info, spzeros(Int, length(params)))
        c_zero = CAtom(param_info, ComplexRational(0,0,1), spzeros(Int, length(params)))
        cumulant_cache = ReducedCumulantList(1)
        ensembles = Ensemble[ss.ensemble for ss in subspaces if ss.ensemble !== nothing]

        qss = new( subspaces, subspace_info, ensembles,                           # Subspaces
                operatortypes, operatortype_info,                                 # Abstract Operators 
                params, param_info,                                               # Variables / Parameters
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
    # Header line
    if length(qspace.params) < 12
        param_str = join([p.param_str for p in qspace.params], ",")
    else
        max_val = maximum(qspace.param_info.param_group_by_index)
        indexes = Int[] 
        for i in 1:max_val
            push!(indexes, findfirst(==(i), qspace.param_info.param_group_by_index))
        end
        param_str = join([qspace.params[i].param_str for i in indexes], ",")
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
    for op in qspace.operatortypes
        println(io, "   - ", string(op))
    end
end


## Test 
#xi, yi, zi = base_operators("i", qs)
#I = base_operators("I", qs)
#alpha, beta = base_operators("params", qs)
function cleanup_terms(terms::Vector{Tuple{T,S}})::Vector{Tuple{T,S}} where {T<:Number,S}
    # 1) sort once by index
    sort!(terms, by = x -> x[2])
    # 2) prealloc output to worst‑case length and scan in one pass
    n = length(terms)
    T0 = typeof(terms[1][1])
    S0 = typeof(terms[1][2])
    cleaned = Vector{Tuple{T0,S0}}(undef, n)
    cnt = 0
    i = 1
    @inbounds while i ≤ n
        sumc, idx = terms[i]           # destructure once
        j = i + 1
        # inner loop: accumulate identical idx
        @inbounds while j ≤ n && terms[j][2] == idx
            sumc += terms[j][1]
            j += 1
        end

        # push nonzero
        if sumc != zero(T0)
            cnt += 1
            cleaned[cnt] = (sumc, idx)
        end

        i = j
    end
    resize!(cleaned, cnt)                # trim unused slots
    return cleaned
end
end # module QSpaces
