module QAlgebra

using Preferences
# === Default Coefficient Preferences ===
const DEFAULT_COEFF_PREFS = Dict(
    :FLIP_IF_FIRST_TERM_NEGATIVE  => true,
    :DO_BRACED => true,
    :EXPAND_CUMULANTS => false,
    :QUADGK_ATOL => 1e-9,
    :QUADGK_RTOL => 1e-7,
    :SAMPLE_ATOL => 1e-9,
    :SAMPLE_RTOL => 1e-7,
    )

"""
    get_default(name::Symbol)

Fetch the persisted default for the printing/simplification preference identified by
`name`. Falls back to the hard-coded setting in `DEFAULT_COEFF_PREFS` if no user
override was stored via `Preferences.jl`.
"""
function get_default(name::Symbol)
    return @load_preference(String(name), DEFAULT_COEFF_PREFS[name])
end
function set_default(name::Symbol, value)
    @set_preferences!(String(name) => value)
end

# --- Module-level global variables (not const) ---
# These are initialized immediately at module load time
FLIP_IF_FIRST_TERM_NEGATIVE  = get_default(:FLIP_IF_FIRST_TERM_NEGATIVE )
DO_BRACED = get_default(:DO_BRACED)
EXPAND_CUMULANTS = get_default(:EXPAND_CUMULANTS)
QUADGK_ATOL = get_default(:QUADGK_ATOL)
QUADGK_RTOL = get_default(:QUADGK_RTOL)
SAMPLE_ATOL = get_default(:SAMPLE_ATOL)
SAMPLE_RTOL = get_default(:SAMPLE_RTOL)

@inline function _update_pref!(name::Symbol, value)
    set_default(name, value)
    @eval $(Symbol(name)) = $value
end

"""
    set_flip_if_first_term_negative(mode::Bool)

Persist the behaviour for handling leading negative coefficients when pretty-printing.
`true` flips the overall sign so the first term is positive; `false` keeps the original
ordering even if the leading term is negative.
"""
function set_flip_if_first_term_negative(mode::Bool)
    _update_pref!(:FLIP_IF_FIRST_TERM_NEGATIVE, mode)
end

"""
    set_do_braced(mode::Bool)

Toggle whether printed expressions bundle terms that share a common coefficient into
braced groups. The choice is stored persistently across Julia sessions.
"""
function set_do_braced(mode::Bool)
    _update_pref!(:DO_BRACED, mode)
end

"""
    set_expand_cumulants(mode::Bool)

Choose whether cumulant expressions are emitted in fully expanded form (`true`) or kept
in the more compact symbolic representation (`false`). The chosen mode is saved using
`Preferences.jl`.
"""
function set_expand_cumulants(mode::Bool)
    _update_pref!(:EXPAND_CUMULANTS, mode)
end

"""
    set_quadgk_atol(value)

Persist the absolute tolerance used by default when evaluating quadrature integrals.
"""
function set_quadgk_atol(value)
    _update_pref!(:QUADGK_ATOL, Float64(value))
end

"""
    set_quadgk_rtol(value)

Persist the relative tolerance used by default when evaluating quadrature integrals.
"""
function set_quadgk_rtol(value)
    _update_pref!(:QUADGK_RTOL, Float64(value))
end

"""
    set_sample_atol(value)

Persist the absolute tolerance used by default when generating ensemble samples.
"""
function set_sample_atol(value)
    _update_pref!(:SAMPLE_ATOL, Float64(value))
end

"""
    set_sample_rtol(value)

Persist the relative tolerance used by default when generating ensemble samples.
"""
function set_sample_rtol(value)
    _update_pref!(:SAMPLE_RTOL, Float64(value))
end



export get_default, set_flip_if_first_term_negative, set_do_braced, set_expand_cumulants,
       set_quadgk_atol, set_quadgk_rtol, set_sample_atol, set_sample_rtol,
       FLIP_IF_FIRST_TERM_NEGATIVE , DO_BRACED, EXPAND_CUMULANTS,
       QUADGK_ATOL, QUADGK_RTOL, SAMPLE_ATOL, SAMPLE_RTOL, pushindex!

include("Helper.jl")
include("QSort.jl")
include("OffsetArrays.jl")
export OffsetArray

include("QIndexes.jl")
using .QIndexes
export QIndex, AbstractIndex, TimeIndex, SampleIndex

include("StringUtils.jl")
using .StringUtils
export symbol2formatted, str2sub, str2sup, brace, braket, indices2str
export int_exponent2str, exponentdag2str

include("QSpaceOps/Sampler.jl")
using .Sampler

export QDistribution, QNormal, QUniform, QEnsembleFunction, pdf
export QInterpolator, build_interpolation_nodes, eval_interpolation, nodes, basis_values, basis_values!
export QIntegrator, integrate_node_funs, normalization_constant, eval_integration
export build_discrete_samples, build_continuous_samples
export EnsembleSamples, AbstractEnsembleSample, DiscreteSamples, ContinuousSamples

include("ParameterGroups.jl")
using .ParameterGroups
export ParameterGroupKind, ParameterGroup, ParameterGroupScalar, ParameterGroupTimeScalar, ParameterGroupTimeFunction, ParameterGroupDistribution, ParameterGroupEnsembleFunction, ParameterGroupEnsembleTimeFunction

include("CFunctions.jl")
using .CFunctions
export CFunction, CAbstract, CIntegral, CCustomType, CCustomTypeIndexed, CParticle, CAtom, CSum, CRational, CProd, CExp, CLog
export CMatrix, CVector, CPower, CAtomIndexed, CCustomTypeIndexed
export define_cabstract, define_ctype, define_cintegral
export reorder, evaluate, stringer, to_stringer, to_string, sort_key
export coeff, var_exponents, expand, substitute
export contains_non_simple_CFunction, Indexed
export list_cabstracts, list_ctypes, list_cintegrals
export where_acting, where_acting!, which_ensemble_acting, which_ensemble_acting!
export which_params_acting, which_params_acting!, separate_by_cond
export ParameterValues, update_t!, resolve_param!, resolve_ensemble_values!, refresh_ensemble_values!
export WhereWhichParamGroup, AbstractIndexParameters

include("Cumulants.jl")
using .Cumulants

include("QSpace.jl")
using .QSpaces
export OperatorSet, max_operator_magnitude
export Ensemble, SubSpace, SubSpaceDefinitions, ConcreteIndexes, SubSpaceIndex 
export OperatorType, OperatorTypeInfo, OperatorDefinitions
export Parameter, ParameterInfo, ParameterDefinitions, set_parameter_group_definition!, SubSpaceIndex, EnsembleIndex
export QSpace, SubSpaceDicts, AbstractOperatorDicts
export get_parameter_group, get_parameter, get_subspace, get_subspace_index, get_ensemble, get_operator_type
export QubitPauli, QubitPM, Ladder
export DiscreteSamples, ContinuousSamples

include("Plotting.jl")

include("QExpressions.jl")
#include("CFunctions_updates.jl")
using .QExpressions
using .QExpressions: diffQEqOrdered
export QObj, QAtom, QComposite, QCompositeN, QMultiComposite, QAbstract, QTerm, QExpr, QCumulant, QAtomProduct, permutation, AbstractQSum, QSum, QIntegral, QInt, ∑, ∫, integral, QCompositeProduct, diffQEq, d_dt
export QCommutator, QExp, QLog, QPower, power, QRoot, root, Expectation
export Dag, Commutator, same_qspace
export is_t_var, is_local, contains_non_simple_QObj, contains_non_simple, contains_abstract, contains_time, max_order_of_terms, contains_which_t_indices, iter_QAtomProducts, iter_QInts
export is_unitary, is_hermitian, substitution_properties_fulfilled
export base_operators, get_parameter, get_operator, get_abstract, QExprLookup 
export string, latex_string
export term
export @define, @define_basics, QExpr2CFunction, Cumulant, cumulant_string

export contains_abstract, are_indices_defined, which_summations_acting, which_summations_to_root, are_all_neq
export Substitution, --> 
export reorder, reorder_full, reorder_time, neq, flatsums, complexsums , Sum2Int
# Preindexing
export QAtomOrdered, QAtomIndexed, QNeutral, OrderbyOperator, OrderedQAtomProduct, OrderedQExpr, OrderedDiffQEq, diffQEqOrdered

include("ConcreteIndexingHelper/Indexing.jl")
using .Indexing
export BinomialCache, _UsedBuf, EnsembleRankWorkspace, MultiEnsembleWorkspace
export index_rank_for_ensemble!, index_rank_for_ensemble_continuum!
export index_number_for_ensemble, index_number_for_ensemble_continuum
export combined_index_for_ensembles!, combined_index_for_ensembles_continuum!

include("ConcreteIndexingHelper/Combinatorics.jl")
using .QCombinatorics
export ensemble_iterator, ensemble_iterator_continuum, multi_ensemble_iterator, multi_ensemble_iterator_continuum

#include("QEqSets.jl")
#using .QEqSets
export diffQEqOrdered
end # module QAlgebra
