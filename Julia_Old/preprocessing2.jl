using Combinatorics
using Base.Threads
include("operator_terms.jl")
include("diff_Eq.jl")
include("indexing.jl")
include("cumulants.jl")
include("sampled_integrals.jl")
include("preprocessing.jl")
include("print_terms.jl")

mutable struct SystemAndDicts
    all_eqs_indexed::Vector{DE_Term_indexed}
    operator_strings_dict::Dict{String,Op_Group_Type}
    how_many_total::Int
    all_combinations_list::Vector{Vector{Vector{Int}}}
    vars_vec::Vector{Vector{String}}
    cumulant_terms::Vector{Cumulant_indexed}
    cumulant_terms_dict::Dict{String,Op_Group_Type}
    lower_order_cumulant_terms::Vector{Full_Cumulant_indexed}
    do_lower_order_cumulants::Bool
    order::Int
    spin_order::Int
    is_pauli::Bool
    reduced_system::Bool

    function SystemAndDicts(all_eqs_indexed::Vector{DE_Term_indexed}, operator_strings_dict::Dict{String,Op_Group_Type}, how_many_total::Int, all_combinations_list::Vector{Vector{Vector{Int}}}, vars_vec::Vector{Vector{String}}, cumulant_terms::Vector{Cumulant_indexed}, cumulant_terms_dict::Dict{String,Op_Group_Type}, order::Int, spin_order::Int, is_pauli::Bool, lower_order_cumulant_terms::Union{Vector{Full_Cumulant_indexed},Bool}=false, reduced_system::Bool=false)
        do_lower_order_cumulants::Bool = true
        if typeof(lower_order_cumulant_terms) == Bool
            lower_order_cumulant_terms = Vector{Full_Cumulant_indexed}()
            do_lower_order_cumulants = false
        end
        new(all_eqs_indexed, operator_strings_dict, how_many_total, all_combinations_list, vars_vec, cumulant_terms, cumulant_terms_dict, lower_order_cumulant_terms, do_lower_order_cumulants, order, spin_order, is_pauli, reduced_system)
    end
end

# Reduce size of terms variable, by removing parts that are no longer needed after a calculation
function reduced_system(sys::SystemAndDicts)
    # empty the all_eqs_indexed, cumulant_terms, lower_order_cumulant_terms vectors
    all_eqs_indexed::Vector{DE_Term_indexed} = []
    cumulant_terms::Vector{Cumulant_indexed} = []
    lower_order_cumulant_terms::Vector{Full_Cumulant_indexed} = []
    all_combinations_list::Vector{Vector{Vector{Int}}} = []
    return SystemAndDicts(all_eqs_indexed, sys.operator_strings_dict, sys.how_many_total, all_combinations_list, sys.vars_vec, cumulant_terms, sys.cumulant_terms_dict, sys.order, sys.spin_order, sys.is_pauli, lower_order_cumulant_terms, true)
end
function reduced_system!(sys::SystemAndDicts)
    # empty the all_eqs_indexed, cumulant_terms, lower_order_cumulant_terms vectors
    sys.all_eqs_indexed = []
    sys.cumulant_terms = []
    sys.lower_order_cumulant_terms = []
    sys.all_combinations_list = []
    return sys
end
## Test 
#sys = reduced_system(system)

function index_to_term(de_term_indexed::DE_Term_indexed, system::SystemAndDicts, samples::SamplesAndWeights)::DE_Term
    # Transform a vector of indexed terms into a vector of terms
    operator_strings_dict = system.operator_strings_dict
    cumulant_terms_dict = system.cumulant_terms_dict
    vars_vec = system.vars_vec
    sample_num = samples.sample_num
    return index_to_term(de_term_indexed, operator_strings_dict, cumulant_terms_dict, vars_vec, sample_num)
end

function prepare_indexed_eqs_from_samples(max_order::Int, samples::SamplesAndWeights; do_lower_order_cumulants::Bool=true, less_spins::Int=0, max_spins::Int=-1, only_prepare::Bool=false, only_output_numbers::Bool=false, printing::Bool=true, do_progreqspace::Bool=true, threaded::Bool=true, is_pauli::Bool=true)::Union{SystemAndDicts,Vector{DE_Term_Multi}}
    # Extract sample information 
    if max_spins < 0
        if less_spins >= max_order
            error("less_spins must be less than max_order")
        end
        max_spins::Int = max_order - less_spins
    else
        less_spins = max_order - max_spins
    end
    sample_num::Int = samples.sample_num
    sample_vars::Dict{String,Vector{ComplexF64}} = samples.sample_vars
    weights_dict::Dict{String,Vector{Float64}} = samples.weights_dict
    locations = samples.locations
    if !only_prepare && printing
        println("---- Generating Indexed Differential Equations of Order ", max_order, " with max. spin order ", max_order - less_spins, " and ", samples.sample_num, " samples ----")
    end
    all_eqs = multi_spin_DE_operators_up_to_order(max_order, less_spins, is_pauli=is_pauli)

    conjugate_to_basis!(all_eqs)
    t0 = time()
    vars_dict, how_many_cumulant, big_tuple = prepare_parse_DE_to_indexes(all_eqs, weights_dict, max_order, max_spins)
    how_many_total = big_tuple[1]
    num_eqs_str, num_cumulants_str = dotted_str(how_many_total), dotted_str(how_many_cumulant)
    num_eqs_str = lpad(num_eqs_str, length(num_cumulants_str))
    if printing
        println("   # of Equations: ", num_eqs_str)
        println("   # of Cumulants: ", num_cumulants_str)
    end
    if only_prepare
        return all_eqs
    end
    if only_output_numbers
        return how_many_total, how_many_cumulant
    end
    if printing
        println("Abstract Term Generation completed (", time2str(time() - t0), ")")
    end
    t0 = time()
    all_eqs_indexed, operator_strings_dict, cumulant_terms_dict = parse_DE_to_indexes(all_eqs, weights_dict, sample_vars, big_tuple, threaded=threaded)
    if printing
        println("Term Generation completed (", time2str(time() - t0), ")")
    end
    t0 = time()
    vars_vec = preprocess_index_inversion(vars_dict)
    multi_determine_indexes::Function = big_tuple[3]
    all_combinations_list::Vector{Vector{Vector{Int}}} = big_tuple[11]
    cumulant_terms = cumulant_terms_from_dict(cumulant_terms_dict, operator_strings_dict, multi_determine_indexes, all_combinations_list, printing=do_progress, max_spins, threaded=threaded)
    if printing
        println("Cumulants constructed (", time2str(time() - t0), ")")
        t0 = time()
    end
    if do_lower_order_cumulants
        lower_order_cumulant_terms = lower_order_full_cumulant_terms(operator_strings_dict, max_order, multi_determine_indexes, all_combinations_list, max_spins, threaded=threaded)
        if printing
            println("Lower Order Full Cumulants constructed (", time2str(time() - t0), ")")
            println("---- Done -------------------------------------------------------------------------------------------")
        end
        system = SystemAndDicts(all_eqs_indexed, operator_strings_dict, how_many_total, all_combinations_list, vars_vec, cumulant_terms, cumulant_terms_dict, max_order, max_spins, is_pauli, lower_order_cumulant_terms)
    else
        if printing
            println("---- Done -------------------------------------------------------------------------------------------")
        end
        system = SystemAndDicts(all_eqs_indexed, operator_strings_dict, how_many_total, all_combinations_list, vars_vec, cumulant_terms, cumulant_terms_dict, max_order, max_spins, is_pauli)
    end
    return system
end

