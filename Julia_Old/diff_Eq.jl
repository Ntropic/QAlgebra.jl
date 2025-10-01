include("operator_terms.jl")
include("print_terms.jl")
##############################################################################################################################
#### Single Qubit ############################################################################################################
##############################################################################################################################

function Trace_of_D_Operator(Op::Union{Term,Vector{Term}}, exp_op::Union{Term,Vector{Term}}, vars::Union{String,Vector{String}}=String[])::Union{Term,Vector{Term}}
    # returns the dissipator operator D[Op] \rho = Op \rho Op^\dagger - 1/2 (Op^\dagger Op \rho + \rho Op^\dagger Op) times the expected operator
    # Tr[exp_op D[Op] \rho] --> use trace property (Tr[AB] = Tr[BA]) to get \rho on the right side
    # Calculate Op^dagger * exp_op * Op - 1/2 (exp_op * Op^\dagger * Op + Op^\dagger * Op * exp_op)
    # Generate the 3 terms and sum them up
    Op_dag = dagger_term(Op)
    Op_dag_X_Op = multiply_terms(Op_dag, Op)

    # 1st Op^dagger * exp_op * Op
    Op_dag_X_exp_op = multiply_terms(Op_dag, exp_op)
    first_term = multiply_terms(Op_dag_X_exp_op, Op)
    # 2nd exp_op * Op^\dagger * Op    and    3rd Op^\dagger * Op * exp_op
    second_term = add_terms(multiply_terms(exp_op, Op_dag_X_Op), multiply_terms(Op_dag_X_Op, exp_op))
    scale_term!(second_term, -0.5)
    # combine terms
    terms = add_terms(first_term, second_term)
    if length(vars) > 0
        terms = add_vars(terms, vars)
    end
    return terms
end
## Test
#Op = make_term("z")
#exp_op = make_term("-x")
#display(Trace_of_D_Operator(Op, exp_op, "test_var")) # should be -2*"-x"

function Trace_of_Hamiltonian_commutator(exp_op::Union{Term,Vector{Term}}, index::String="")::Union{Term,Vector{Term}}
    # consists of delta_s term, g term and drive term (\sqrt{kappa}, beta)
    # Use Tr(exp_op * [H, rho]) = Tr([exp_op, H] * rho) = Tr([exp_op, H0] * rho) + Tr([exp_op, H1] * rho) + ...
    # for the parts H = H0 + H1 + ...
    is_pauli = exp_op.is_pauli
    index_str::String = ""
    if length(index) > 0
        index_str = "_" * index
    end
    z = make_term("z" * index, 1 / 2, is_pauli=is_pauli)
    delta_s_terms = add_vars(commutator_terms(exp_op, z, -1.0im), raw"\Delta" * index_str)
    a_sigma_plus = multiply_terms(make_term("-", is_pauli=is_pauli), sigma_plus(index, is_pauli=is_pauli))  # a_- * sig_+
    a_dag_sigma_minus = multiply_terms(make_term("+", is_pauli=is_pauli), sigma_minus(index, is_pauli=is_pauli))  # a_+ * sig_-
    g_terms = add_vars(commutator_terms(exp_op, add_terms(a_sigma_plus, a_dag_sigma_minus), -1.0im), "g" * index_str)
    brace_term = [make_term("+", 1.0im, [raw"\sqrt{\kappa}", raw"\beta"], is_pauli=is_pauli), make_term("-", -1.0im, [raw"\sqrt{\kappa}", raw"\beta^{*}"], is_pauli=is_pauli)]
    drive_terms = commutator_terms(exp_op, brace_term, -1.0im)
    terms = add_terms(delta_s_terms, g_terms, drive_terms)
end
## Test
#exp_op = make_term("-")
#display(Trace_of_Hamiltonian_commutator(exp_op)) 

function time_derivative_expectation_of_operator(exp_op::Term, index::String="")::DE_Term
    # Liouvillian = -i [H, \rho] \rho + 1/2*\gamma D[\sigma_z] \rho + \Gamma D[\sigma_-] + \kappa D[a] \rho
    is_pauli = exp_op.is_pauli
    which_ind = get_indexes(exp_op)
    de_term::DE_Term = DE_Term(exp_op, which_ind[1])
    term_ham = Trace_of_Hamiltonian_commutator(exp_op, index)
    term_gamma = scale_term(Trace_of_D_Operator(make_term("z" * index, is_pauli=is_pauli), exp_op, raw"\gamma"), 0.5)
    term_Gamma = Trace_of_D_Operator(sigma_minus(index, is_pauli=is_pauli), exp_op, raw"\Gamma")
    term_kappa = Trace_of_D_Operator(make_term("-", is_pauli=is_pauli), exp_op, raw"\kappa")
    de_term.terms = add_terms(term_ham, term_gamma, term_Gamma, term_kappa)
    return de_term
end
## Test
#exp_op = make_term("-")
#display(time_derivative_expectation_of_operator(exp_op))


struct DE_Term_indexed # modified DE_Term, but with indexes instead of strings
    exp_index::Int                                            # operator index
    linear_terms::Vector{Indexed_Term}
    cumulant_terms::Vector{Indexed_Term}
    constant_terms::Vector{Constant_Term}
    #linear_terms_summed::Vector{Indexed_Term}
    #cumulant_term_summed::Vector{Indexed_Term}
    clamped::Bool   # is this term clamped from -1.0 to 1.0?
end

function Base.show(io::IO, ::MIME"text/plain", x::DE_Term_indexed)
    println(io, "DE_Term_indexed")
    println(io, " - exp_index: " * string(x.exp_index))
    # print vector elements of linear_terms as table with columns (operator index, conjugate, variables index, coefficient)
    #linear_terms = ["Op. Ind." "Conj?" "Var. Ind." "Coeff."; [string(i) string(c) string(v) string(coeff) for (i, c, v, coeff) in x.linear_terms]]
    println(io, " - linear_terms: ")
    print(io, pretty_table_from_indexed_terms(x.linear_terms, ["i(Op.)", "Conj.", "i(Var.)", "Coeff."]))
    println(io, " - cumulant_terms: ")
    print(io, pretty_table_from_indexed_terms(x.cumulant_terms, ["i(Cum.)", "Conj.", "i(Var.)", "Coeff."]))
    println(io, " - constant_terms: ")
    print(io, pretty_table_from_constant_terms(x.constant_terms, ["i(Var.).", "Coeff."]))
end



function get_indexes(exp_op::Union{Term,Vector{Term}})::Tuple{Vector{String},Bool,Bool}
    # Get the indexes of the operators and return an array of the indexes
    spin_indices = exp_op.spin_indices
    # get only the unique and sorted indexes into a single joined string
    indexes = string.(sort(unique(spin_indices)))
    spin_op::Bool = length(exp_op.spin_indices) > 0 ? true : false
    cavity_op::Bool = length(exp_op.bosons) > 0 ? true : false
    return indexes, spin_op, cavity_op
end
## test 
#exp_op = make_term("zi*xj*yh+-")
#display(get_indexes(exp_op))

function get_delta_and_g_str(index::String="", index_delta::Bool=true, index_g::Bool=True)::Tuple{String,String}
    # Get the delta and g strings
    index_str::String = ""
    if length(index) > 0
        index_str = "_" * index
    end
    delta_str = index_delta ? raw"\Delta" * index_str : raw"\Delta"    # add index to delta if index_Delta is true
    g_str = index_g ? "g" * index_str : "g"    # add index to g if index_g is true
    return delta_str, g_str
end
## test
#display(get_delta_and_g_str("i", true, true))

function Trace_of_L_i_operator(exp_op::Term, index::String="", index_delta::Bool=true, index_g::Bool=true)::Vector{Term}
    is_pauli = exp_op.is_pauli
    delta_str, g_str = get_delta_and_g_str(index, index_delta, index_g)
    z = make_term("z" * index, 1 / 2, is_pauli=is_pauli)
    delta_s_terms = add_vars(commutator_terms(exp_op, z, -1.0im), delta_str)
    a_sigma_plus = multiply_terms(make_term("-", is_pauli=is_pauli), sigma_plus(index, is_pauli=is_pauli))  # a_- * sig_+
    a_dag_sigma_minus = multiply_terms(make_term("+", is_pauli=is_pauli), sigma_minus(index, is_pauli=is_pauli))  # a_+ * sig_-
    g_terms = add_vars(commutator_terms(exp_op, add_terms(a_sigma_plus, a_dag_sigma_minus), -1.0im), g_str)
    term_gamma = scale_term(Trace_of_D_Operator(make_term("z" * index, is_pauli=is_pauli), exp_op, raw"\gamma"), 0.5)
    term_Gamma = Trace_of_D_Operator(sigma_minus(index, is_pauli=is_pauli), exp_op, raw"\Gamma")
    return add_terms(delta_s_terms, g_terms, term_gamma, term_Gamma)
end
## Test
#exp_op = make_term("-xi*yj")
#display(Trace_of_L_i_operator(exp_op, "h"))

# Complete Liouvillian
function time_derivative_expectation_of_operator_multispin(exp_op::Term, index_delta::Bool=true, index_g::Bool=true)::DE_Term_Multi
    # Liouvillian = -i [H, \rho] \rho + 1/2*\gamma D[\sigma_z] \rho + \Gamma D[\sigma_-] + \kappa D[a] \rho# get all indexes of exp_op
    # if i is in indexes change
    # in multi spin case, we construct the hamiltonian in decomposed way (see paper of this project)
    # so that the list all_terms has the structure:
    #               [ sum_i , sum_i \neq ind, i=ind[0] + ... + i=ind[-1] + term_kappa+drive_terms ]
    # discuss the different cases:
    # case 1 s + c: [       , sum_i \neq ind, i=ind[0] + ... + i=ind[-1] + term_kappa+drive_terms ]
    # case 2 s + _: [       ,               , i=ind[0] + ... + i=ind[-1]                          ]
    # case 3 _ + c: [ sum_i ,               ,                            + term_kappa+drive_terms ]
    which_ind, s, c = get_indexes(exp_op)  # s=has spin, c=has cavity
    is_pauli = exp_op.is_pauli
    which_case::Int = 4
    if s
        which_case -= 2
    end
    if c
        which_case -= 1
    end
    # Determine sum index
    sum_ind::String = "i"
    i = "i"
    if i in which_ind
        for i in ["h", "g", "f", "e", "d", "c", "b", "a"]
            if !(i in which_ind)
                sum_ind = i
                break
            end
        end
    end
    de_term::DE_Term_Multi = DE_Term_Multi(exp_op, which_ind, sum_ind)

    delta_str, g_str = get_delta_and_g_str(sum_ind, index_delta, index_g)

    if which_case == 3          #### do_g_term (sum_i)
        a_sigma_plus = multiply_terms(make_term("-", is_pauli=is_pauli), sigma_plus(sum_ind, is_pauli=is_pauli))  # a_- * sig_+
        a_dag_sigma_minus = multiply_terms(make_term("+", is_pauli=is_pauli), sigma_minus(i, is_pauli=is_pauli))  # a_+ * sig_-
        sum_i = add_vars(commutator_terms(exp_op, add_terms(a_sigma_plus, a_dag_sigma_minus), -1.0im), g_str)
        de_term.sum_i_terms = combine_terms(sum_i)
    elseif which_case == 1      ####  do sum_i \neq ind term(s)
        de_term.sum_i_neq_j_terms = combine_terms(Trace_of_L_i_operator(exp_op, sum_ind, index_delta, index_g))
    end

    curr_terms::Vector{Term} = Term[]
    if which_case in [1, 2]  #### do terms with: i=ind[0] + ... + i=ind[-1] 
        for i in which_ind
            append!(curr_terms, Trace_of_L_i_operator(exp_op, i))
        end
    end
    # add spin independent terms
    brace_term = add_terms(make_term("+", 1.0im, [raw"\sqrt{\kappa}", raw"\beta"], is_pauli=is_pauli), make_term("-", -1.0im, [raw"\sqrt{\kappa}", raw"\beta^{*}"], is_pauli=is_pauli))
    drive_terms = commutator_terms(exp_op, brace_term, -1.0im)
    term_kappa = Trace_of_D_Operator(make_term("-", is_pauli=is_pauli), exp_op, raw"\kappa")
    append!(curr_terms, add_terms(term_kappa, drive_terms))
    de_term.non_sum_terms = combine_terms(curr_terms)
    return de_term
end
## Test 
#exp_op = make_term("-")
#display(time_derivative_expectation_of_operator_multispin(exp_op))

###### To Matrix
# In case it wasn't separated into different orders
function make_terms_matrix(term_vec::Union{Vector{Term},Vector{DE_Term},Vector{DE_Term_Multi}})::Union{Vector{Vector{Term}},Vector{Vector{DE_Term}},Vector{Vector{DE_Term_Multi}}}
    # Make a matrix of terms from a vector of terms
    type = typeof(term_vec[1])
    vecvec::Vector{Vector{type}} = Vector{Vector{type}}()
    curr_order::Int = 0
    curr_max_term::Int = 0
    for term in term_vec
        if typeof(term) == Term
            curr_order, _, _ = how_many_operators_in_term(term)
        else
            curr_order, _, _ = how_many_operators_in_term(term.exp_op)
        end
        if curr_order > curr_max_term
            # generate subvectors for all orders up to curr_order
            for i in curr_max_term+1:curr_order
                push!(vecvec, Vector{type}())
            end
            curr_max_term = curr_order
        end
        push!(vecvec[curr_order], term)
    end
    return vecvec
end
#all_eqs_flat = Vector{DE_Term}()
#for eqs in all_eqs
#    for eq in eqs
#        push!(all_eqs_flat, eq)
#    end
#end
#all_eq_vec = make_terms_matrix(all_eqs_flat)
