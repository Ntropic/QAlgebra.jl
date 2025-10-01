
include("diff_Eq_solver.jl")
include("combinatorics.jl")


struct DensityMatrix  # Simply a Density Matrix, no target states or anything else!
    sample_indexes::Vector{Int}
    operator_indexes::Vector{Int}
    cumulant_indexes::Vector{Int}
end

struct MatrixGenerator  # for every term in the densitymatrices, the corresponding matrix is stored and the operator names are stored
    constant_matrix::Matrix{ComplexF64}   # 
    operator_matrices::Vector{Matrix{ComplexF64}}
    operator_minus::Vector{Bool}
    operator_num::Vector{Bool}
    operator_names::Vector{String}
    cumulant_matrices::Vector{Matrix{ComplexF64}}
    cumulant_minus::Vector{Bool}  # minus terms have only a - for the cavity, <->+<+> => x and  i<-> - i<+> => y
    cumulant_num::Vector{Bool}    # num terms (+-) are limited to abs = 1
    cumulant_names::Vector{String}
    dim::Int
end

function Fidelity_from_expectation_values(sigma_mat::Matrix{T}, sqrt_rho::Matrix{T}) where {T<:ComplexF64}
    # Calculate abs(Tr(sqrt(sqrt(rho)*sigma*sqrt(rho))))^2

    F = abs(LinearAlgebra.tr(sqrt(sqrt_rho * sigma_mat * sqrt_rho)))^2
    if F > 1.0 # get rid of a little bit of noise, 
        F = 1.0
    end
    return F
end
function Infidelity_from_expectation_values(sigma_mat::Matrix{T}, sqrt_rho::Matrix{T}) where {T<:ComplexF64}
    return 1 - Fidelity_from_expectation_values(sigma_mat, sqrt_rho,)
end

function num_i_geq_j(nums::Vector{Int}, i::Int, j::Int)
    # check if number i appears more often (or equally many times) in nums, as number j
    return sum(nums .== i) >= sum(nums .== j)  
end
function conditional_num_i_geq_j(nums::Vector{Int}, i::Int, j::Int, do_check::Bool=false)
    if do_check
        return num_i_geq_j(nums, i, j)
    else
        return true
    end
end
function spin_combinations_separated_indexes(max_spins::Int, is_pauli::Bool=false)::Vector{Tuple{String,Vector{Int}, Bool}}
    # all combinations of 'x,y,z' for a given number of spins with (greater equal ordering)
    # returns a Vector of tuples, where the first element is the string of the spin operators and the second element is the corresponding index (which spins are the operators referring to) and finally a Bool to indicate whether we need to consider the conjugate of the operator
    spin_types::Vector{String} = []
    if is_pauli
        spin_types = ["i", "x", "y", "z"]
    else
        spin_types = ["i", "p", "m", "z"]
    end
    spin_combinations::Vector{Tuple{String,Vector{Int}, Bool}} = []

    # all combinations of number_of_spins [x,y,z]
    for i_s in mrange([length(spin_types) for i in 1:max_spins])
        # how_many 1's in i_s 
        do_conjugate::Bool = false
        if !conditional_num_i_geq_j(i_s, 2, 3, !is_pauli)
            do_conjugate = true
            # flip the p and m terms 
            for i in 1:length(i_s)
                if i_s[i] == 2
                    i_s[i] = 3
                elseif i_s[i] == 3
                    i_s[i] = 2
                end
            end
        end
        n_is::Int = count(i -> i == 1, i_s)
        # Sort the string x before y before z and generate an indexing off the sorting 
        order::Vector{Int} = sortperm(i_s)[1+n_is:end]
        #println(order, " ", i_s)
        curr_comb::String = ""
        for i in 1:length(order)
            curr_comb *= spin_types[i_s[order[i]]]
        end
        push!(spin_combinations, (curr_comb, order, do_conjugate))
    end
    return spin_combinations
end
function tuple_vec_to_dict(tuple_vec::Vector{Tuple{String,Int}})::Dict{String,Vector{Int}}
    dict::Dict{String,Vector{Int}} = Dict()
    for tup in tuple_vec
        if haskey(dict, tup[1])
            push!(dict[tup[1]], tup[2])
        else
            dict[tup[1]] = [tup[2]]
        end
    end
    return dict
end
function all_indexed_spin_operators(max_number_of_spins::Int, is_pauli::Bool=false)::Tuple{Vector{Vector{Int}},Dict{String,Vector{Int}},Vector{String}}
    # Returns the different possible index combinations for the operators, 
    # then a dictionary, with the index combinations as values and the operator as key (negative signs indicate conjugates of operators), 
    # finally a list of the operators
    index_combinations::Vector{Vector{Int}} = Vector{Vector{Int}}()
    new_combinations::Vector{Tuple{String,Int}} = []

    my_combinations::Vector{Tuple{String,Vector{Int}, Bool}} = spin_combinations_separated_indexes(max_number_of_spins, is_pauli)
    index_list = [c[2] for c in my_combinations]
    index_combinations = unique(index_list)
    # for every combination, find the index in index_combinations and append it to the new_combinations
    for comb in my_combinations
        index = findfirst(x -> x == comb[2], index_combinations)
        if comb[3] == true 
            index = -index  # sign specifies 
        end
        push!(new_combinations, (comb[1], index))
    end
    spin_index_combinations = tuple_vec_to_dict(new_combinations)
    string_combinations::Vector{String} = unique([s[1] for s in new_combinations])
    #string_combinations::Vector{String} = [k for k in keys(spin_index_combinations)]
    return index_combinations, spin_index_combinations, string_combinations
end
## Test
#all_indexed_spin_operators(2, false)
function all_indexed_spin_and_cavity_operators(order::Int, is_pauli::Bool; do_minus::Bool=false)
    # create also all ordered combinations of + and - for the cavity (with number of + <= number of -)
    # Returns the different possible index combinations for the operators, 
    # then a dictionary, with the index combinations as values and the operator as key (negative signs indicate conjugates of operators), 
    # finally a list of the operators
    less_spins = 1
    spin_order = order - less_spins
    index_combinations, spin_index_combinations, spins_list = all_indexed_spin_operators(spin_order, is_pauli)
    # for every spin index combination, add all combinations of -, +-, -- etc. for the cavity
    cavity_spin_index_combinations = Dict{String,Vector{Int}}()
    string_combinations::Vector{String} = []
    if do_minus
        cavity_combinations = ["", "-", "+-"]
    else
        cavity_combinations = ["", "+-"]
    end
    for c in cavity_combinations
        for spin_comb in spins_list
            spin_index_pre = spin_index_combinations[spin_comb]
            #spin_index = abs.(spin_index_pre)
            #conjugates = spin_index_pre .< 0
            curr_spin_order = length(spin_comb)
            new_string = c * spin_comb
            cavity_spin_index_combinations[new_string] = spin_index_pre
            push!(string_combinations, new_string)
        end
    end
    return index_combinations, cavity_spin_index_combinations, string_combinations
end
## Test 
#index_combinations, cavity_spin_index_combinations, string_combinations = all_indexed_spin_and_cavity_operators(2)


####### Prepare ############################################################################################################
# Take system variable, sample variable, order  
function make_mat(curr_combinations::Vector{Vector{Int}}, curr_string::String, spin_order::Int)::Tuple{Vector{Matrix{ComplexF64}},Vector{Bool},Vector{Bool}}
    x::Matrix{ComplexF64} = Matrix{ComplexF64}([0.0 1.0; 1.0 0.0])
    y::Matrix{ComplexF64} = Matrix{ComplexF64}([0.0 -im; im 0.0])
    z::Matrix{ComplexF64} = Matrix{ComplexF64}([1.0 0.0; 0.0 -1.0])
    I::Matrix{ComplexF64} = Matrix{ComplexF64}([1.0 0.0; 0.0 1.0])
    minus::Matrix{ComplexF64} = Matrix{ComplexF64}([0.0 1.0; 0.0 0.0])
    plus::Matrix{ComplexF64} = Matrix{ComplexF64}([0.0 0.0; 1.0 0.0])
    n::Matrix{ComplexF64} = -2 * z
    Ipz::Matrix{ComplexF64} = I + z
    spin_matrix_list::Vector{Matrix{ComplexF64}} = [I, x, y, z]
    cavity_matrix_list::Vector{Matrix{ComplexF64}} = [Ipz, minus, n]
    int_v::Int = Int('v')

    curr_matrices::Vector{Matrix{ComplexF64}} = []
    operator_minus::Vector{Bool} = []
    operator_num::Vector{Bool} = []
    curr_cavity_index::Int = 1
    n_m::Int = count(==('-'), curr_string)
    n_p::Int = count(==('+'), curr_string)
    curr_minus = false
    curr_num = false
    if n_m > 0
        if n_p > 0
            curr_cavity_index = 3
            curr_num = true
        else
            curr_cavity_index = 2
            curr_minus = true
        end
    end
    factor = 1 / 2^(spin_order + 1)
    for curr_comb in curr_combinations
        curr_mat_indexes::Vector{Int} = ones(Int, spin_order)
        for i in 1:length(curr_comb)
            c = curr_comb[i]
            s = curr_string[n_p+n_m+i]
            # 2 if x, 3 if y, 4 if z
            curr_index = Int(s) - int_v
            curr_mat_indexes[c] = curr_index
        end
        curr_mat_vec::Vector{Matrix{ComplexF64}} = [spin_matrix_list[i] for i in curr_mat_indexes]
        curr_mat::Matrix{ComplexF64} = kron(cavity_matrix_list[curr_cavity_index], curr_mat_vec...)

        push!(curr_matrices, curr_mat * factor)
        push!(operator_minus, curr_minus)
        push!(operator_num, curr_num)
    end
    return curr_matrices, operator_minus, operator_num
end
function Density_Matrix_Generator(system::SystemAndDicts, samples::SamplesAndWeights; threaded::Bool=true, do_minus::Bool=false)::Tuple{Vector{DensityMatrix},MatrixGenerator}
    # Generate density matrix structs for different spin combinations from sampled spins :
    order::Int = system.order
    spin_order::Int = system.spin_order
    if order == spin_order
        spin_order = order - 1
    end
    is_pauli::Bool = system.is_pauli
    sample_num::Int = samples.sample_num
    index_combinations, cavity_spin_index_combinations, string_combinations = all_indexed_spin_and_cavity_operators(order, is_pauli, do_minus=do_minus)
    multi_determine_indexes = determine_multi_indexes_gen(sample_num, spin_order)[1]
    # separate string combinations into cumulant and linear terms 
    linear_combinations::Vector{String} = []
    cumulant_combinations::Vector{String} = []
    linear_first_index::Vector{Int} = []
    cumulant_first_index::Vector{Int} = []
    linear_op_spin_orders::Vector{Vector{Int}} = []
    cumulant_op_spin_orders::Vector{Vector{Int}} = []
    linear_index_combinations::Vector{Vector{Vector{Int}}} = []
    cumulant_index_combinations::Vector{Vector{Vector{Int}}} = []
    linear_how_many::Vector{Int} = []
    cumulant_how_many::Vector{Int} = []
    for s in string_combinations # check if key is in system.operator_strings_dict or system.cumulant_terms_dict
        curr_combinations_indexes = cavity_spin_index_combinations[s]
        curr_combinations = [index_combinations[i] for i in curr_combinations_indexes]
        if haskey(system.operator_strings_dict, s)
            push!(linear_combinations, s)
            op_spin_orders = system.operator_strings_dict[s].op_spin_orders
            first_index = system.operator_strings_dict[s].first_index
            how_many = system.operator_strings_dict[s].how_many
            push!(linear_first_index, first_index)
            push!(linear_op_spin_orders, op_spin_orders)
            push!(linear_index_combinations, curr_combinations)
            push!(linear_how_many, how_many)
        elseif haskey(system.cumulant_terms_dict, s)
            push!(cumulant_combinations, s)
            op_spin_orders = system.cumulant_terms_dict[s].op_spin_orders
            first_index = system.cumulant_terms_dict[s].first_index
            how_many = system.cumulant_terms_dict[s].how_many
            push!(cumulant_first_index, first_index)
            push!(cumulant_op_spin_orders, op_spin_orders)
            push!(cumulant_index_combinations, curr_combinations)
            push!(cumulant_how_many, how_many)
        elseif length(s) > 0
            error("Error: $(s) is not in the system (checked in operator_strings_dict and cumulant_terms_dict)")
        end
    end
    how_many_linear::Int = sum([length(i) for i in linear_index_combinations])
    how_many_cumulant::Int = sum([length(i) for i in cumulant_index_combinations])

    # Determine the number of terms via mgreater_range
    sample_index_combinations::Vector{Vector{Int}} = collect(mgreater_range(sample_num, spin_order))
    number_of_density_matrices::Int = length(sample_index_combinations)
    density_matrices::Vector{DensityMatrix} = Vector{DensityMatrix}(undef, number_of_density_matrices)
    @usethreads threaded for ind in 1:number_of_density_matrices
        #   for ind in 1:number_of_density_matrices  
        i::Vector{Int} = sample_index_combinations[ind]
        linear_terms::Vector{Int} = Vector{Int}(undef, how_many_linear)
        cumulant_terms::Vector{Int} = Vector{Int}(undef, how_many_cumulant)
        # iterate over all linear strings and cumulant strings and generate the density matrices
        counter::Int = 1
        for (curr_first_index, curr_spin_order, curr_combinations, curr_how_many) in zip(linear_first_index, linear_op_spin_orders, linear_index_combinations, linear_how_many)
            # generate the density matrix for the current sample index combination and the current linear string
            for comb in curr_combinations
                curr_comb = i[comb]
                inserted_indexes = indexes_to_index_vector(curr_comb, curr_spin_order)
                curr_index = multi_determine_indexes(inserted_indexes) + curr_first_index - 1
                if curr_index - curr_first_index > curr_how_many
                    error("Linear Error: $(curr_index - curr_first_index) > $(curr_how_many)   -   $(curr_spin_order)  -  $(curr_comb)")
                end
                linear_terms[counter] = curr_index
                counter += 1
            end
        end
        counter = 1
        for (curr_first_index, curr_spin_order, curr_combinations, curr_how_many) in zip(cumulant_first_index, cumulant_op_spin_orders, cumulant_index_combinations, cumulant_how_many)
            # generate the density matrix for the current sample index combination and the current cumulant string
            for comb in curr_combinations
                curr_comb = i[comb]
                inserted_indexes = indexes_to_index_vector(curr_comb, curr_spin_order)
                curr_index = multi_determine_indexes(inserted_indexes) + curr_first_index - 1
                if curr_index - curr_first_index > curr_how_many
                    error("Cumulant Error: $(curr_index - curr_first_index) > $(curr_how_many)   -   $(curr_spin_order)  -  $(curr_comb)")
                end
                cumulant_terms[counter] = curr_index
                counter += 1
            end
        end
        # store the density matrix in the density matrix struct
        density_matrix = DensityMatrix(i, linear_terms, cumulant_terms)
        # store the density matrix in the density matrix struct
        density_matrices[ind] = density_matrix
    end

    # generate_matrix
    linear_matrices::Vector{Matrix{ComplexF64}} = Vector{Matrix{ComplexF64}}(undef, how_many_linear)
    cumulant_matrices::Vector{Matrix{ComplexF64}} = Vector{Matrix{ComplexF64}}(undef, how_many_cumulant)
    operator_minus::Vector{Bool} = Vector{Bool}(undef, how_many_linear)
    cumulant_minus::Vector{Bool} = Vector{Bool}(undef, how_many_cumulant)
    operator_names::Vector{String} = Vector{String}(undef, how_many_linear)
    cumulant_names::Vector{String} = Vector{String}(undef, how_many_cumulant)
    operator_num::Vector{Bool} = Vector{Bool}(undef, how_many_linear)
    cumulant_num::Vector{Bool} = Vector{Bool}(undef, how_many_cumulant)
    # for every index combination generate the corresponding matrix using I (identity) matrices for the missing indexes and the x,y,z matrices for the present indexes using kronecker products
    counter2 = 1
    for (curr_combinations, curr_string) in zip(linear_index_combinations, linear_combinations)
        curr_matrices, curr_minuses, curr_nums = make_mat(curr_combinations, curr_string, spin_order)
        linear_matrices[counter2:counter2+length(curr_matrices)-1] = curr_matrices
        operator_minus[counter2:counter2+length(curr_matrices)-1] = curr_minuses
        operator_num[counter2:counter2+length(curr_matrices)-1] = curr_nums
        operator_names[counter2:counter2+length(curr_matrices)-1] = [curr_string for i in 1:length(curr_matrices)]
        counter2 += length(curr_matrices)
    end
    counter2 = 1
    for (curr_string, curr_combinations) in zip(cumulant_combinations, cumulant_index_combinations)
        curr_matrices, curr_minuses, curr_nums = make_mat(curr_combinations, curr_string, spin_order)
        cumulant_matrices[counter2:counter2+length(curr_matrices)-1] = curr_matrices
        cumulant_minus[counter2:counter2+length(curr_matrices)-1] = curr_minuses
        cumulant_num[counter2:counter2+length(curr_matrices)-1] = curr_nums
        cumulant_names[counter2:counter2+length(curr_matrices)-1] = [curr_string for i in 1:length(curr_matrices)]
        counter2 += length(curr_matrices)
    end
    dim::Int = 2^order
    I::Matrix{ComplexF64} = Matrix{ComplexF64}([1.0 0.0; 0.0 1.0])
    z::Matrix{ComplexF64} = Matrix{ComplexF64}([1.0 0.0; 0.0 -1.0])
    Ipz::Matrix{ComplexF64} = I + z
    constant_matrix::Matrix{ComplexF64} = kron(Ipz, [I for i in 1:spin_order]...) * 1 / 2^(spin_order + 1)
    matrix_generator::MatrixGenerator = MatrixGenerator(constant_matrix, linear_matrices, operator_minus, operator_num, operator_names, cumulant_matrices, cumulant_minus, cumulant_num, cumulant_names, dim)

    return density_matrices, matrix_generator
end
## Test 
#density_matrices, matrix_generator = Density_Matrix_Generator(3, system, samples)


##### Construct from Preparations ############################################################################################################

# Function to calculate the density matrix from solution, density_matrix and matrix_generator
function make_density_matrix(density_matrix_indexes::DensityMatrix, matrix_generator::MatrixGenerator, solution_at_T::DifferentialSolutionAtT{T}) where {T}
    # generate density matrix from solution and density_matrix_indexes
    rho::Matrix{T} = matrix_generator.constant_matrix
    # now iterate over all other terms and add them to the density matrix
    for i in 1:length(density_matrix_indexes.operator_indexes)
        curr_val = solution_at_T.vals[density_matrix_indexes.operator_indexes[i]]
        curr_mat = matrix_generator.operator_matrices[i]
        if matrix_generator.operator_num[i] # curr_num
            if abs(curr_val) > 1.0
                curr_val /= abs(curr_val)       # normalize
            end
            rho += curr_mat * curr_val
        elseif matrix_generator.operator_minus[i]  # curr_minus
            if abs(curr_val) > 1.0
                curr_val /= abs(curr_val)       # normalize
            end
            # mat *  val
            rho += curr_mat * curr_val
            # mat^\dagger * conj(val)
            rho += curr_mat' * conj(curr_val)
        else
            rho += curr_mat * curr_val
        end
    end
    # now add the cumulants 
    for i in 1:length(density_matrix_indexes.cumulant_indexes)  #cumulant_minus
        curr_val = solution_at_T.cumulant_vals[density_matrix_indexes.cumulant_indexes[i]]
        curr_mat = matrix_generator.cumulant_matrices[i]
        if matrix_generator.cumulant_num[i] # curr_num
            if abs(curr_val) > 1.0
                curr_val /= abs(curr_val)       # normalize
            end
            rho += curr_mat * curr_val
        elseif matrix_generator.cumulant_minus[i]  # curr_minus
            if abs(curr_val) > 1.0
                curr_val /= abs(curr_val)       # normalize
            end
            # mat *  val
            rho += curr_mat * curr_val
            # mat^\dagger * conj(val)
            rho += curr_mat' * conj(curr_val)
        else
            rho += curr_mat * curr_val
        end
    end
    return rho
end
## Test 
#rho = make_density_matrix(density_matrices[1], matrix_generator, solution)
#LinearAlgebra.tr(rho)




##### To check the inner workings 
function make_density_matrix_printing(density_matrix_indexes::DensityMatrix, matrix_generator::MatrixGenerator, solution_at_T::DifferentialSolutionAtT{T}) where {T}
    # generate density matrix from solution and density_matrix_indexes
    rho::Matrix{T} = matrix_generator.constant_matrix
    # now iterate over all other terms and add them to the density matrix
    for i in 1:length(density_matrix_indexes.operator_indexes)
        curr_val = solution_at_T.vals[density_matrix_indexes.operator_indexes[i]]
        curr_mat = matrix_generator.operator_matrices[i]
        # print the current term string and value
        println("$(matrix_generator.operator_names[i]): $(curr_val)")
        if matrix_generator.operator_num[i] # curr_num
            if abs(curr_val) > 1.0
                curr_val /= abs(curr_val)       # normalize
            end
            rho += curr_mat * curr_val
        elseif matrix_generator.operator_minus[i]  # curr_minus
            if abs(curr_val) > 1.0
                curr_val /= abs(curr_val)       # normalize
            end
            # mat *  val
            rho += curr_mat * curr_val
            # mat^\dagger * conj(val)
            rho += curr_mat' * conj(curr_val)
        else
            rho += curr_mat * curr_val
        end
    end
    # now add the cumulants 
    for i in 1:length(density_matrix_indexes.cumulant_indexes)  #cumulant_minus
        curr_val = solution_at_T.cumulant_vals[density_matrix_indexes.cumulant_indexes[i]]
        curr_mat = matrix_generator.cumulant_matrices[i]
        # print the current term string and value
        println("$(matrix_generator.cumulant_names[i]): $(curr_val)")
        if matrix_generator.cumulant_num[i] # curr_num
            if abs(curr_val) > 1.0
                curr_val /= abs(curr_val)       # normalize
            end
            rho += curr_mat * curr_val
        elseif matrix_generator.cumulant_minus[i]  # curr_minus
            if abs(curr_val) > 1.0
                curr_val /= abs(curr_val)       # normalize
            end
            # mat *  val
            rho += curr_mat * curr_val
            # mat^\dagger * conj(val)
            rho += curr_mat' * conj(curr_val)
        else
            rho += curr_mat * curr_val
        end
    end
    return rho
end