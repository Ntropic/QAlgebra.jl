include("preprocessing2.jl")
include("distribution_sampling.jl")
include("indexing.jl")

# Helper scripts
function haskeys(dict::Dict, keys::Array)::Bool
    # checks if a dictionary has all the keys in the array
    for key in keys
        if !haskey(dict, key)
            return false
        end
    end
    return true
end
# generate arguments for a function from a dictionary and an argument_name list 
function generate_args(dict::Dict, arg_names::Array)::Array
    args::Array = []
    for arg_name in arg_names
        push!(args, dict[arg_name])
    end
    return args
end

#############################################################################################################
####### Initial States ######################################################################################
function homogenize_and_normalize_parameter(phi::Union{Int64,Float64,ComplexF64,Vector{Float64},Vector{ComplexF64}}, coeff_num::Int, normalize::Bool=true, name::String="phi")::Vector{ComplexF64}
    # phi is a parameter, turn into a vector if necessary, make complex and pad zeros. normalize if amplitude too large unless bool param is false
    phi2::Vector{ComplexF64} = zeros(ComplexF64, coeff_num)
    if isa(phi, Vector)
        if length(phi) != coeff_num
            error("phi (has length $(length(phi))) must have the same length as the number of coefficients (has length $coeff_num)")
        end
        if !isa(phi, Vector{ComplexF64})
            phi2 = phi .+ 0.0im    # make sure phi is complex
        end
    else
        phi2[1] = phi + 0.0im
    end
    if normalize
        abs_phi = abs.(phi2)
        max_abs_phi = maximum(abs_phi)
        if max_abs_phi > 1.0
            println("Warning: Max. abs. value of $(name) is greater than 1.0, normalizing to 1.0.")
            phi2 = phi2 / max_abs_phi
        end
    end
    return phi2
end
function single_excitation_state_gen_distributed(N::Int, samples::SamplesAndWeights, A::Union{Int64,Float64,ComplexF64,Vector{Float64},Vector{ComplexF64}}, phi::Union{Int64,Float64,ComplexF64,Vector{Float64},Vector{ComplexF64}}; non::Bool=false)::Function
    # if none, then calculate the single non excitation state instead
    # Generates the state 1/\sqrt{N} sum_{i=1}^N \sigma_x^i |0...0>
    # phi specifies the phases in the space of function basis (transformaed via coeffs2vals)
    #           see samples.order_comb to associate basis function terms with phis vector elements
    # where N is the number of spins
    # and n_z is the number of spins in the z-direction in the operator
    # Calculates <\prod_{i=1}^n \sigma_z^{k_i}>:
    coeffs2vals = samples.coeffs2vals
    sample_num::Int = size(coeffs2vals, 1)
    coeff_num::Int = size(coeffs2vals, 2)
    # check if phi is vector subtype -> if it is check length
    phi::Vector{ComplexF64} = homogenize_and_normalize_parameter(phi, coeff_num, true, "phi")
    A::Vector{ComplexF64} = homogenize_and_normalize_parameter(A, coeff_num, true, "A")

    z_pref = (-1 + 2 / N)
    z_amp = z_pref .* A

    z_amp_by_index::Vector{ComplexF64} = coeffs2vals * z_amp
    max_phase_amp = sqrt.(1 .- abs2.(z_amp_by_index))
    phase_by_index::Vector{ComplexF64} = coeffs2vals * phi
    x_phase_by_index::Vector{Float64} = real(phase_by_index) .* max_phase_amp
    y_phase_by_index::Vector{Float64} = imag(phase_by_index) .* max_phase_amp

    function generator(op_indexes::Vector{Vector{Int}})::ComplexF64
        op_exp::ComplexF64 = 1.0 + im * 0.0
        # how many times does z occur?
        for xi in op_indexes[1]
            op_exp *= x_phase_by_index[xi]
        end
        for yi in op_indexes[2]
            op_exp *= y_phase_by_index[yi]
        end
        for zi in op_indexes[3]
            op_exp *= z_amp_by_index[zi]
        end
        return op_exp
    end
    function generator_non(op_indexes::Vector{Vector{Int}})::ComplexF64
        op_exp::ComplexF64 = 1.0 + im * 0.0
        # how many times does z occur?
        for xi in op_indexes[1]
            op_exp *= x_phase_by_index[xi]
        end
        for yi in op_indexes[2]
            op_exp *= y_phase_by_index[yi]
        end
        for zi in op_indexes[3]
            op_exp *= -z_amp_by_index[zi]
        end
        return op_exp
    end
    if non
        return generator_non
    else
        return generator
    end
end

# This function can be made more elaborate to account for g and \Delta dependent states 
function single_excitation_state_gen(N::Int; non::Bool=false)::Function
    # if none, then calculate the single non excitation state instead
    # Generates the state 1/\sqrt{N} sum_{i=1}^N \sigma_x^i |0...0>
    # where N is the number of spins
    # and n_z is the number of spins in the z-direction in the operator
    # Calculates <\prod_{i=1}^n \sigma_z^{k_i}>:
    function generator(op_indexes::Vector{Vector{Int}})::ComplexF64
        op_exp::ComplexF64 = 0.0 + im * 0.0
        # how many times does z occur?
        n_x::Int = length(op_indexes[1])
        n_y::Int = length(op_indexes[2])
        n_z::Int = length(op_indexes[3])
        op_exp = (-1)^n_z * (1 - 2 * n_z / N)
        if n_x + n_y > 0
            op_exp = 0.0 + im * 0.0
        end
        return op_exp
    end
    function generator_non(op_indexes::Vector{Vector{Int}})::ComplexF64
        op_exp::ComplexF64 = 0.0 + im * 0.0
        # how many times does z occur?
        n_x::Int = length(op_indexes[1])
        n_y::Int = length(op_indexes[2])
        n_z::Int = length(op_indexes[3])
        op_exp = (1 - 2 * n_z / N)
        if n_x + n_y > 0
            op_exp = 0.0 + im * 0.0
        end
        return op_exp
    end
    if non
        return generator_non
    else
        return generator
    end
end

function gibbs_state_gen(p::Float64)::Function
    # generator, that returns a function that computes the expectation value for a specified number of fcreation and annihilation operators
    # p is the occupation probability of the ground state
    # Theta is the phase of the coherent state
    x::Float64 = 0
    y::Float64 = 0
    z::Float64 = 1 - 2 * p
    function generator(op_indexes::Vector{Vector{Int}})::ComplexF64
        # if any x or y in string return 0.0
        op_exp::ComplexF64 = 0.0 + im * 0.0
        n_x::Int = length(op_indexes[1])
        n_y::Int = length(op_indexes[2])
        n_z::Int = length(op_indexes[3])
        op_exp = z^n_z
        if n_x + n_y > 0
            op_exp = 0.0 + im * 0.0
        end
        return op_exp
    end
    return generator
end

function coherent_state_gen(alpha::Float64, Theta::Float64)::Function
    # generator, that returns a function that computes the expectation value for a specified number of fcreation and annihilation operators
    # alpha is the amplitude of the coherent state
    # Theta is the phase of the coherent state
    abs_alpha::Float64 = abs(alpha)
    e_itheta::ComplexF64 = exp(im * Theta)
    function generator(cavity_operators::Vector{Int})::ComplexF64
        creation::Int = cavity_operators[1]
        annihilation::Int = cavity_operators[2]
        return abs_alpha^(creation + annihilation) * e_itheta^(creation - annihilation)
    end
    return generator
end

function operators2initialstates(system::SystemAndDicts, samples::SamplesAndWeights, spin_type::String, cavity_type::String, spin_param::Dict, cavity_param::Dict)::Vector{ComplexF64}
    # Generate initial conditions for the Equations
    # spin_type is either "gibbs" or "single_excitation", "single_non_excitation", "functional_single_excitation", "functional_single_non_excitation" (update with new functions)
    # cavity_type currently_only_supports "coherent" (update with new functions)
    initial_values::Vector{ComplexF64} = Vector{ComplexF64}(undef, system.how_many_total)
    operator_strings_dict::Dict{String,Op_Group_Type} = system.operator_strings_dict
    all_combinations_list = system.all_combinations_list
    if cavity_type == "coherent"
        coherent_args::Vector{String} = ["alpha", "Theta"]
        if !haskeys(cavity_param, coherent_args)
            error("cavity_param must have keys alpha and Theta for coherent states")
        end
        cavity_fun = coherent_state_gen(generate_args(cavity_param, coherent_args)...)
    else
        error("cavity_type unknown")
    end
    if spin_type == "gibbs"
        if !haskey(spin_param, "p")
            error("spin_param must have key p for gibbs states")
        end
        spin_fun = gibbs_state_gen(spin_param["p"])
    elseif spin_type == "single_excitation" || spin_type == "single_non_excitation"
        if !haskey(spin_param, "N")
            error("spin_param must have key N for single (non) excitation states")
        end
        if !haskey(spin_param, "phi") && !haskey(spin_param, "A")
            error("spin_param must have key phi or A for single (non) excitation states")
        end
        coeffs2vals = samples.coeffs2vals
        if haskey(spin_param, "phi")
            phi = spin_param["phi"]
        else
            phi = zeros(ComplexF64, size(coeffs2vals, 2))
        end
        if haskey(spin_param, "A")
            A = spin_param["A"]
        else
            A = zeros(ComplexF64, size(coeffs2vals, 2))
            A[1] = 1.0
        end
        if spin_type == "single_non_excitation"
            spin_fun = single_excitation_state_gen_distributed(spin_param["N"], samples, A, phi; non=true)
        else
            spin_fun = single_excitation_state_gen_distributed(spin_param["N"], samples, A, phi)
        end
    else
        error("spin_type unknown. Options are: gibbs, single_excitation, single_non_excitation")
    end
    for (opstr, op_group_type) in operator_strings_dict
        cavity_operators = [count(c -> c == '+', opstr), count(c -> c == '-', opstr)]
        op_spin_orders = op_group_type.op_spin_orders
        first_index = op_group_type.first_index
        curr_how_many = op_group_type.how_many
        all_comb = collect(generate_all_combinations_of_combinations(op_spin_orders, all_combinations_list))
        @threads for i in 1:curr_how_many
            op_indexes = indexes_to_index_vector(all_comb[i], op_spin_orders)
            initial_values[first_index+i-1] = spin_fun(op_indexes) * cavity_fun(cavity_operators)
        end
    end
    return initial_values
end
## Test 
#term_str = ["z", "zz", "zzz", "zzzz"]
#spin_param = Dict("N" => 10^6)
#cavity_param = Dict("alpha" => 1.0, "Theta" => 0.0)
#initial_values = term_str2states(term_str, "single_excitation", "coherent", spin_param, cavity_param)

#############################################################################################################################################
#### Phase & Amplitude Selectors ############################################################################################################
#############################################################################################################################################

# A function that takes samples as an input aswell as orders or weighted orders and returns an index vector 
# use arbitrarily many order and weighted orders as input via args... 
# The function returns a vector of indices that can be used to access the samples for initial state generation
# order 0 is constant, 1, 3, 5... are sin and 2, 4, 6... are cos, make a verctor for the different dimensions of the samples (e.g. g, Delta, ...)
function weighted_vector(samples::SamplesAndWeights, args...)
    order_comb::Vector{Vector{Float64}} = samples.order_comb
    # for every arg in args check if it is a tuple, vector or single value and then use it to generate the index vector
    index_vec::Vector{Float64} = zeros(length(order_comb))
    for (i, arg) in enumerate(args)
        if isa(arg, Tuple)
            weight = arg[1]
            order = arg[2]
        elseif isa(arg, Vector)
            weight = 1.0
            order = arg
        end
        if !isa(order, Vector{Int})
            if isa(order, Vector) # try making it int 
                order_old = copy(order)
                order = convert(Vector{Int}, order_old)
                if order != order_old
                    error("Order has to be a vector of integers. Got: $order_old")
                end
            else
                error("Order has to be a vector of integers. Got: $order")
            end
        end
        # find order in order_comb
        ind = findfirst(x -> x == order, order_comb)
        if ind == nothing
            error("Order not found in order_comb. Got: $order")
        end
        index_vec[ind] = weight
    end
    return index_vec
end
## Test 
# weighted_vector(samples, (1.0, [0,1]), [1,0])

function even_odd_vector(samples::SamplesAndWeights, type::Symbol=:even, type2::Symbol=:real; rng::Union{MersenneTwister,Bool}=false)::Vector{ComplexF64}
    # type takes :even (default) or :odd to generate even or odd or :complex or :real (default) to specify
    if type2 == :even || type2 == :odd
        type, type2 = type2, type
    end
    if type != :even && type != :odd
        error("Type has to be :even or :odd. Got: $type")
    end
    if type2 != :real && type2 != :complex
        error("Type2 has to be :real or :complex. Got: $type2")
    end
    order_comb::Vector{Vector{Float64}} = samples.order_comb
    do_rng = false
    if !isa(rng, Bool)
        do_rng = true
    end
    index_vec::Vector{ComplexF64} = zeros(ComplexF64, length(order_comb))
    for (i, order) in enumerate(order_comb)
        orders_odd::Vector{Bool} = [x % 2 == 1 for x in order]
        any_odd::Bool = any(orders_odd)
        if type == :even
            if !any_odd
                if !do_rng
                    index_vec[i] = 1.0
                else
                    index_vec[i] = rand(rng)
                    if type2 == :complex
                        index_vec[i] += im * rand(rng)
                    end
                end
            end
        elseif type == :odd
            if any_odd
                if !do_rng
                    index_vec[i] = 1.0
                else
                    index_vec[i] = rand(rng)
                    if type2 == :complex
                        index_vec[i] += im * rand(rng)
                    end
                end
            end
        end
    end
    # Normalize
    index_vec = index_vec / sum(index_vec)
    return index_vec
end
## Test
#even_odd_vector(samples, :even, rng=MersenneTwister(1234))