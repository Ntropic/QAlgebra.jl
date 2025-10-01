include("preprocessing.jl")
include("preprocessing2.jl")
#########################################################################################################################
#### Functions to generate system parameters and pulse variables ########################################################
#########################################################################################################################

function basis_wurst_pulse_value(t::Float64, T::Float64, n::Int, Lambda, phi0::Float64=0.0)
    # Lambda is the bandwidth of the pulse
    phi = Lambda * (t - t^2 / T) * 4 * pi + phi0
    val = (1 - abs(sin(pi * (t - T / 2) / T))^n) * exp(im * phi)
    return val
end
function basis_wurst_pulse_values(t::Vector{Float64}, T::Float64, n::Int, Lambda::Float64, phi0::Float64=0.0)
    # Vector of values
    vals = basis_wurst_pulse_value.(t, T, n, Lambda, phi0)
    return vals
end

function fourier_from_pulse_values(vals::Union{Vector{Float64},Vector{ComplexF64}}, t::Vector{Float64}, omegas::Vector{Float64})::Vector{ComplexF64}
    # For every omega in omegas compute the fourier coefficient
    coeffs::Vector{ComplexF64} = Vector{ComplexF64}(undef, length(omegas))
    for (i, omega) in enumerate(omegas)
        coeffs[i] = sum(vals .* exp.(-im * omega * t))
    end
    return coeffs
end

# Generators and Parameters
mutable struct Pulse_Param_Struct
    Amplitudes
    omegas::Vector{Float64}
    phi0::Float64
    T::Float64
    n::Int
    # initialize with default values
    function Pulse_Param_Struct(; Amplitudes=[1.0], omegas::Vector{Float64}=[1.0], phi0::Float64=0.0, T::Float64=1.0, n::Int=10)
        new(Amplitudes, omegas, phi0, T, n)
    end
end

# Create random initial pulse_parameters
function Random_Pulse_Param_Wurst(max_Amplitude=0.1, T::Float64=10, n::Int=10; omega_min::Float64=1.0, omega_max::Float64=-1.0, how_many::Int=1, rng=MersenneTwister(1234), normalize_omega::Bool=true) 
    omegas::Vector{Float64} = []
    if how_many == -1
        error("how_many must be specified")
    end
    if omega_max == -1
        if how_many > 1
            omegas = collect(range(omega_min, omega_min * (how_many - 1), length=how_many))
        else
            omegas = [omega_min]
        end
    else
        if how_many > 1
            omegas = collect(range(omega_min, omega_max, length=how_many))
        else
            omegas = [omega_min]
        end
    end
    if normalize_omega
        omegas = omegas/T
    end
    Amplitudes::Vector{Float64} = [(rand(rng) * 2 - 1) * max_Amplitude for i in 1:how_many]
    phi0::Float64 = 0.0
    return Pulse_Param_Struct(Amplitudes=Amplitudes, omegas=omegas, phi0=phi0, T=T, n=n)
end

function wurst_pulse_generator(pulse_params::Pulse_Param_Struct)::Function
    # Generate a Wurst pulse generating function from the parameters
    pp = pulse_params
    function fun(t::Float64)
        if t < 0 || t > pp.T
            return 0.0 + 0.0im
        end
        val::ComplexF64 = pp.Amplitudes[1] * basis_wurst_pulse_value(t, pp.T, pp.n, pp.omegas[1], pp.phi0)
        for i in 2:length(pp.Amplitudes)
            val += pp.Amplitudes[i] * basis_wurst_pulse_value(t, pp.T, pp.n, pp.omegas[i], pp.phi0)
        end
        return val
    end
    return fun
end
## Test
#pulse_param = Pulse_Param_Struct(Amplitudes=[1.0], omegas=[1.0])
#wurst_fun = wurst_pulse_generator(pulse_param)
#wurst_fun(0.5)

##### Fourier Basis (Sin Base) ##########################################################################################
function Random_Pulse_Param_Sin(max_Amplitude=0.1, T::Float64=10.0, n::Int=10; omega_min::Float64=-1.0, omega_max::Float64=-1.0, how_many::Int=1, rng=MersenneTwister(1234))
    Amplitudes = rand(rng, how_many) * 2 * max_Amplitude .- max_Amplitude
    if omega_min < 0.0
        omega_min = pi / T
        if omega_max < 0.0
            if how_many < 1
                error("how_many must be specified if omega_max is not specified")
            end
            omega_max = omega_min * how_many
        end
    elseif omega_max < 0.0
        omega_max = omega_min * how_many
    end
    omegas::Vector{Float64} = collect(range(omega_min, omega_max, length=how_many))
    phi0::Float64 = 0.0
    return Pulse_Param_Struct(Amplitudes=Amplitudes, omegas=omegas, phi0=phi0, T=T)
end

function sin_pulse_generator(pulse_params::Pulse_Param_Struct)::Function
    # Generate a sin pulse generating function from the parameters
    pp = pulse_params
    function fun(t::Float64)
        if t < 0 || t > pp.T
            return 0.0 + 0.0im
        end
        factor = (1 - abs(sin(pi * (t - pp.T / 2) / pp.T))^pp.n)
        val = pp.Amplitudes[1] * exp(-1im * (pp.omegas[1] * t + pp.phi0))
        for i in 2:length(pp.Amplitudes)
            val += pp.Amplitudes[i] * exp(-1im * (pp.omegas[i] * t + pp.phi0))
        end
        return val * factor
    end
    return fun
end



#########################################################################################################################
#### System Parameters at times t #######################################################################################
#########################################################################################################################

function num_args(f::Function)
    mthds = methods(f)
    for m in mthds
        sig = m.sig
        return length(sig.parameters) - 1
    end
end

function expand_and_multiply_array(array::Array, added::Vector)::Array
    # get dimensions of array
    dims = size(array)
    # create new array with added dimension
    new_dims::Vector{Int} = Vector{Int}(undef, length(dims) + 1)
    for i in 1:length(dims)
        new_dims[i] = dims[i]
    end
    new_dims[length(dims)+1] = length(added)
    new_array = Array{eltype(array)}(undef, Tuple(new_dims))
    # fill new array
    colons = ntuple(_ -> Colon(), length(dims))
    for i in 1:length(added)
        new_array[colons..., i] = array .* added[i]
    end
    return new_array
end
#M = Matrix{ComplexF64}([1. 2. 3.; 4. 5. 6.])
#added = [1, 2, 3]
#display(expand_and_multiply_array(matrix, added))
function value_generator(pulse_param, val_dict::Dict, system::SystemAndDicts)::Tuple{Function,Vector{Int}}
    vars_vec::Vector{Vector{String}} = system.vars_vec
    return value_generator(pulse_param, val_dict, vars_vec)
end
#function value_generator_remove_zeros(pulse_param, val_dict::Dict, system::SystemAndDicts)::Tuple{Function,Vector{Int}}
#    vars_vec::Vector{Vector{String}} = system.vars_vec
#    value_gen, constant_indexes = value_generator(pulse_param, val_dict, vars_vec)
#    # remove zero terms from system 
#    remove_indexes = find_zero_parameters(param_generator, constant_value_indexes; how_many_eps=how_many_eps)
#    for 
#    return value_gen, system
#end

function value_generator(pulse_param, val_dict::Dict, vars_vec::Vector{Vector{String}})::Tuple{Function,Vector{Int}}
    # Generates a generator function that constructs values of the variables (as Vector{ComplexF64}) at a given time t
    # vars_vec is --> inverse_vars_vec of preprocess_index_inversion function
    # user just needs to call generator with the time t 
    # returns a function that takes a time t and returns a vector of values
    T::Float64 = pulse_param.T
    values::Vector{ComplexF64} = ones(ComplexF64, length(vars_vec))
    curr_values::Vector{ComplexF64} = ones(ComplexF64, length(vars_vec))
    # store where function handles are, the value of those indexes without the function and the indexes to which they apply
    func_handles::Vector{Function} = Function[]
    func_args::Vector{Bool} = Bool[]
    func_handle_indexes::Vector{Int} = Int[]
    index_conjugations::Vector{Bool} = Bool[]
    index_sqrt::Vector{Bool} = Bool[]
    do_conj::Bool = false
    do_sqrt::Bool = false
    # if val_dict doesn't have key "", add it as value 0.0
    if !haskey(val_dict, "")
        val_dict[""] = 1.0 + 0.0im
    end
    # Start constructing
    for (i, vars) in enumerate(vars_vec)
        # get values from dictionary
        for j in 1:length(vars)
            curr_var = vars[j]
            # Determine if conjugate
            do_conj = false
            if occursin("^{*}", curr_var)
                do_conj = true
                curr_var = curr_var[1:end-4]  # remove string elements
            end
            do_sqrt = false
            if occursin("\\sqrt{", curr_var)
                do_sqrt = true
                curr_var = curr_var[7:end-1] # remove string elements
            end
            # Determine if function handle
            # check ig curr_var is in val_dict
            if !haskey(val_dict, curr_var)
                error("Variable $curr_var not in val_dict")
            end
            curr_val = val_dict[curr_var]
            if typeof(curr_val) <: Function
                # store function handle
                push!(func_handles, curr_val)
                push!(func_handle_indexes, i)
                push!(index_conjugations, do_conj)
                push!(index_sqrt, do_sqrt)
                if num_args(curr_val) == 1
                    push!(func_args, false)
                else
                    push!(func_args, true) # pass more arguments
                end
            else
                if do_sqrt
                    curr_val = sqrt(curr_val)
                end
                if do_conj
                    curr_val = conj(curr_val)
                end
                values[i] *= curr_val
            end
        end
    end
    curr_values = copy(values)
    # Construct the generator function
    function generator(t::Float64)
        # update values
        if t > T
            error("t > T")
        end
        # take values and copy to curr_values
        # update function handles
        curr_values = copy(values)
        curr_val::ComplexF64 = 0.0
        for (fun, args, ind, do_conj, do_sqrt) in zip(func_handles, func_args, func_handle_indexes, index_conjugations, index_sqrt)
            #curr_val::ComplexF64 = 0.0
            if !args
                curr_val = fun(t)
            else
                curr_val = fun(t, pulse_param)
            end
            if do_sqrt
                curr_val = sqrt(curr_val)
            end
            if do_conj
                curr_val = conj(curr_val)
            end
            curr_values[ind] *= curr_val
        end
        return curr_values
    end
    # go through the indexes, if an index has no function handle, then it is a constant value
    constant_value_indexes::Vector{Int} = Int[]
    for i in 1:length(values)
        if !(i in func_handle_indexes)
            push!(constant_value_indexes, i)
        end
    end
    return generator, constant_value_indexes
end
## Test
#how_many::Int = 10
#max_A::Float64 =  0.1
#T::Float64 = 10.0
#n = 10
#pulse_param::Pulse_Param_Struct = Random_Pulse_Param_Wurst(how_many, max_A, T, n, rng=rng)
#val_dict::Dict = Dict("\\beta" => wurst_pulse_generator(pulse_param), "\\Delta" => 0.1, "g" => [0.2,0.21,0.22,0.23,0.24,0.25], "\\gamma" => 0.1, "\\Gamma" => 0.1, "\\kappa" => 0.1)
#param_generator = value_generator(pulse_param, val_dict, vars_vec)
