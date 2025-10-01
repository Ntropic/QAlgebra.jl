using QuadGK
using LinearAlgebra
using Random
using Base.Threads
using FileIO
include("combinatorics.jl")
include("distribution_sampling.jl")

#### Fourier and weighted Fourier Integrals ####
function fourier_fun_gen(order::Int, range::Vector{Float64})
    # orders are [const, sin_1, cos_1, sin_2, cos_2, ...]
    if order == 0
        return x -> 1.0
    end
    fourier_order = (order + 1) ÷ 2
    do_sin = order % 2 == 1
    if do_sin
        return x -> sin(2 * pi * fourier_order * x / (range[2] - range[1]))
    else
        return x -> cos(2 * pi * fourier_order * x / (range[2] - range[1]))
    end
end
function fourier_exp2val(x::Float64, order::Int, range::Vector{Float64})
    if order == 0
        return 1.0
    end
    fourier_order = (order + 1) ÷ 2
    do_sin = order % 2 == 1
    if do_sin
        return sin(2 * pi * fourier_order * x / (range[2] - range[1]))
    else
        return cos(2 * pi * fourier_order * x / (range[2] - range[1]))
    end
end
function fourier_exp2vals(x::Vector{Float64}, order::Int, range::Vector{Float64})
    if order == 0
        return ones(Float64, length(x))
    end
    fourier_order = (order + 1) ÷ 2
    do_sin = order % 2 == 1
    if do_sin
        return sin.(2 * pi * fourier_order * x / (range[2] - range[1]))
    else
        return cos.(2 * pi * fourier_order * x / (range[2] - range[1]))
    end
end
# Test 
#fig = Figure()
#ax = Axis(fig[1, 1])
#x = range(-1.0, 1.0, length=1000)
#for order in 1:3
#    f = fourier_fun_gen(order, [-1.0, 1.0])
#    lines!(ax, x, f.(x), label="Order: $order")
#end
#fig

function integrate_weighted_fouriers(fun::Function, range::Vector{Float64}, order::Int; reltol::Float64=1e-6, abstol::Float64=1e-10)
    if length(range) != 2
        error("range must be a vector of length 2")
    end
    p = fourier_fun_gen(order, range)
    norm = quadgk(x -> p(x) * fun(x), range[1], range[2], rtol=reltol, atol=abstol)[1]
    return norm
end
function integrate_to_order(fun::Function, range::Vector{Float64}, max_order::Int; reltol::Float64=1e-6, abstol::Float64=1e-10)
    norms::Vector{Float64} = Vector{Float64}(undef, max_order + 1)
    for i in 0:max_order
        norms[i+1] = integrate_weighted_fouriers(fun, range, i, reltol=reltol, abstol=abstol)
    end
    # remove close to zero 
    for i in 1:length(norms)
        if abs(norms[i]) < abstol
            norms[i] = 0.0
        end
    end
    return norms
end
## Test with gaussian fun
#fun = x -> exp(-x^2)
#integrate_to_order(fun, [-4.0, 4.0], 6, 0.0, reltol=10^-15)

#### Generalize to products of basis functions
function integrate_weighted_fouriers(fun::Function, range::Vector{Float64}, order::Int, order2::Int; reltol::Float64=1e-6, abstol::Float64=1e-10)
    if length(range) != 2
        error("range must be a vector of length 2")
    end
    p = fourier_fun_gen(order, range)
    q = fourier_fun_gen(order2, range)
    norm = quadgk(x -> p(x) * q(x) * fun(x), range[1], range[2], rtol=reltol, atol=abstol)[1]
    return norm
end
function integrate_to_order_mat(fun::Function, range::Vector{Float64}, max_order::Int; reltol::Float64=1e-6, abstol::Float64=1e-10)
    norms::Matrix{Float64} = Matrix{Float64}(undef, max_order + 1, max_order + 1)
    for i in 0:max_order
        norms[i+1, i+1] = integrate_weighted_fouriers(fun, range, i, i, reltol=reltol, abstol=abstol)
        for j in i+1:max_order
            norms[i+1, j+1] = integrate_weighted_fouriers(fun, range, i, j, reltol=reltol, abstol=abstol)
            norms[j+1, i+1] = norms[i+1, j+1]
        end
    end
    # remove close to zero 
    for i in 1:max_order+1
        for j in 1:max_order+1
            if abs(norms[i, j]) < abstol
                norms[i, j] = 0.0
            end
        end
    end
    return norms
end

#### Matrix Based Linear Regression ####
function fit_matrix(coeffs2vals::Matrix{Float64}, weight_vector::Vector{Float64}=Vector{Float64}(); ridge_lambda::Float64=1e-13)
    if length(weight_vector) == 0
        div_mat = coeffs2vals' * coeffs2vals - ridge_lambda*Diagonal(ones(size(coeffs2vals, 2)))
        return div_mat \ (coeffs2vals')                 # use pseudoinverse
        #return inv(coeffs2vals' * coeffs2vals) * coeffs2vals'
    else
        weights::Diagonal{Float64, Vector{Float64}} = Diagonal(weight_vector)
        div_mat = (coeffs2vals' * weights * coeffs2vals) - ridge_lambda*Diagonal(ones(size(coeffs2vals, 2)))
        return div_mat \ (coeffs2vals' * weights)
        #return inv(coeffs2vals' * weights * coeffs2vals) * coeffs2vals' * weights
    end
end
function n_d_fourier_val(location::Vector{Float64}, orders::Vector{Int}, ranges::Vector{Vector{Float64}})
    # Generate value of fourier with orders 
    # for each direction (product)
    res = 1.0
    for i in 1:length(location)
        res *= fourier_exp2val(location[i], orders[i], ranges[i])
    end
    return res
end

function make_coeffs2vals_matrix(locations::Vector{Vector{Float64}}, order_comb::Vector{Vector{Int}}, ranges::Vector{Vector{Float64}})
    # Generate values 
    n_orders::Int = length(order_comb)
    n_locations::Int = length(locations)
    normalisation::Float64 = 1#/(2/pi * sqrt(n_locations)) # normalisation factor
    coeffs2vals::Matrix{Float64} = zeros(n_locations, n_orders) # coefficients for each order to values at locations (via matrix*vector : coeffs2vals*coeffs = values)
    for i in 1:n_locations
        for j in 1:n_orders
            coeffs2vals[i, j] = n_d_fourier_val(locations[i], order_comb[j], ranges) * normalisation
        end
    end
    return coeffs2vals
end

# generate fit matrix from location vectors and order vectors (for polynomial orders in the directions)
function n_d_fourier_fit(locations::Vector{Vector{Float64}}, order_comb::Vector{Vector{Int}}, ranges::Vector{Vector{Float64}}, weight_vector::Vector{Float64}=Vector{Float64}())::Tuple{Matrix{Float64},Matrix{Float64}}
    # Generate values 
    coeffs2vals = make_coeffs2vals_matrix(locations, order_comb, ranges)
    # for every coefficient determine the norm 
    # use linear regression to find coefficients from values, determine regression matrix via 
    vals2coeffs::Matrix{Float64} = fit_matrix(coeffs2vals, weight_vector)
    return vals2coeffs, coeffs2vals
end

#### Position Sampling ####
function make_positions_vec(fun_vec::Vector{Function}, range_vec::Vector{Vector{Float64}}, n_sample_vec::Vector{Int}; rng::MersenneTwister=MersenneTwister(1234), reltol::Float64=1e-15, abstol::Float64=1e-15, on_borders::Bool=false, type::Symbol=:equidistant)
    sample_pos::Vector{Vector{Float64}} = []
    dir_num = length(fun_vec)
    for i in 1:dir_num
        if type==:equidistant
            curr_pos, curr_center = equidistant_sample(range_vec[i], n_sample_vec[i], on_borders=on_borders, reorder=true)
        elseif type==:random
            curr_pos, curr_center = sample_from_distribution(fun_vec[i], range_vec[i], n_sample_vec[i], reltol=reltol, abstol=abstol, rng=rng, reorder=true, threaded=true)
        elseif type == :equiprobable
            curr_pos, curr_center = equidistant_sample_from_distribution(fun_vec[i], range_vec[i], n_sample_vec[i], on_borders=on_borders, reorder=true)
        else 
            error("type must be :equidistant, :random or :equiprobable")
        end
        push!(sample_pos, curr_pos)
    end
    return sample_pos
end

function sample_grid_positions(fun_vec::Vector{Function}, range_vec::Vector{Vector{Float64}}, n_sample_vec::Vector{Int}; rng::MersenneTwister=MersenneTwister(1234), reltol::Float64=1e-6, abstol::Float64=1e-10, on_borders::Bool=false, type::Symbol=:equidistant)
    # Generate n samples for each direction in range_vec
    dir_num = length(fun_vec)
    if dir_num != length(range_vec)
        error("fun_vec, range_vec must have the same length")
    end
    if dir_num != length(n_sample_vec)
        error("fun_vec, n_sample_vec must have the same length")
    end
    sample_pos = make_positions_vec(fun_vec, range_vec, n_sample_vec, rng=rng, reltol=reltol, abstol=abstol, on_borders=on_borders, type=type)
    return sample_pos
end

function real_i_s(i)
    # reorder i 
    rest = (i - 1) % 2
    val = i ÷ 2
    if rest == 1
        return val
    else
        return -val
    end
end
function cubed_sample_grid(fun_vec::Vector{Function}, range_vec::Vector{Vector{Float64}}, orders_vec::Union{Vector{Int},Int}=Vector{Int}(); on_borders::Bool=true, reltol::Float64=1e-15, abstol::Float64=1e-15, extra_points::Union{Int,Vector{Int}}=0, type::Symbol=:equidistant, rng::MersenneTwister=MersenneTwister(1234))
    # Generate an n-dimensional sample grid 
    # type can be :equidistant or :random or :equiprobable or :randomprob
    dir_num = length(fun_vec)
    if typeof(orders_vec) == Int
        orders_vec = [orders_vec for i in 1:dir_num]
    end
    if dir_num != length(range_vec)
        error("fun_vec, range_vec must have the same length")
    end
    n_sample_vec = 2 * orders_vec .+ 1 # zeroth order  and sin and cos terms each 
    n_sample_vec_extra = n_sample_vec .+ extra_points
    if any(orders_vec .> n_sample_vec .- 1)
        error("orders_vec must be smaller than n_sample_vec - 1")
    end
    locations::Vector{Vector{Float64}} = []
    if type == :random || type == :randomprob
        total_n::Int = prod(n_sample_vec_extra)
        if type == :randomprob
            curr_locations = sample_from_distributions(fun_vec, range_vec, total_n, reltol=reltol, abstol=abstol, rng=rng, threaded=true, reorder=false)[1]
        else # :random
            curr_locations = sample_from_intervals(range_vec, total_n, rng=rng)
        end
        for counter in 1:total_n
            push!(locations, [curr_locations[i][counter] for i in 1:dir_num])
        end
    else
        sample_pos = sample_grid_positions(fun_vec, range_vec, n_sample_vec_extra, rng=rng, reltol=reltol, abstol=abstol, on_borders=on_borders, type=type)
        # make a grid of all combinations
        for i_s in mrange(n_sample_vec_extra, startpoint=1)
            push!(locations, [sample_pos[i][i_s[i]] for i in 1:dir_num])
        end
    end

    # make orders vector
    order_comb::Vector{Vector{Int}} = []
    for i_s in mrange(n_sample_vec, startpoint=1)
        push!(order_comb, i_s .- 1)
    end
    return locations, order_comb
end


function basis_function_integrals(fun_vec::Vector{Function}, range_vec::Vector{Vector{Float64}}, orders_vec::Vector{Vector{Int}}; reltol::Float64=1e-6, abstol=1e-10)
    integrals = Vector{Vector{Float64}}()
    dir_num = length(fun_vec)
    for i in 1:dir_num
        max_order_i = maximum([ord[i] for ord in orders_vec])
        curr_integrals = integrate_to_order(x -> fun_vec[i](x), range_vec[i], max_order_i, reltol=reltol, abstol=abstol)
        push!(integrals, curr_integrals)
    end
    return integrals
end
function basis_function_integrals(fun_vec::Vector{Function}, range_vec::Vector{Vector{Float64}}, max_orders::Vector{Int}; reltol::Float64=1e-6, abstol=1e-10)
    integrals = Vector{Vector{Float64}}()
    dir_num = length(fun_vec)
    for i in 1:dir_num
        max_order_i = max_orders[i]
        curr_integrals = integrate_to_order(x -> fun_vec[i](x), range_vec[i], max_order_i, reltol=reltol, abstol=abstol)
        push!(integrals, curr_integrals)
    end
    return integrals
end
function basis_function_integrals_mat(fun_vec::Vector{Function}, range_vec::Vector{Vector{Float64}}, orders_vec::Vector{Vector{Int}}; reltol::Float64=1e-6, abstol=1e-10)
    integrals = Vector{Matrix{Float64}}()
    dir_num = length(fun_vec)
    for i in 1:dir_num
        max_order_i = maximum([ord[i] for ord in orders_vec])
        curr_integrals = integrate_to_order_mat(x -> fun_vec[i](x), range_vec[i], max_order_i, reltol=reltol, abstol=abstol)
        push!(integrals, curr_integrals)
    end
    return integrals
end
function basis_function_integrals_mat(fun_vec::Vector{Function}, range_vec::Vector{Vector{Float64}}, max_orders::Vector{Int}; reltol::Float64=1e-6, abstol=1e-10)
    integrals = Vector{Matrix{Float64}}()
    dir_num = length(fun_vec)
    for i in 1:dir_num
        max_order_i = max_orders[i]
        curr_integrals = integrate_to_order_mat(x -> fun_vec[i](x), range_vec[i], max_order_i, reltol=reltol, abstol=abstol)
        push!(integrals, curr_integrals)
    end
    return integrals
end

function weights_and_errors(fun_vec::Union{Vector{Function},Vector{Vector{Float64}}}, locations::Vector{Vector{Float64}}, order_comb::Vector{Vector{Int}}, range_vec::Vector{Vector{Float64}}; fun_return::Bool=true, reltol::Float64=1e-15, weighted::Bool=true)
    if length(fun_vec) != length(range_vec)
        error("fun_vec and range_vec must have the same length")
    end
    dir_num = length(order_comb[1])
    weighted_integrals::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    if typeof(fun_vec) == Vector{Function}
        weighted_integrals = basis_function_integrals(fun_vec, range_vec, order_comb, reltol=1e-6)
        # location_fit_weights are given by probabilities at the locations
        if !weighted
            location_fit_weights = ones(length(locations))   # the integrals themselves should do this weighting for us. 
        else
            location_fit_weights::Vector{Float64} = zeros(length(locations))
            for i in 1:length(locations)
                location_fit_weights[i] = prod([fun_vec[j](locations[i][j]) for j in 1:dir_num])   # weight the locations importance by their probability amplitude!
            end
        end
    else
        weighted_integrals = fun_vec
        location_fit_weights = ones(length(locations))  # weight the locations importance equally 
    end
    integrals_vector = zeros(length(order_comb))
    for i in 1:length(order_comb)
        integrals_vector[i] = prod([weighted_integrals[j][order_comb[i][j]+1] for j in 1:length(order_comb[i])])
    end
    if weighted
        vals2coeffs, coeffs2vals = n_d_fourier_fit(locations, order_comb, range_vec, location_fit_weights)
    else
        vals2coeffs, coeffs2vals = n_d_fourier_fit(locations, order_comb, range_vec)
    end
    # Determine the weights for the HEOM's
    weights::Vector{Float64} = vec(integrals_vector' * vals2coeffs)
    val_estimator = coeffs2vals * vals2coeffs #  element wise error is err=vals-vals_estimator*vals 
    cov_mat = inv(coeffs2vals' * coeffs2vals) # covariance matrix of the coefficients so that cov = sigma^2*cov_mat where sigma^2 = sum(err.^2)
    dof = length(locations) - length(order_comb)
    cov_mat = cov_mat ./ dof # normalize by degrees of freedom
    function get_error(values::Union{Vector{Float64},Vector{ComplexF64}})
        err_2 = abs(values - val_estimator * values) .^ 2
        sigma_2 = sum(err_2) / dof
        #covariance = inv(coeffs2vals' * coeffs2vals) * sigma_2
        return sqrt.(sigma_2), sqrt.(diag(cov_mat)) # std of samples and paramaters respectively
    end
    if fun_return
        return weights, get_error, vals2coeffs, coeffs2vals, integrals_vector
    else
        return weights, val_estimator, vals2coeffs, coeffs2vals, integrals_vector, cov_mat, order_comb
    end
end

# Add variance matrices 
function Ex_and_Varx(weighted_integrals::Vector{Vector{Float64}}, weighted_integrals_mat::Vector{Matrix{Float64}}, locations::Vector{Vector{Float64}}, order_comb::Vector{Vector{Int}}, range_vec::Vector{Vector{Float64}}; reltol::Float64=1e-6, abstol::Float64=1e-10, eps::Float64=1e-12, loadsave::Bool=false, printing::Bool=false)
    # Construct vals2coeffs to determine the expansion coefficients of the basis functions,
    # weights to directly compute the expectation value as a linear combination,
    #and G_ij, so that coeff' * G_ij * coeff = Var(x) to calculate the variance of x around the expectation value
    if length(weighted_integrals) != length(range_vec)
        error("weighted_integrals and range_vec must have the same length")
    end
    if length(weighted_integrals) != length(weighted_integrals_mat)
        error("weighted_integrals and weighted_integrals_mat must have the same length")
    end
    dir_num = length(order_comb[1])
    width_vec = [(range_vec[i][2] - range_vec[i][1]) / 2 for i in 1:dir_num]

    max_orders::Vector{Int} = [maximum([ord[i] * 2 for ord in order_comb]) for i in 1:dir_num]
    # make a hash out of weighted_integrals, locations, order_comb, range_vec
    if loadsave
        hash_code = hash((weighted_integrals, weighted_integrals_mat, locations, order_comb, range_vec))
        #check if file Cached/Ex_and_Varx_****.jld2 exists, if yes, load it, else calculate and save it 
        file_name = joinpath("Cached", "Ex_and_Varx_$(hash_code).jld2")
        if isfile(file_name)
            if printing
                println("Loading from file ", file_name)
            end
            return load(file_name, "vals2coeffs"), load(file_name, "E_weights"), load(file_name, "G_ij")
        end
    end
    # Create all unique combinations of order_comb terms (order_comb[i] + order_comb[j]) 
    n_comb::Int = length(order_comb)

    integrals_vector = zeros(length(order_comb))
    for i in 1:length(order_comb)
        integrals_vector[i] = prod([weighted_integrals[j][order_comb[i][j]+1] for j in 1:dir_num])
    end
    # remove components with abs(integrals) < eps   ### removed this as it messes with the standard deviation calculation
    #keep_indexes = findall(abs.(integrals_vector) .> eps)
    #red_integrals_vector = integrals_vector[keep_indexes]
    #red_n_comb = length(keep_indexes)
    #red_order_comb = order_comb[keep_indexes]

    vals2coeffs = n_d_fourier_fit(locations, order_comb[1:n_comb], range_vec)[1]
    #red_vals2coeffs = vals2coeffs[keep_indexes, :]
    # Determine the weights for the HEOM's
    E_weights::Vector{Float64} = vec(integrals_vector[1:n_comb]' * vals2coeffs) # weights for the expectation value
    F_ij::Matrix{Float64} = Matrix{Float64}(undef, n_comb, n_comb)
    @threads for i in 1:n_comb
        for j in i:n_comb
            curr_order1::Vector{Int} = order_comb[i] .+ 1
            curr_order2::Vector{Int} = order_comb[j] .+ 1
            curr_integral::Float64 = prod([weighted_integrals_mat[k][curr_order1[k], curr_order2[k]] for k in 1:dir_num])
            F_ij[i, j] = curr_integral
            F_ij[j, i] = F_ij[i, j]
        end
    end
    G_ij::Matrix{Float64} = F_ij - integrals_vector[1:n_comb] * integrals_vector[1:n_comb]' # c' * G_ij * c = Var(x)
    if loadsave
        if printing
            println("Saving to file ", file_name)
        end
        save(file_name, "vals2coeffs", vals2coeffs, "E_weights", E_weights, "G_ij", G_ij)
    end
    return vals2coeffs, E_weights, G_ij
end

# Combine variables in structs 
struct SamplesAndWeights
    locations::Vector{Vector{Float64}}
    prob_fun::Vector{Function}
    range_vec::Vector{Vector{Float64}}
    vals2coeffs::Matrix{Float64}
    coeffs2vals::Matrix{Float64}
    weights_dict::Dict{String,Vector{Float64}}
    order_comb::Vector{Vector{Int}}
    integrals_vector::Vector{Float64}
    sample_vars::Dict{String,Vector{ComplexF64}}
    sample_num::Int
    var_strs::Vector{String}
    has_std_dict::Bool
    order_to_exp_std_dict::Dict{Vector{Int}, Tuple{Matrix{Float64}, Vector{Float64}, Matrix{Float64}}}
    function SamplesAndWeights(locations::Vector{Vector{Float64}}, prob_fun::Vector{Function}, range_vec::Vector{Vector{Float64}}, vals2coeffs::Matrix{Float64}, coeffs2vals::Matrix{Float64}, weights_dict::Dict{String,Vector{Float64}}, order_comb::Vector{Vector{Int}}, integrals_vector::Vector{Float64}, sample_vars::Dict{String,Vector{ComplexF64}}, sample_num::Int, var_strs::Vector{String}, order_to_exp_std_dict::Union{Dict{Vector{Int},Tuple{Matrix{Float64}, Vector{Float64}, Matrix{Float64}}},Bool}=false)
        # check if bool?
        has_std_dict::Bool = true
        if typeof(order_to_exp_std_dict) == Bool
            has_std_dict = false
            order_to_exp_std_dict = Dict{Vector{Int},Tuple{Matrix{Float64}, Vector{Float64}, Matrix{Float64}}}()
        end
        new(locations, prob_fun, range_vec, vals2coeffs, coeffs2vals, weights_dict, order_comb, integrals_vector, sample_vars, sample_num, var_strs, has_std_dict, order_to_exp_std_dict)
    end
end

function Ex_and_Stdx(vals::Vector{T}, samples::SamplesAndWeights, operator_orders::Vector{Int}; get_cache::Bool=false) where T
    if sum(operator_orders) == 0
        if get_cache
            return reshape([1.0], 1, 1), [1.0], reshape([0.0], 1, 1)
        end
        return vals[1], 0.0 + 0.0im
    end
    if haskey(samples.order_to_exp_std_dict, operator_orders)
        vals2coeffs, E_weights, G_ij = samples.order_to_exp_std_dict[operator_orders]
        if get_cache
            return vals2coeffs, E_weights, G_ij
        end
    else
        error("Operator orders [$(operator_orders)] not found in samples.order_to_exp_std_dict.")
    end
    coeffs::Vector{ComplexF64} = vals2coeffs * vals
    Ex::ComplexF64 = E_weights' * vals
    Varx::ComplexF64 = coeffs' * G_ij * coeffs
    return Ex, sqrt(Varx)
end

function all_operator_order_comb_weights(max_op_order::Int, fun_vec::Vector{Function}, order_comb::Vector{Vector{Int}}, locations::Vector{Vector{Float64}}, range_vec::Vector{Vector{Float64}}; reltol::Float64=1e-6, abstol::Float64=1e-10, eps::Float64=1e-12, max_cum_order::Int=-1, loadsave::Bool=false, printing::Bool=false)::Dict{Vector{Int}, Tuple{Matrix{Float64}, Vector{Float64}, Matrix{Float64}}}
    # create order_combinations
    # max_op_order is the order of operators for which the weights std and vars are computed
    # fun_vec is probability distributions for each variable
    # order_comb is the order of polynomials of our integration
    # locations are the locations of the sample points, 
    # range_vec is the range of the variables
    # reltol and abstol are the tolerances for the quadrature
    # eps is the tolerance for removing orders
    # real_vals is a boolean that indicates if the values are real or complex
    # max_cum_order is the maximum polynomial order for which the computation is done (special cases: -1=max order of order_comb, 0=max order of order_comb*operator order)

    order_to_exp_std_dict::Dict{Vector{Int}, Tuple{Matrix{Float64}, Vector{Float64}, Matrix{Float64}}} = Dict() # Empty dict for the weights by order combination
    # Precompute integrals
    more_order_comb::Vector{Vector{Int}} = Vector{Vector{Int}}()
    n_comb = length(order_comb)
    for i in 1:n_comb
        for j in i:n_comb
            new_comb::Vector{Int} = order_comb[i] + order_comb[j]
            if !in(new_comb, more_order_comb)
                push!(more_order_comb, new_comb)
            end
        end
    end

    all_integrals::Vector{Vector{Float64}} = basis_function_integrals(fun_vec, range_vec, more_order_comb, reltol=reltol, abstol=abstol)
    normalized_all_integrals = [all_integrals[i] ./ all_integrals[i][1] for i in 1:length(all_integrals)]
    all_integrals_mat::Vector{Matrix{Float64}} = basis_function_integrals_mat(fun_vec, range_vec, more_order_comb, reltol=reltol, abstol=abstol)
    normalized_all_integrals_mat = [all_integrals_mat[i] ./ all_integrals_mat[i][1, 1] for i in 1:length(all_integrals)]
    # how many in total: length of all_vec_sum_n(max_op_order, 3)
    all_vecs::Vector{Vector{Int}} = all_vec_sum_n(max_op_order, 3)
    all_vecs_unique::Vector{Vector{Int}} = unique(sort.(all_vecs))
    how_many::Int = length(all_vecs)
    how_many_unique::Int = length(all_vecs_unique)
    values_list::Vector{Tuple{Matrix{Float64}, Vector{Float64}, Matrix{Float64}}} = []
    for counter in 1:how_many_unique
        operator_orders = all_vecs_unique[counter]
        sum_operator_orders = sum(operator_orders)
        if sum_operator_orders > 0
            all_locations, all_orders, all_range_vec = reduced_tensor_grid(operator_orders, fun_vec, order_comb, locations, range_vec; loc_ind=true, max_cum_order=max_cum_order, reltol=reltol, abstol=abstol)
            all_orders_integrals = repeat(normalized_all_integrals, sum_operator_orders)
            all_orders_integrals_mat = repeat(normalized_all_integrals_mat, sum_operator_orders)
            push!(values_list, Ex_and_Varx(all_orders_integrals, all_orders_integrals_mat, all_locations, all_orders, all_range_vec; reltol=reltol, abstol=abstol, eps=eps, loadsave=loadsave, printing=printing))
            if !isapprox(sum(values_list[end][2]), 1.0, atol=1e-8)
                println("Warning: sum of weights for operator_orders $operator_orders is not close to 1.0, sum-1=", 1.0-sum(values_list[end][2]))
            end
        end
    end
    # check if sum of weights for every values list entry are close to 1
    for (counter, operator_orders) in enumerate(all_vecs)
        #sum_operator_orders = sum(operator_orders)#
        if sum(operator_orders) > 0
            # find the index of the unique operator_orders
            sorted_operator_orders = sort(operator_orders)
            # find the index of the unique operator_orders which is equal to sorted_operator_orders
            index = findfirst(x -> x == sorted_operator_orders, all_vecs_unique)
            order_to_exp_std_dict[operator_orders] = values_list[index-1]
        end
    end
    return order_to_exp_std_dict
end



#### Tensor Grids --- For Fidelities ####
function reduced_tensor_grid(operator_order::Int, fun_vec::Vector{Function}, order_comb::Vector{Vector{Int}}, locations::Vector{Vector{Float64}}, range_vec::Vector{Vector{Float64}}; max_cum_order::Int=0, reltol::Float64=1e-6, abstol::Float64=1e-10, loc_ind::Bool=true)
    # For higher orders such as multi qubit fidelities, we need to generate a tensor grid of locations, but mgreater_range based
    dir_num = length(fun_vec)
    orders_vec = [maximum([order_comb[i][j] for i in 1:length(order_comb)]) for j in 1:dir_num]
    if max_cum_order == 0
        max_cum_order = maximum(orders_vec) * operator_order
    elseif max_cum_order == -1
        max_cum_order = maximum(orders_vec)
    end
    order_comb_len = length(order_comb)
    all_orders::Vector{Vector{Int}} = []
    curr_orders::Vector{Int} = []
    for i_s in mgreater_range(order_comb_len, operator_order)
        curr_orders = []
        for i in i_s
            append!(curr_orders, order_comb[i])
        end
        sum_curr_order = sum(curr_orders)
        if sum_curr_order <= max_cum_order
            push!(all_orders, curr_orders)
        end
    end
    # all_integrals
    all_integrals::Vector{Vector{Float64}} = basis_function_integrals(fun_vec, range_vec, order_comb, reltol=reltol, abstol=abstol)
    normalized_all_integrals = [all_integrals[i] ./ all_integrals[i][1] for i in 1:length(all_integrals)] # normalizes the integrals to 1, so that int dx p(x) = 1
    all_orders_integrals = repeat(normalized_all_integrals, operator_order)
    # Create all combinations of locations for all directions 
    all_locations::Vector{Vector{Float64}} = []
    len_locations = length(locations)
    curr_locations::Vector{Float64} = []
    for i_s in mgreater_range(len_locations, operator_order; offdiag=0)
        curr_locations = []
        for i in i_s
            append!(curr_locations, locations[i])
        end
        push!(all_locations, curr_locations)
    end
    locations_ind::Vector{Vector{Int}} = []
    if loc_ind
        for i_s in mgreater_range(len_locations, operator_order; offdiag=0)
            push!(locations_ind, i_s)
        end
    end
    all_range_vec = repeat(range_vec, operator_order)
    if !loc_ind
        return all_locations, all_orders, all_orders_integrals, all_range_vec
    else
        return all_locations, all_orders, all_orders_integrals, all_range_vec, locations_ind
    end
end

# Generalize reduced_tensor_grid script to include different orders for x y and z operators
# for each we would use a reduced tensor grid of the fidelity type 
function reduced_tensor_grid(operator_orders::Vector{Int}, fun_vec::Vector{Function}, order_comb::Vector{Vector{Int}}, locations::Vector{Vector{Float64}}, range_vec::Vector{Vector{Float64}}; max_cum_order::Int=0, reltol::Float64=1e-6, abstol::Float64=1e-10, loc_ind::Bool=true, do_int::Bool=false)
    dir_num = length(fun_vec)
    sum_operator_order = sum(operator_orders)
    orders_vec = [maximum([order_comb[i][j] for i in 1:length(order_comb)]) for j in 1:dir_num]
    if max_cum_order == 0
        max_cum_order = maximum(orders_vec) * sum_operator_order
    elseif max_cum_order == -1
        max_cum_order = maximum(orders_vec)
    end
    order_comb = sort(order_comb, by=x -> sum(x))
    order_comb_len = length(order_comb)
    all_orders::Vector{Vector{Int}} = []
    curr_orders::Vector{Int} = []
    for i_s in multi_mgreater_range(order_comb_len, operator_orders)
        curr_orders = []
        for i_d in i_s
            for i in i_d
                append!(curr_orders, order_comb[i])
            end
        end
        sum_curr_order = sum(curr_orders)
        if sum_curr_order <= max_cum_order
            push!(all_orders, curr_orders)
        end
    end
    # Create all combinations of locations for all directions 
    all_locations::Vector{Vector{Float64}} = []
    len_locations = length(locations)
    curr_locations::Vector{Float64} = []
    for i_s in multi_mgreater_range(len_locations, operator_orders; offdiag=0)
        curr_locations = []
        for i_d in i_s
            for i in i_d
                append!(curr_locations, locations[i])
            end
        end
        push!(all_locations, curr_locations)
    end
    locations_ind::Vector{Vector{Vector{Int}}} = []
    if loc_ind
        for i_s in multi_mgreater_range(len_locations, operator_orders; offdiag=0)
            push!(locations_ind, i_s)
        end
    end
    # remove elements with too high orders 

    all_range_vec = repeat(range_vec, sum_operator_order)

    if do_int
        all_integrals::Vector{Vector{Float64}} = basis_function_integrals(fun_vec, range_vec, order_comb, reltol=reltol, abstol=abstol)
        normalized_all_integrals = [all_integrals[i] ./ all_integrals[i][1] for i in 1:length(all_integrals)]
        all_orders_integrals = repeat(normalized_all_integrals, sum_operator_order)
        if !loc_ind
            return all_locations, all_orders, all_orders_integrals, all_range_vec
        else
            return all_locations, all_orders, all_orders_integrals, all_range_vec, locations_ind
        end
    else
        if !loc_ind
            return all_locations, all_orders, all_range_vec
        else
            return all_locations, all_orders, all_range_vec, locations_ind
        end
    end
end

#### Single Script to make samples and weights
function prepare_samples_and_weights(prob_fun::Vector{Function}, dim_names::Vector{String}, which_dims::Vector{String}, N::Int, range_vec::Vector{Vector{Float64}}, orders_vec::Union{Vector{Int},Int}=Vector{Int}(); prepare_analysis::Int=-1, extra_points::Union{Int,Vector{Int}}=0, on_borders::Bool=false, type::Symbol=:equidistant, max_comb_order::Int=0, loadsave::Bool=false, printing::Bool=false, reltol::Float64=1e-15, abstol::Float64=1e-15, rng::MersenneTwister=MersenneTwister(1234), weighted::Bool=false)
    # --- Samples and Weights --------------------------------------------------------------------
    # prob_fun_vec is the probabilitiy distribution of the samples for different dimensions of the system
    # dim_names is the names of the dimensions (andd needs to be the same length as prob_fun_vec)
    # which_dims is the dimensions for which the weights are to be prepared
    # N is the number of spins in the system
    # range_vec is the range of the samples for each dimension [in which interval to sample the grid]
    # order_vec is the order of the samples [to which order are the basis functions expanded]
    # Additional parameters:
    # prepare_analysis is an integer that indicates if and to which order (Int) of analysis to be prepared 
    # on_borders is a boolean to sample the grid on the borders
    # extra_points is the number of extra points to be added to the sample grid
    # type is a Symbol (:equidistant, :equiprobable, :random, :randomprob)
    # equidistant is a boolean to sample the grid equidistantly or (if false) probabilistically
    # max_comb_order is the maximum order of the combined basis functions in the analytis functions order_comb[i]+order_comb[j]
    #       max_comb_order = 0 ==> max_comb_order = maximum(orders_vec) * operator_order
    #       max_comb_order = -1 ==> max_comb_order = maximum(orders_vec)
    # weighted is a boolean to weight the samples by the probability distribution for a weighted fit (or not)
    dim_num = length(prob_fun)
    if dim_num != length(range_vec)
        error("prob_fun and range_vec must have the same length")
    end
    locations, order_comb = cubed_sample_grid(prob_fun, range_vec, orders_vec; on_borders=on_borders, extra_points=extra_points, reltol=reltol, abstol=abstol, type=type)
    sample_vars::Dict{String,Vector{ComplexF64}} = Dict()
    for (i, dim) in enumerate(dim_names)
        sample_vars[dim] = [locations[j][i] for j in 1:length(locations)]
    end
    weights_dict::Dict{String,Vector{Float64}} = Dict()
    # if empty string isn'T in which_dims, then add it! 
    if !any([x == "" for x in which_dims])
        push!(which_dims, "")
    end
    vals2coeffs::Matrix{Float64} = zeros(length(order_comb), length(locations))
    coeffs2vals::Matrix{Float64} = zeros(length(locations), length(order_comb))
    integrals_vector::Vector{Float64} = zeros(length(order_comb))
    for dim in which_dims
        g_fun_vec::Vector{Function} = copy(prob_fun)
        if length(dim) > 0
            if !(dim in dim_names)
                error("Dimension ", dim, " not found in dim_names")
            end
            which_ind = findfirst(x -> x == dim, dim_names)
            g_fun_vec[which_ind] = x -> x * N * prob_fun[which_ind](x)
        end
        if dim == ""
            weights, vals2coeffs, coeffs2vals, integrals_vector = weights_and_errors(g_fun_vec, locations, order_comb, range_vec, fun_return=true, reltol=reltol, weighted=weighted)[[1, 3, 4, 5]]
        else
            weights = weights_and_errors(g_fun_vec, locations, order_comb, range_vec, fun_return=true, reltol=reltol)[1]
        end
        weights_dict[dim] = weights
    end
    sample_num = length(locations)
    if prepare_analysis > 0
        order_to_exp_std_dict = all_operator_order_comb_weights(prepare_analysis, prob_fun, order_comb, locations, range_vec, max_cum_order=max_comb_order, loadsave=loadsave, printing=printing, reltol=reltol, abstol=abstol)
        return SamplesAndWeights(locations, prob_fun, range_vec, vals2coeffs, coeffs2vals, weights_dict, order_comb, integrals_vector, sample_vars, sample_num, dim_names, order_to_exp_std_dict)
    else
        return SamplesAndWeights(locations, prob_fun, range_vec, vals2coeffs, coeffs2vals, weights_dict, order_comb, integrals_vector, sample_vars, sample_num, dim_names)
    end
end



function get_weights(samples::SamplesAndWeights, spin_orders::Union{Int,Vector{Int}}=1; normalize::Bool=false)
    # if spin_orders is integer 
    if typeof(spin_orders) == Int
        error("spin_orders must be a vector of integers")
        #if normalize
        #    return samples.weights_dict[""] ./ maximum(samples.weights_dict[""])
        #end
        #return samples.weights_dict[""]
    end
    # check if spin orders is in 
    if samples.has_std_dict
        # check for key 
        if spin_orders == [0, 0, 0]
            return [1.0]
        end
        if haskey(samples.order_to_exp_std_dict, spin_orders)
            weights = samples.order_to_exp_std_dict[spin_orders][2]
            if normalize
                return weights ./ maximum(weights)
            end
            return weights
        end
        println("Spin orders ", spin_orders, " not found in order_to_exp_std_dict, construct weights from weight products as an approximation.")
    end
    weights = samples.weights_dict[""]
    if sum(spin_orders) == 1
        if normalize
            return weights ./ maximum(weights)
        end
        return weights
    end
    # more than 1 spin operator
    all_weights = zeros(samples.sample_num)
    for (i, combs) in enumerate(multi_mgreater_range_flat(samples.sample_num, curr_spin_orders))
        all_weights[i] = prod([weights[c] for c in combs])
    end
    weights = all_weights
    if normalize
        return weights ./ maximum(weights)
    end
    return weights
end