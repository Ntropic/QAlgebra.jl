using QuadGK
using Random
using Roots
using Base.Threads

macro usethreads(multithreaded, expr::Expr)
    ex = quote
        if $multithreaded
            Threads.@threads $expr
        else
            $expr
        end
    end
    esc(ex)
end


#############################################################################################################
####### Distribution Sampling ###############################################################################
# A function to sample values from distributions, taking a function handle as input, normalizing the distribution given by it, via integration over bounds (specified via input)
# then pick random numbers between 0 and 1 and use root finding to find the corresponding point in the distribution.
# add an input argument to specify that function has been normalized
function sample_from_distribution(fun::Function, bounds::Vector{Float64}, n::Int; rng::MersenneTwister=MersenneTwister(1234), normalized::Bool=false, reltol::Float64=1e-15, abstol::Float64=1e-15, reorder::Bool=false, threaded::Bool=true, presort::Bool=true)::Tuple{Vector{Float64},Float64}
    # normalize the function
    norm::Float64 = 1.0
    if !normalized
        norm = quadgk(fun, bounds[1], bounds[2]; rtol=reltol)[1]
    end
    #println("norm = ", norm)
    # sample n random numbers between 0 and 1
    random_numbers::Vector{Float64} = rand(rng, n)
    # sort the random numbers
    if presort
        sort!(random_numbers)
    end
    cdf::Function = x -> quadgk(x -> fun(x) / norm, bounds[1], x; rtol=reltol)[1]
    # find the root of cdf(x) - random_number = 0
    values::Vector{Float64} = zeros(Float64, n)
    center_pos = find_zero(x -> cdf(x) - 1 / 2, (bounds[1], bounds[2]))
    @usethreads threaded for i in 1:n
        values[i] = find_zero(x -> cdf(x) - random_numbers[i], (bounds[1], bounds[2]))
    end
    if reorder
        values = sort(values, by=x -> abs(x - center_pos))
    end
    return values, center_pos
end
# Vectorized variant 
function sample_from_distributions(funs::Vector{Function}, bounds::Vector{Vector{Float64}}, n::Int; rng::MersenneTwister=MersenneTwister(1234), normalized::Bool=false, reltol::Float64=1e-15, abstol::Float64=1e-15, reorder::Bool=false, threaded::Bool=true)::Tuple{Vector{Vector{Float64}},Vector{Float64}}
    values::Vector{Vector{Float64}} = []
    center_pos::Vector{Float64} = []
    for (fun, bound) in zip(funs, bounds)
        curr_values, curr_center_pos = sample_from_distribution(fun, bound, n; rng=rng, normalized=normalized, reltol=reltol, abstol=abstol, reorder=false, threaded=threaded, presort=false)
        push!(values, curr_values)
        push!(center_pos, curr_center_pos)
    end
    return values, center_pos
end

function sample_from_interval(bounds::Vector{Float64}, n::Int; rng::MersenneTwister=MersenneTwister(1234))
    return rand(rng, n) .* (bounds[2] - bounds[1]) .+ bounds[1]
end
function sample_from_intervals(bounds::Vector{Vector{Float64}}, n::Int; rng::MersenneTwister=MersenneTwister(1234))
    return [rand(rng, n) .* (b[2] - b[1]) .+ b[1] for b in bounds]
end

# not a random sample, but equidistant in the distribution
function equidistant_sample_from_distribution(fun::Function, bounds::Vector{Float64}, n::Int; on_borders::Bool=false, normalized::Bool=false, reltol::Float64=1e-6, reorder::Bool=false, border_dist::Float64=1.0)::Tuple{Vector{Float64},Float64}
    # samples equidistant with respect to the cdf height
    # if reorder, then reorder the samples by distace from the center
    # normalize the function
    norm::Float64 = 1.0
    if !normalized
        norm = quadgk(fun, bounds[1], bounds[2]; rtol=reltol)[1]
    end
    #println("norm = ", norm)
    # sample n random numbers between 0 and 1
    if on_borders
        nums = range(0.0, 1.0, length=n)
    else
        delta = 1.0 / n * border_dist
        nums = range(delta / 2, 1.0 - delta / 2, length=n)
    end
    cdf::Function = x -> quadgk(x -> fun(x) / norm, bounds[1], x; rtol=reltol)[1]
    # find the root of cdf(x) - random_number = 0
    values::Vector{Float64} = zeros(Float64, n)
    center_pos = find_zero(x -> cdf(x) - 1 / 2, (bounds[1], bounds[2]))
    if !on_borders
        for i in 1:n
            values[i] = find_zero(x -> cdf(x) - nums[i], (bounds[1], bounds[2]))
        end
    else
        for i in 2:n-1
            values[i] = find_zero(x -> cdf(x) - nums[i], (bounds[1], bounds[2]))
        end
        values[1] = bounds[1]
        values[n] = bounds[2]
    end
    if reorder
        values = sort(values, by=x -> abs(x - center_pos))
    end
    return values, center_pos
end
## Test
#fun = x -> exp(-abs(x))
#bounds = [-1.0, 1.0]
#n = 10000
#values, center_pos = sample_from_distribution(fun, bounds, 10; reorder = true)
#plt = histogram(values, bins=100, label="Sampled Values", normalize=true)

function equidistant_sample(bounds::Vector{Float64}, n::Int; on_borders::Bool=false, reorder::Bool=false, border_dist::Float64=1.0)::Tuple{Vector{Float64},Float64}
    # samples equidistant in interval [bounds[1], bounds[2]], n equidistant points, if on_borders, then include the borders, otherwise not on borders
    # if reorder, then reorder the samples by distace from the center
    center_pos = (bounds[1] + bounds[2]) / 2
    if on_borders
        values = collect(range(bounds[1], bounds[2], length=n))
    else
        delta = (bounds[2] - bounds[1]) / (n + 1) * border_dist     # full border distance
        #delta = (bounds[2] - bounds[1]) / n * border_dist          # half border distance
        values = collect(range(bounds[1] + delta / 2, bounds[2] - delta / 2, length=n))
    end
    if reorder
        values = sort(values, by=x -> abs(x - center_pos))
    end
    return values, center_pos
end

# Sample to accuracy
function sample_gaussians_to_accuracy(fun::Function, x0s::Vector{Float64}, sigmas::Vector{Float64}, reference_result::Float64, max_error::Float64; max_num::Int=-1, rng::MersenneTwister=MersenneTwister(1234), return_samples::Bool=true)
    # sample points from a product of gaussian distributions with means x0s[i] and sigmas sigmas[i] until 
    # the difference between the weighted sum of fun values at the sampled points and the reference_result is smaller than max_error
    sum_of_values::Float64 = 0.0
    number_of_values::Int = 0
    curr_value::Float64 = -10^100 # sum_of_values / number_of_values
    dims::Int = length(x0s)
    x::Vector{Float64} = zeros(Float64, dims)
    y::Float64 = 0.0
    if return_samples
        x_s::Vector{Vector{Float64}} = []
        y_s::Vector{Float64} = []
        vals::Vector{Float64} = []
    end
    while number_of_values != max_num && abs(curr_value - reference_result) > max_error
        # sample a point
        x = randn(rng, dims) .* sigmas .+ x0s
        # evaluate the function at the point
        y = fun(x)
        # update the sum and number of values
        sum_of_values += y
        number_of_values += 1
        curr_value = sum_of_values / number_of_values
        if return_samples
            push!(x_s, x)
            push!(y_s, y)
            push!(vals, curr_value)
        end
    end
    if !return_samples
        return number_of_values, curr_value
    else
        return number_of_values, curr_value, x_s, y_s, vals
    end
end

function samples_4_samples_of_samples(fun::Function, x0s::Vector{Float64}, sigmas::Vector{Float64}, how_many::Int; rng=MersenneTwister(1234), bounds::Union{Vector{Vector{Float64}},Bool}=false)
    dim::Int = length(x0s)
    x_s::Vector{Vector{Float64}} = []
    if typeof(bounds) == Bool
        if bounds == false
            bounds = [[-Inf, Inf] for i in 1:dim]
        end
    end
    for i in 1:how_many
        repeat = true
        while repeat
            curr_x = randn(rng, dim) .* sigmas .+ x0s
            if all([curr_x[i] > bounds[i][1] && curr_x[i] < bounds[i][2] for i in 1:dim])
                push!(x_s, curr_x)
                repeat = false
            end
        end
    end
    y_s::Vector{Float64} = zeros(Float64, how_many)
    if dim == 1
        for i in 1:how_many
            y_s[i] = fun(x_s[i][1])
        end
    else
        for i in 1:how_many
            y_s[i] = fun(x_s[i])
        end
    end
    return x_s, y_s
end
function samples_of_samples(x_s::Vector{Vector{Float64}}, y_s::Vector{Float64}, how_many_subsets::Int, how_many_samples_each::Int; rng=MersenneTwister(1234))
    # make how many subsets of how_many_samples_each samples from the samples in x_s and y_s and construct the sum of values/how_many_samples_each
    n = length(x_s)
    if how_many_samples_each > n
        error("how_many_samples_each > n")
    end
    vals::Vector{Float64} = zeros(Float64, how_many_subsets)
    for i in 1:how_many_subsets
        indices = rand(rng, 1:n, how_many_samples_each)
        vals[i] = sum(y_s[indices]) / how_many_samples_each
    end
    return vals
end
function samples_of_samples(fun::Function, x0s::Vector{Float64}, sigmas::Vector{Float64}, how_many_subsets::Int, how_many_samples_each::Int; bounds::Union{Vector{Vector{Float64}},Bool}=false, rng=MersenneTwister(1234))
    # make how many subsets of how_many_samples_each samples from the samples in x_s and y_s and construct the sum of values/how_many_samples_each
    vals::Vector{Float64} = zeros(Float64, how_many_subsets)
    for i in 1:how_many_subsets
        _, y_s = samples_4_samples_of_samples(fun, x0s, sigmas, how_many_samples_each; rng=rng, bounds=bounds)
        vals[i] = sum(y_s) / how_many_samples_each
    end
    return vals
end

function weighted_sample_sum(fun::Function, pdf::Function, locations::Vector{Float64})::Float64
    # sample from pdf and evaluate fun at the sampled points, then sum up the values weighted by the pdf
    # locations is a vector of locations where to sample from the pdf
    # pdf is a function that takes a vector of locations and returns a vector of pdf values at those locations
    # fun is a function that takes a vector of locations and returns a vector of function values at those locations
    # returns the weighted sum of fun values
    pdf_vals = pdf(locations)
    fun_vals = fun(locations)
    return sum(pdf_vals .* fun_vals) / sum(pdf_vals)
end
function weighted_sample_sum(fun::Function, pdf::Vector{Function}, locations::Vector{Vector{Float64}})::Float64
    # sample from pdf and evaluate fun at the sampled points, then sum up the values weighted by the pdf
    # locations is a vector of locations where to sample from the pdf
    # pdf is a vector of functions that take a vector of locations and return a vector of pdf values at those locations
    # fun is a function that takes a vector of locations and returns a vector of function values at those locations
    # returns the weighted sum of fun values
    dim = length(pdf)
    n = length(locations)
    summer::Float64 = 0.0
    pdf_vals = zeros(Float64, dim)
    cum_prob = 0.0
    for i in 1:n
        for j in 1:dim
            pdf_vals[j] = pdf[j](locations[i][j])
        end
        if dim == 1
            summer += fun(locations[i][1]) * pdf_vals[1]
        else
            summer += fun(locations[i]) * prod(pdf_vals)
        end
        cum_prob += prod(pdf_vals)
    end
    return summer / cum_prob
end

function better_weighted_sample_sum(fun::Function, pdf::Function, locations::Vector{Float64}, borders::Vector{Float64})
    # use integral over intervals for probabilities
    min_val = borders[1]
    max_val = borders[2]
    interval_borders = [min_val]
    for i in 1:length(locations)-1
        push!(interval_borders, (locations[i] + locations[i+1]) / 2)
    end
    push!(interval_borders, max_val)
    # get the pdf weights
    interval_weights = zeros(length(interval_borders) - 1)
    for i in 1:length(interval_borders)-1
        interval_weights[i] = quadgk(pdf, interval_borders[i], interval_borders[i+1])[1]
    end
    # get the function values 
    fun_values = fun.(locations)
    # get the weighted sum
    weighted_sum = sum(fun_values .* interval_weights)
    return weighted_sum
end

function better_weighted_sample_sum(fun::Function, pdf::Vector{Function}, locations::Vector{Vector{Float64}}, borders::Vector{Float64})
    dim = length(pdf)
    if dim > 1
        error("dim > 1 is not implemented yet")
    end
    locs = [l[1] for l in locations]
    return better_weighted_sample_sum(fun, pdf[1], locs, borders)
end