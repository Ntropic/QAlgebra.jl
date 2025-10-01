using Combinatorics

# Function to calculate the sum of binomials for each S and R
function calculateBinomialSums(n::Int, m::Int)::Array{Int,3}
    # Initialize the result matrix
    bin_array::Array{Int,3} = zeros(Int, n, n, m - 1)
    bin_cache::Array{Int,2} = zeros(Int, n, m)
    for R in 1:m
        for ind in 1:n
            bin_cache[ind, R] = binomial(ind + R - 1, R - 1)
        end
    end
    # Iterate over each S and R to calculate the sum of binomials
    for R in 2:m                  # curr position in vector
        for j in 0:n-1       #prev_digit
            curr_n = n - j
            for S in 0:(curr_n-1)        # curr digit
                # Store the result in the matrix
                bin_array[S+1+j, j+1, m-R+1] = sum([bin_cache[curr_n-s, R] for s in 1:S])
            end
        end
    end
    return bin_array
end
## Example usage
#n = 3
#m = 3
#bin_array = calculateBinomialSums(n, m)

function determine_indexes_gen(n::Int, m::Int)::Function
    # Generate a function to determine the indexes of a sequence vector::Vector{Int}
    # n is the maximum digit value and m the maximum length of the sequences
    # Initialize the result matrix
    bin_arr::Array{Int,3} = calculateBinomialSums(n, m)
    function determine_indexes(vector::Vector{Int})::Int
        len = length(vector)
        if len == 0
            return 1
        elseif len == 1
            return vector[1]
        end
        sequence_shift = m - len
        index::Int = bin_arr[vector[1], 1, sequence_shift+1]
        for i in 2:len-1
            index += bin_arr[vector[i], vector[i-1], i+sequence_shift]
        end
        max_value::Int = bin_arr[n, 1, sequence_shift+1] + 1
        return index + vector[len] - vector[len-1] + 1
    end
    return determine_indexes
end
## Test 
#n = 1000
#m = 3
#determine_indexes = determine_indexes_gen(n, m)
#determine_indexes([90, 91, 92])

function determine_multi_indexes_gen_one_out(n::Int, m::Int)::Function
    # Generate a function to determine the indexes of a sequence vector::Vector{Int}
    # n is the maximum digit value and m the maximum length of the sequences
    # Initialize the result matrix
    bin_arr::Array{Int,3} = calculateBinomialSums(n, m)
    function determine_indexes_and_max_val(vector::Vector{Int})::Tuple{Int,Int}
        len = length(vector)
        if len == 0
            return 1, 1
        elseif len == 1
            return vector[1], n
        end
        sequence_shift = m - len
        index::Int = bin_arr[vector[1], 1, sequence_shift+1]
        for i in 2:len-1
            index += bin_arr[vector[i], vector[i-1], i+sequence_shift]
        end
        max_value::Int = bin_arr[n, 1, sequence_shift+1] + 1
        return index + vector[len] - vector[len-1] + 1, max_value
    end
    function determine_combined_indexes(vector::Vector{Vector{Int}})::Int
        local_indexes::Vector{Int} = zeros(Int, 3)
        max_values::Vector{Int} = zeros(Int, 3)
        for i in 1:3
            local_indexes[i], max_values[i] = determine_indexes_and_max_val(vector[i])
        end
        #println(local_indexes, max_values)
        return (local_indexes[1] - 1) * max_values[2] * max_values[3] + (local_indexes[2] - 1) * max_values[3] + local_indexes[3]
    end
    return determine_combined_indexes
end
## Test 
#n = 5
#m = 3
#determine_indexes = determine_multi_indexes_gen_one_out(n, m)
#inds::Vector{Vector{Int}} = [[], [1], [1, 4]]
#determine_indexes(inds)

function determine_multi_indexes_gen(n::Int, m::Int)::Tuple{Function,Function}
    # Generate a function to determine the indexes of a sequence vector::Vector{Int}
    # n is the maximum digit value and m the maximum length of the sequences
    # Initialize the result matrix
    bin_arr::Array{Int,3} = calculateBinomialSums(n, m)
    max_values::Vector{Int} = [1, n]
    append!(max_values, [bin_arr[n, 1, m-len+1] + 1 for len in 2:m])
    function determine_indexes(vector::Vector{Int})::Int
        len = length(vector)
        if len == 0
            return 1
        elseif len == 1
            return vector[1]
        end
        sequence_shift = m - len
        index::Int = bin_arr[vector[1], 1, sequence_shift+1]
        for i in 2:len-1
            index += bin_arr[vector[i], vector[i-1], i+sequence_shift]
        end

        return index + vector[len] - vector[len-1] + 1
    end
    function determine_combined_indexes(vector::Vector{Vector{Int}})::Int
        local_indexes::Vector{Int} = zeros(Int, 3)
        for i in 1:3
            local_indexes[i] = determine_indexes(vector[i])
        end
        max_val2::Int = max_values[length(vector[2])+1]
        max_val3::Int = max_values[length(vector[3])+1]
        #println(local_indexes, max_values)
        return (local_indexes[1] - 1) * max_val2 * max_val3 + (local_indexes[2] - 1) * max_val3 + local_indexes[3]
    end
    function determine_combined_indexes_with_zero(vector::Vector{Vector{Int}}, where_zero::Int)::Vector{Int}
        lens::Vector{Int} = [length(v) for v in vector]
        local_indexes::Vector{Int} = [1, 1, 0]
        for i in 1:3
            if i != where_zero
                local_indexes[i] = determine_indexes(vector[i])
            else
                lens[i] += 1
            end
        end
        zero_indexes::Vector{Int} = Vector{Int}(undef, n)
        curr_var::Vector{Int} = Vector{Int}(undef, length(vector[where_zero]) + 1)
        for i in 1:n
            curr_var = [i, vector[where_zero]...]
            zero_indexes[i] = determine_indexes(sort!(curr_var))
        end
        max_val2::Int = max_values[lens[2]+1]
        max_val3::Int = max_values[lens[3]+1]
        #println(local_indexes, max_values)
        prefix::Int = (local_indexes[1] - 1) * max_val2 * max_val3 + (local_indexes[2] - 1) * max_val3 + local_indexes[3]
        pref::Int = max_val3
        if where_zero == 1
            pref *= max_val2
            return prefix .+ zero_indexes .* pref .- pref

        elseif where_zero == 2
            return prefix .+ zero_indexes .* pref .- pref
        else
            return prefix .+ zero_indexes
        end
    end
    return determine_combined_indexes, determine_combined_indexes_with_zero
end
#n = 5
#m = 3
#determine_indexes, determine_combined_indexes_with_zero = determine_multi_indexes_gen(n, m)
## Test combined_with_zero generation
#initial_indexes = Vector{Vector{Int}}([[5], [1], []])
#iterative_index_list::Vector{Int} = []
#for (i, xyz) in enumerate(multi_mgreater_range(n, [1,1,1]))
#    # check if the index xyz contains initial_index[2] elements 
#    if issubset(initial_indexes[1], xyz[1]) && issubset(initial_indexes[2], xyz[2]) && issubset(initial_indexes[3], xyz[3])
#        push!(iterative_index_list, i)
#    end
#end
#method_derived_indexes::Vector{Int} = determine_combined_indexes_with_zero(initial_indexes, 3)
#for i in 1:n 
#    println(method_derived_indexes[i], " ", iterative_index_list[i])
#end