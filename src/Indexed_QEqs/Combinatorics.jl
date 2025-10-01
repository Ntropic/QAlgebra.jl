module  QCombinatorics
export generate_ensembles, generate_multi_ensembles, ensemble_iterator, multi_ensemble_iterator
# ==============================> Index Generation <=================================================================================
# --- combinations helper ---
function combinations(v::Vector{Int}, k::Int)
    result = Vector{Vector{Int}}()
    n = length(v)

    function backtrack(start::Int, chosen::Vector{Int})
        if length(chosen) == k
            push!(result, copy(chosen))
            return
        end
        for i in start:n
            push!(chosen, v[i])
            backtrack(i+1, chosen)
            pop!(chosen)
        end
    end

    backtrack(1, Int[])
    return result
end

# --- ensemble generator ---
function generate_ensembles(block_sizes::Vector{Int}, max_n::Int)
    result = Vector{Vector{Vector{Int}}}()

    function backtrack(block_idx::Int, used::Vector{Int}, curr::Vector{Vector{Int}})
        if block_idx > length(block_sizes)
            push!(result, deepcopy(curr))
            return
        end

        bsize = block_sizes[block_idx]
        available = setdiff(1:max_n, used)
        for combo in combinations(collect(available), bsize)
            push!(curr, combo)
            append!(used, combo)
            backtrack(block_idx+1, used, curr)
            deleteat!(used, length(used)-bsize+1:length(used))
            pop!(curr)
        end
    end

    backtrack(1, Int[], Vector{Vector{Int}}())
    return result
end
function generate_multi_ensembles(block_sizes_list::Vector{Vector{Int}}, max_ns::Vector{Int})
    n_ensembles = length(block_sizes_list)
    @assert length(max_ns) == n_ensembles

    # precompute all ensembles for each case
    all_ensembles = [generate_ensembles(block_sizes_list[i], max_ns[i]) for i in 1:n_ensembles]

    result = Vector{Vector{Vector{Vector{Int}}}}()

    function backtrack(idx::Int, curr::Vector{Vector{Vector{Int}}})
        if idx > n_ensembles
            push!(result, deepcopy(curr))
            return
        end
        for ens in all_ensembles[idx]
            push!(curr, ens)
            backtrack(idx+1, curr)
            pop!(curr)
        end
    end

    backtrack(1, Vector{Vector{Vector{Int}}}())
    return result
end

# =========> Index Generation as Iterator

"""
    EnsembleIterator(block_sizes::Vector{Int}, max_n::Int)

Lazy iterator over all ensembles of given block sizes from 1:max_n.
Yields `Vector{Vector{Int}}` like the old `generate_ensembles`.
"""
struct EnsembleIterator
    block_sizes::Vector{Int}
    max_n::Int
end

Base.IteratorSize(::Type{EnsembleIterator}) = Base.SizeUnknown()
Base.eltype(::Type{EnsembleIterator}) = Vector{Vector{Int}}

function Base.iterate(it::EnsembleIterator, state=(1, Int[], Vector{Vector{Int}}()))
    block_sizes, max_n = it.block_sizes, it.max_n
    block_idx, used, curr = state

    # backtracking state machine
    while true
        if block_idx > length(block_sizes)
            return deepcopy(curr), (block_idx, copy(used), copy(curr))  # yield
        end

        bsize = block_sizes[block_idx]
        available = setdiff(1:max_n, used)

        # if no candidates left, stop
        isempty(available) && return nothing

        # try combinations lazily
        for combo in Iterators.product(combinations(collect(available), bsize))
            new_used = vcat(used, combo...)
            new_curr = [curr...; combo]
            return deepcopy(new_curr), (block_idx+1, new_used, new_curr)
        end

        return nothing
    end
end

function ensemble_iterator(block_sizes::Vector{Int}, max_n::Int)
    Channel() do ch
        function backtrack(block_idx::Int, used::Vector{Int}, curr::Vector{Vector{Int}})
            if block_idx > length(block_sizes)
                put!(ch, deepcopy(curr))
                return
            end
            bsize = block_sizes[block_idx]
            available = setdiff(1:max_n, used)
            for combo in combinations(collect(available), bsize)
                push!(curr, combo)
                append!(used, combo)
                backtrack(block_idx+1, used, curr)
                deleteat!(used, length(used)-bsize+1:length(used))
                pop!(curr)
            end
        end
        backtrack(1, Int[], Vector{Vector{Int}}())
    end
end

"""
    multi_ensemble_iterator(block_sizes_list::Vector{Vector{Int}}, max_ns::Vector{Int})

Lazy iterator over ensembles for each (block_sizes, max_n) pair.
Yields `Vector{Vector{Vector{Int}}}` like `generate_multi_ensembles`.
"""
function multi_ensemble_iterator(block_sizes_list::Vector{Vector{Int}}, max_ns::Vector{Int})
    n_ensembles = length(block_sizes_list)
    @assert length(max_ns) == n_ensembles

    Channel() do ch
        function backtrack(idx::Int, curr::Vector{Vector{Vector{Int}}})
            if idx > n_ensembles
                put!(ch, deepcopy(curr))
                return
            end
            for ens in ensemble_iterator(block_sizes_list[idx], max_ns[idx])
                push!(curr, ens)
                backtrack(idx+1, curr)
                pop!(curr)
            end
        end
        backtrack(1, Vector{Vector{Vector{Int}}}())
    end
end
    
end