module  QCombinatorics
export ensemble_iterator, ensemble_iterator_continuum, multi_ensemble_iterator, multi_ensemble_iterator_continuum

using ..Indexing: MultiEnsembleWorkspace

"""
# Combination iterators

Ensemble index combinations are generated lazily via `Channel`s. Typical usage:

```julia
mws = MultiEnsembleWorkspace([4, 3], BitVector([false, true]))
it = multi_ensemble_iterator([[2], [1, 1]], mws)
for combo in it
    @show combo
end
```

The iterator respects the continuum mask stored in the workspace and visits all
sorted block configurations without allocating intermediate arrays of
combinations.
"""

# ==============================> Helpers <=================================================================================

function _for_each_combination_discrete(available::Vector{Int}, k::Int, callback::Function)
    block = Vector{Int}(undef, k)
    last = length(available)
    function rec(pos::Int, start_idx::Int)
        if pos > k
            callback(block)
            return
        end
        remaining = k - pos + 1
        max_start = last - remaining + 1
        @inbounds for idx in start_idx:max_start
            block[pos] = available[idx]
            rec(pos + 1, idx + 1)
        end
    end
    k == 0 ? callback(Int[]) : rec(1, 1)
    return nothing
end

function _for_each_combination_discrete(callback::Function, available::Vector{Int}, k::Int)
    _for_each_combination_discrete(available, k, callback)
end

function _for_each_combination_continuum(values::Vector{Int}, k::Int, callback::Function)
    block = Vector{Int}(undef, k)
    last = length(values)
    function rec(pos::Int, start_idx::Int)
        if pos > k
            callback(block)
            return
        end
        @inbounds for idx in start_idx:last
            block[pos] = values[idx]
            rec(pos + 1, idx)
        end
    end
    k == 0 ? callback(Int[]) : rec(1, 1)
    return nothing
end

function _for_each_combination_continuum(callback::Function, values::Vector{Int}, k::Int)
    _for_each_combination_continuum(values, k, callback)
end

# enumerate an ensemble's blocks and dispatch to callback with a fresh copy
function _enumerate_ensemble(block_sizes::Vector{Int}, max_n::Int, is_continuum::Bool, callback::Function)
    if isempty(block_sizes)
        callback(Vector{Vector{Int}}())
        return
    end
    values = collect(1:max_n)
    blocks = Vector{Vector{Int}}()
    used = Int[]
    nblocks = length(block_sizes)

    function rec(block_idx::Int)
        if block_idx > nblocks
            callback(copy.(blocks))
            return
        end

        bsize = block_sizes[block_idx]
        if bsize == 0
            push!(blocks, Int[])
            rec(block_idx + 1)
            pop!(blocks)
            return
        end

        if is_continuum
            _for_each_combination_continuum(values, bsize) do block
                push!(blocks, copy(block))
                rec(block_idx + 1)
                pop!(blocks)
            end
        else
            available = setdiff(values, used)
            isempty(available) && return
            _for_each_combination_discrete(available, bsize) do block
                push!(blocks, copy(block))
                append!(used, block)
                rec(block_idx + 1)
                deleteat!(used, length(used)-bsize+1:length(used))
                pop!(blocks)
            end
        end
    end

    rec(1)
    return nothing
end

# allow do-block syntax that places callback first
function _enumerate_ensemble(callback::Function, block_sizes::Vector{Int}, max_n::Int, is_continuum::Bool)
    _enumerate_ensemble(block_sizes, max_n, is_continuum, callback)
end

# =========> Index Generation as Iterator <=================================================================================

function ensemble_iterator(block_sizes::Vector{Int}, max_n::Int)
    Channel() do ch
        _enumerate_ensemble(block_sizes, max_n, false) do blocks
            put!(ch, blocks)
        end
    end
end
function ensemble_iterator_continuum(block_sizes::Vector{Int}, max_n::Int)
    Channel() do ch
        _enumerate_ensemble(block_sizes, max_n, true) do blocks
            put!(ch, blocks)
        end
    end
end

"""
    multi_ensemble_iterator(block_sizes_list::Vector{Vector{Int}}, max_ns::Vector{Int}[, as_continuum::BitVector])

Lazy iterator over ensembles for each `(block_sizes, max_n)` pair. When
`as_continuum` is provided, it marks which ensembles allow repeated indices
(`i ≤ j` ordering). Use the `MultiEnsembleWorkspace` overload to pull both
`max_ns` and the continuum mask from an existing indexing workspace. The iterator
returns `Vector{Vector{Vector{Int}}}` objects that can be fed directly into the
ranking routines in `Indexing.jl`.
"""
function multi_ensemble_iterator(block_sizes_list::Vector{Vector{Int}}, max_ns::Vector{Int}, as_continuum::BitVector)
    n_ensembles = length(block_sizes_list)
    @assert length(max_ns) == n_ensembles
    @assert length(as_continuum) == n_ensembles

    Channel() do ch
        curr = Vector{Vector{Vector{Int}}}()
        function backtrack(idx::Int)
            if idx > n_ensembles
                put!(ch, deepcopy(curr))
                return
            end
            _enumerate_ensemble(block_sizes_list[idx], max_ns[idx], as_continuum[idx]) do blocks
                push!(curr, blocks)
                backtrack(idx + 1)
                pop!(curr)
            end
        end
        backtrack(1)
    end
end

function multi_ensemble_iterator(block_sizes_list::Vector{Vector{Int}}, max_ns::Vector{Int})
    return multi_ensemble_iterator(block_sizes_list, max_ns, falses(length(block_sizes_list)))
end

function multi_ensemble_iterator_continuum(block_sizes_list::Vector{Vector{Int}}, max_ns::Vector{Int})
    return multi_ensemble_iterator(block_sizes_list, max_ns, trues(length(block_sizes_list)))
end

function multi_ensemble_iterator(block_sizes_list::Vector{Vector{Int}}, workspace::MultiEnsembleWorkspace)
    n = length(block_sizes_list)
    max_ns = Vector{Int}(undef, n)
    @inbounds for i in 1:n
        max_ns[i] = workspace.rank_ws[i].bin_cache.max_n
    end
    return multi_ensemble_iterator(block_sizes_list, max_ns, workspace.as_continuum)
end
    
end
