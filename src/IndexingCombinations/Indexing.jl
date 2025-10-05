module Indexing 

export BinomialCache, _UsedBuf, EnsembleRankWorkspace, MultiEnsembleWorkspace

using Base.Threads
#### Define Structs >===================================================<
"""
    BinomialCache(max_n::Int, max_k::Int)

Cache of binomial coefficients C(n,k) for
    n ∈ [max_n - max_k, max_n]
    k ∈ [0, max_k], with k ≤ n

All values are stored as `Int`.  
⚠ Beware: `Int` will overflow if `C(n,k)` exceeds machine integer range.
"""
struct BinomialCache
    max_n::Int
    max_k::Int
    min_n::Int
    vals::Matrix{Int}  # (max_k+1, nrows)
    function BinomialCache(max_n::Int, max_k::Int)
        max_n < 0 && throw(ArgumentError("max_n must be ≥ 0"))
        max_k < 0 && throw(ArgumentError("max_k must be ≥ 0"))

        min_n = max_n - max_k
        nrows = max_n - min_n + 1
        vals  = zeros(Int, max_k+1, nrows)  # row = k+1, col = n-min_n+1

        # fill top row (n = max_n) using multiplicative formula
        col_top = nrows
        len_top = min(max_k, max_n) + 1
        vals[1, col_top] = 1
        @inbounds for k in 1:len_top-1
            prev = vals[k, col_top]
            vals[k+1, col_top] = (prev * (max_n - (k-1))) ÷ k
        end

        # fill downward
        @inbounds for n in (max_n-1):-1:min_n
            col     = n - min_n + 1
            col_abv = col + 1
            lenn    = min(max_k, n) + 1
            for k in 0:lenn-1
                vals[k+1, col] = (vals[k+1, col_abv] * (n + 1 - k)) ÷ (n + 1)
            end
        end

        return BinomialCache(max_n, max_k, min_n, vals)
    end
end
function BinomialCache(max_n::Int)::BinomialCache
    return BinomialCache(max_n, max_n)
end

@inline function _get(bin_cache::BinomialCache, n::Int, k::Int)
    bin_cache.vals[k+1, n - bin_cache.min_n + 1]
end

# accessors
(bin_cache::BinomialCache)(n::Int, k::Int) = _get(bin_cache, n, k)
Base.getindex(bin_cache::BinomialCache, n::Int, k::Int) = _get(bin_cache, n, k)


mutable struct _UsedBuf
    data::Vector{Int}  # capacity = bin_cache.max_k
    len::Int
end
@inline _usedbuf(cap::Int) = _UsedBuf(Vector{Int}(undef, cap), 0)
mutable struct EnsembleRankWorkspace
    bin_cache::BinomialCache
    used::_UsedBuf
    block_sizes::Vector{Int}
    ranks_by_blk::Vector{Int}
end
function EnsembleRankWorkspace(bin_cache::BinomialCache)
    # max #blocks ≤ max_k (each block must have ≥1 element)
    max_blocks = bin_cache.max_k
    EnsembleRankWorkspace(bin_cache, _usedbuf(bin_cache.max_k), Vector{Int}(undef, max_blocks), Vector{Int}(undef, max_blocks))
end

######## Multi-ensemble workspace ########

mutable struct MultiEnsembleWorkspace
    rank_ws::Vector{EnsembleRankWorkspace}
    sizes::Vector{Int}
    ranks::Vector{Int}
end

function MultiEnsembleWorkspace(bin_caches::Vector{BinomialCache})
    n = length(bin_caches)
    ensemble_work_space = Vector{EnsembleRankWorkspace}(undef, n)
    @inbounds for i in 1:n
        ensemble_work_space[i] = EnsembleRankWorkspace(bin_caches[i])
    end
    MultiEnsembleWorkspace(ensemble_work_space, Vector{Int}(undef, n), Vector{Int}(undef, n))
end
function MultiEnsembleWorkspace(max_ns::Vector{Int})::MultiEnsembleWorkspace
    # sort unique to not regenerate bin_caches so that we get the values and where they are 
    unique_max_ns = sort(unique(max_ns))
    bin_caches::Vector{BinomialCache} = BinomialCache[BinomialCache(max_n) for max_n in unique_max_ns]
    ordered_bin_caches::Vector{BinomialCache} = BinomialCache[]
    for max_n in max_ns 
        index = findfirst(x -> x == (max_n), unique_max_ns)
        push!(ordered_bin_caches, bin_caches[index]
    end
    return MultiEnsembleWorkspace(ordered_bin_caches)
end
function MultiEnsembleWorkspace(max_ns::Vector{Int}; max_threads::Int=1)::Vector{MultiEnsembleWorkspace}
    unique_max_ns = sort(unique(max_ns))
    bin_caches::Vector{BinomialCache} = BinomialCache[BinomialCache(max_n) for max_n in unique_max_ns]
    ordered_bin_caches::Vector{BinomialCache} = BinomialCache[]
    for max_n in max_ns 
        index = findfirst(x -> x == (max_n), unique_max_ns)
        push!(ordered_bin_caches, bin_caches[index]
    end
    return [MultiEnsembleWorkspace(ordered_bin_caches) for _ in 1:max_threads]
end


##### Finding indexes ################################################################################
function index_number_for_ensemble(block_sizes::Vector{Int}, bin_cache::BinomialCache)::Int
    # @assert sum(block_sizes) <= bin_cache.max_k "Cumulative block sizes exceed BinomialCache's max k value."
    if isempty(block_sizes)
        return 1 
    else
        curr_n = bin_cache.max_n
        total_indexes = bin_cache(curr_n, block_sizes[1])
        for i in 2:length(block_sizes)
            curr_n -= block_sizes[i-1]
            total_indexes *= bin_cache(curr_n, block_sizes[i])
        end
        return total_indexes 
    end
end

function bubble_insert_return_index!(used_indexes::Vector{Int}, val::Int, current_index::Int)::Int
    #pos = searchsortedfirst(used_indexes, val) # only makes sense for large vectors because of internal interval halving
    n = length(used_indexes)
    while current_index <= n && used_indexes[current_index] < val
        current_index += 1
    end
    insert!(used_indexes, current_index, val)
    return current_index+1
end


@inline function _insert_return_next!(ub::_UsedBuf, val::Int, current_index::Int)::Int
    @inbounds begin
        n = ub.len
        j = n
        # shift right until we find the insertion spot
        while j >= current_index && ub.data[j] > val
            ub.data[j+1] = ub.data[j]
            j -= 1
        end
        ub.data[j+1] = val
        ub.len = n + 1
        return j + 2 # new starting index
    end
end

@inline function index_rank_by_block!(blocked_indexes::Vector{Int}, ensemble_work_space::EnsembleRankWorkspace)::Int
    bc = ensemble_work_space.bin_cache
    ub = ensemble_work_space.used
    remaining_in_block   = length(blocked_indexes)
    total_index::Int     = 0
    prev_index::Int      = 0
    curr_index::Int      = 1

    # pointer to first used > prev_index  (so used_≤prev = up_prev-1)
    up_prev::Int = 1
    @inbounds while up_prev <= ub.len && ub.data[up_prev] <= prev_index
        up_prev += 1
    end

    @inbounds for idx in blocked_indexes  # sorted ascending
        # r = how many remain to pick after this position
        r = remaining_in_block - 1

        # counts before inserting idx
        # used_≤prev and used_≤idx via two advancing pointers
        up_b = up_prev
        @inbounds while up_b <= ub.len && ub.data[up_b] <= idx
            up_b += 1
        end
        used_le_prev = up_prev - 1
        used_le_b    = up_b   - 1
        U            = ub.len                  # total used so far (before inserting idx)
        A            = bc.max_n - U            # total free elements in compressed space

        # compressed coordinates: f(x) = x - used_≤x
        a = prev_index - used_le_prev
        b = idx        - used_le_b

        # boundary-term sum over unused candidates in (prev_index, idx):
        # S = Σ_{c unused} C(A - f(c), r) = C(A - a, r+1) - C(A - b + 1, r+1)
        total_index += bc(A - a, r + 1) - bc(A - b + 1, r + 1)

        # now insert the actual idx so later positions see it as used
        curr_index = _insert_return_next!(ub, idx, curr_index)

        # advance pointers for next iteration
        remaining_in_block -= 1
        prev_index = idx
        up_prev = up_b
        @inbounds while up_prev <= ub.len && ub.data[up_prev] <= prev_index
            up_prev += 1
        end
    end
    return total_index
end

@inline function index_rank_for_ensemble!(blocked_indexes::Vector{Vector{Int}}, ensemble_work_space::EnsembleRankWorkspace)::Int
    bc = ensemble_work_space.bin_cache
    if isempty(blocked_indexes); return 1; end

    nb = length(blocked_indexes)
    ensemble_work_space.used.len = 0

    @inbounds for i in 1:nb
        ensemble_work_space.block_sizes[i] = length(blocked_indexes[i])
    end

    @inbounds for i in 1:nb
        ensemble_work_space.ranks_by_blk[i] = index_rank_by_block!(blocked_indexes[i], ensemble_work_space)
    end

    curr_max::Int       = bc.max_n
    @inbounds for i in 1:nb
        curr_max -= ensemble_work_space.block_sizes[i]
    end

    combs_to_right::Int = 1
    total_rank::Int     = 0
    @inbounds for k in nb:-1:1
        bs = ensemble_work_space.block_sizes[k]
        r  = ensemble_work_space.ranks_by_blk[k]
        curr_max       += bs
        total_rank     += r * combs_to_right
        combs_to_right *= bc(curr_max, bs)
    end

    return total_rank + 1
end

"""
    index_number_for_ensemble(blocked_indexes::Vector{Vector{Int}}, bin_cache::BinomialCache) -> Int

Return the total number of possible configurations for an ensemble
described by `blocked_indexes`, using `bin_cache`.
"""
@inline function index_number_for_ensemble(blocked_indexes::Vector{Vector{Int}}, bin_cache::BinomialCache)::Int
    nb = length(blocked_indexes)
    if nb == 0
        return 1
    end

    curr_n        = bin_cache.max_n
    total_indexes = bin_cache(curr_n, length(blocked_indexes[1]))

    @inbounds for i in 2:nb
        curr_n       -= length(blocked_indexes[i-1])
        total_indexes *= bin_cache(curr_n, length(blocked_indexes[i]))
    end

    return total_indexes
end

"""
    combined_index_for_ensembles!(ensembles, mws) -> Int

- `ensembles :: Vector{Vector{Vector{Int}}}`  (each ensemble = blocks; blocks sorted)
- `mws :: MultiEnsembleWorkspace`             (owns bin_caches inside)

Returns 1-based mixed-radix index where the **last ensemble varies fastest**:
total = 1 + Σᵢ (rankᵢ - 1) * Π_{j>i} sizeⱼ
"""
@inline function combined_index_for_ensembles!(ensembles::Vector{Vector{Vector{Int}}},  multi_ensemble_workspace::MultiEnsembleWorkspace)::Int
    ne = length(ensembles)

    @inbounds for i in 1:ne
        ensemble_work_space = multi_ensemble_workspace.rank_ws[i]
        multi_ensemble_workspace.ranks[i] = index_rank_for_ensemble!(ensembles[i], ensemble_work_space)
        multi_ensemble_workspace.sizes[i] = index_number_for_ensemble(ensembles[i], ensemble_work_space.bin_cache)
    end

    total::Int  = 1
    factor::Int = 1
    @inbounds for i in ne:-1:1
        total  += (multi_ensemble_workspace.ranks[i] - 1) * factor
        factor *= multi_ensemble_workspace.sizes[i]
    end
    return total
end

end