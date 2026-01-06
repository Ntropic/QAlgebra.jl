export sorted_push!, sorted_push,
       sorted_push_unique!, sorted_push_unique,
       sorted_append!, sorted_append,
       sorted_append_unique!, sorted_append_unique,
       sorted_append_unify!, sorted_append_unify,
       sorted_push_unify_branches, sorted_append_unify_branches

using ComplexRationals

#####################
# Internal utilities #
#####################

# Clamp first_n into [1, n] (for n > 0).
@inline function _clamp_first_n(first_n::Int, n::Int)
    first_n <= 1 ? 1 : (first_n > n ? n : first_n)
end

# Equality induced by lt: a == b  ⇔  !lt(a, b) && !lt(b, a)
@inline _equal_via_lt(a, b, lt) = !lt(a, b) && !lt(b, a)

# Find insertion position for x in dest, starting search at first_n.
# Returns the index at which x should be inserted to preserve sort order.
@inline function _sorted_search_insert_position(dest::Vector{T}, x::T,
                                                lt, first_n::Int) where {T}
    n = length(dest)
    n == 0 && return 1

    i0 = _clamp_first_n(first_n, n)

    @inbounds for i in i0:n
        if lt(x, dest[i])
            return i
        end
    end

    return n + 1
end

# Find either an equal element (via lt) or an insertion position.
# Returns (found_equal::Bool, idx::Int).
@inline function _sorted_search_equal_or_insert(dest::Vector{T}, x::T,
                                                lt, first_n::Int) where {T}
    n = length(dest)
    n == 0 && return false, 1

    i0 = _clamp_first_n(first_n, n)
    pos = n + 1

    @inbounds for i in i0:n
        di = dest[i]
        if _equal_via_lt(di, x, lt)
            return true, i
        elseif lt(x, di)
            pos = i
            break
        end
    end

    return false, pos
end

# Internal: push with unification.
# If an equal element exists, unify(existing, x) is stored in its place,
# otherwise x is inserted. Returns a tuple `(idx, inserted)` where `idx`
# is the position of the resulting element and `inserted` indicates
# whether a new slot was created.
@inline function _sorted_push_unify!(dest::Vector{T}, x::T;
                                     lt=isless, unify, first_n::Int=1) where {T}
    n = length(dest)
    found, pos = _sorted_search_equal_or_insert(dest, x, lt, first_n)

    if found
        dest[pos] = unify(dest[pos], x)
        return pos, false
    end

    push!(dest, x)
    @inbounds for i in n:-1:pos
        dest[i + 1] = dest[i]
    end
    dest[pos] = x

    return pos, true
end

struct _BranchState{T}
    weight::ComplexRational
    values::Vector{T}
    first_n::Int
end

function _branch_insert_one(state::_BranchState{T}, x::T, lt, unify) where {T}
    vec = state.values
    found, pos = _sorted_search_equal_or_insert(vec, x, lt, state.first_n)

    if found
        updates = unify(vec[pos], x)
        new_states = Vector{_BranchState{T}}()
        for (factor, new_val) in updates
            updated = copy(vec)
            updated[pos] = new_val
            push!(new_states, _BranchState(state.weight * factor, updated, max(1, pos - 1)))
        end
        return new_states
    end

    n = length(vec)
    updated = Vector{T}(undef, n + 1)
    @inbounds for i in 1:pos-1
        updated[i] = vec[i]
    end
    @inbounds for i in pos:n
        updated[i + 1] = vec[i]
    end
    updated[pos] = x
    return [_BranchState(state.weight, updated, pos)]
end

function _branch_insert(states::Vector{_BranchState{T}}, x::T, lt, unify) where {T}
    out = Vector{_BranchState{T}}()
    for state in states
        append!(out, _branch_insert_one(state, x, lt, unify))
    end
    return out
end

_branch_results(states::Vector{_BranchState{T}}) where {T} =
    [(state.weight, state.values) for state in states]

"""
    sorted_push_unify_branches(dest::Vector{T}, x::T;
                               lt=isless, first_n::Int=1, unify) -> Vector{Tuple{ComplexRational, Vector{T}}}

Insert `x` into the sorted vector `dest`, branching whenever `x` unifies with an
existing element. The `unify(existing, x)` function must return a collection of
`(weight::ComplexRational, value::T)` tuples; each tuple represents one outcome
that replaces the existing element. The returned value is a vector of
`(weight, vector)` tuples, one per branch, where `vector` is the resulting data
and `weight` is the accumulated product of branch weights along that path.
"""
function sorted_push_unify_branches(dest::Vector{T}, x::T;
                                    lt=isless, first_n::Int=1, unify) where {T}
    initial = _BranchState(ComplexRational(1, 0, 1), copy(dest), first_n)
    states = _branch_insert([initial], x, lt, unify)
    return _branch_results(states)
end

"""
    sorted_append_unify_branches(dest::Vector{T}, src::AbstractVector{<:T};
                                 lt=isless, first_n::Int=1, unify)

Append the sorted sequence `src` into `dest`, branching whenever a unification
produces multiple outcomes. Each branch tracks its own insertion index to seed
the next search, mirroring the behaviour of [`sorted_append!`](@ref). Returns
the same `(weight, vector)` collection as [`sorted_push_unify_branches`](@ref).
"""
function sorted_append_unify_branches(dest::Vector{T}, src::AbstractVector{<:T};
                                      lt=isless, first_n::Int=1, unify) where {T}
    states = [_BranchState(ComplexRational(1, 0, 1), copy(dest), first_n)]
    for x in src
        states = _branch_insert(states, x, lt, unify)
        isempty(states) && break
    end
    return _branch_results(states)
end


############################
# sorted_push / sorted_push!
############################

"""
    sorted_push!(dest::Vector{T}, x::T;
                 lt=isless, first_n::Int=1) -> Int

Insert `x` into sorted vector `dest` in place, preserving sort order
(allowing duplicates).

The search for the insertion position starts at index `first_n`. This is
useful when performing multiple pushes where insertion positions are known
to be non-decreasing (e.g. when inserting elements from another sorted
vector in ascending order).

Returns the insertion index of `x` in `dest` after insertion.

**Precondition:** if you pass `first_n > 1`, you must ensure that `x` is
greater or equal (with respect to `lt`) than all elements in
`dest[1:first_n-1]`, otherwise the result may not be sorted.
"""
function sorted_push!(dest::Vector{T}, x::T;
                      lt=isless, first_n::Int=1) where {T}
    n = length(dest)
    pos = _sorted_search_insert_position(dest, x, lt, first_n)

    push!(dest, x)
    @inbounds for i in n:-1:pos
        dest[i + 1] = dest[i]
    end
    dest[pos] = x

    return pos
end

"""
    sorted_push(dest::Vector{T}, x::T;
                lt=isless, first_n::Int=1) -> Vector{T}

Pure version of [`sorted_push!`](@ref). Returns a new sorted vector
containing `dest` with `x` inserted, allowing duplicates.
"""
function sorted_push(dest::Vector{T}, x::T;
                     lt=isless, first_n::Int=1) where {T}
    out = copy(dest)
    sorted_push!(out, x; lt=lt, first_n=first_n)
    return out
end


##################################
# sorted_push_unique / sorted_push_unique!
##################################

"""
    sorted_push_unique!(dest::Vector{T}, x::T;
                        lt=isless, first_n::Int=1) -> Int

Insert `x` into sorted vector `dest` in place if it is not already present.

Equality is defined by `lt` as

    a and b are equal  ⇔  !lt(a, b) && !lt(b, a).

The search starts at `first_n`. If `first_n > 1`, you must ensure there is
no element equal to `x` (under the `lt`-based equality) in `dest[1:first_n-1]`
if you rely on the uniqueness guarantee.

Returns the index of the element equal to `x` in `dest` after the call:

  - If `x` was already present, this is the index of the existing element.
  - If `x` was inserted, this is its insertion index.
"""
function sorted_push_unique!(dest::Vector{T}, x::T;
                             lt=isless, first_n::Int=1) where {T}
    n = length(dest)
    found, pos = _sorted_search_equal_or_insert(dest, x, lt, first_n)

    if found
        return pos
    end

    push!(dest, x)
    @inbounds for i in n:-1:pos
        dest[i + 1] = dest[i]
    end
    dest[pos] = x

    return pos
end

"""
    sorted_push_unique(dest::Vector{T}, x::T;
                       lt=isless, first_n::Int=1) -> Vector{T}

Pure unique insert: return a new sorted vector containing `dest` with `x`
inserted if it is not already present (according to the `lt`-based equality).
"""
function sorted_push_unique(dest::Vector{T}, x::T;
                            lt=isless, first_n::Int=1) where {T}
    out = copy(dest)
    sorted_push_unique!(out, x; lt=lt, first_n=first_n)
    return out
end


##############################
# sorted_append / sorted_append!
##############################

"""
    sorted_append!(dest::Vector{T}, src::AbstractVector{<:T};
                   lt=isless) -> Vector{T}

In-place merge of two sorted sequences `dest` and `src` into `dest`,
preserving all duplicates.

This is implemented by repeatedly calling [`sorted_push!`](@ref) for each
element of `src`, using the returned insertion index to seed `first_n` for
the next insertion.

- Both `dest` and `src` must already be sorted according to `lt`.
- `src` must not alias `dest`.

Returns `dest`.
"""
function sorted_append!(dest::Vector{T}, src::AbstractVector{<:T};
                        lt=isless) where {T}
    isempty(src) && return dest

    first_n = 1
    @inbounds for x in src
        first_n = sorted_push!(dest, x; lt=lt, first_n=first_n)
    end

    return dest
end

"""
    sorted_append(dest::Vector{T}, src::AbstractVector{<:T};
                  lt=isless) -> Vector{T}

Pure merge of two sorted vectors, allowing duplicates.
Returns a new vector containing all elements from `dest` and `src`
in sorted order.
"""
function sorted_append(dest::Vector{T}, src::AbstractVector{<:T};
                       lt=isless) where {T}
    out = copy(dest)
    sorted_append!(out, src; lt=lt)
    return out
end


######################################
# sorted_append_unique / sorted_append_unique!
######################################

"""
    sorted_append_unique!(dest::Vector{T}, src::AbstractVector{<:T};
                          lt=isless) -> Vector{T}

In-place merge of two sorted sequences `dest` and `src` into `dest`, while
removing duplicates according to the equality implied by `lt`:

    a and b are equal  ⇔  !lt(a, b) && !lt(b, a).

This is implemented via [`sorted_append_unify!`](@ref) with a unifier
that keeps the existing representative.

- Both `dest` and `src` must be sorted according to `lt`.
- `src` must not alias `dest`.
- If `src` is empty, `dest` is returned unchanged (even if it contains
  duplicates).
- If `src` is non-empty, duplicates across `dest ∪ src` are removed.

Returns `dest`.
"""
function sorted_append_unique!(dest::Vector{T}, src::AbstractVector{<:T};
                               lt=isless) where {T}
    isempty(src) && return dest

    # Unifier that keeps the existing representative.
    sorted_append_unify!(dest, src, (existing, _new) -> existing; lt=lt)
    return dest
end

"""
    sorted_append_unique(dest::Vector{T}, src::AbstractVector{<:T};
                         lt=isless) -> Vector{T}

Pure merge of two sorted vectors while removing duplicates according to
the `lt`-based equality. Returns a new sorted vector with all unique
elements from `dest` and `src`.
"""
function sorted_append_unique(dest::Vector{T}, src::AbstractVector{<:T};
                              lt=isless) where {T}
    out = copy(dest)
    sorted_append_unique!(out, src; lt=lt)
    return out
end


########################################
# sorted_append_unify / sorted_append_unify!
########################################

"""
    sorted_append_unify!(dest::Vector{T}, src::AbstractVector{<:T}, unify::Function;
                         lt=isless) -> Vector{T}

In-place merge of two sorted sequences `dest` and `src` into `dest`, while
*unifying* equal elements (according to `lt`) using the function `unify`.

**Requirements on `unify`:**

- It must have the signature `unify(existing, new)`.
- For all elements `a, b` in the same equivalence class, the value returned
  by repeatedly applying `unify` must remain equal (under `lt`) to all
  inputs; i.e. it must not break the sorting order.

Implementation details:

- `src` must not alias `dest`.
- If `dest` is non-empty, it is copied, cleared, and rebuilt by repeatedly
  calling an internal `_sorted_push_unify!` for all elements from the old
  `dest` and then from `src`.
- If `dest` is empty, only elements from `src` are unified and inserted.

Returns `dest`.
"""
function sorted_append_unify!(dest::Vector{T}, src::AbstractVector{<:T}, unify::Function;
                              lt=isless) where {T}
    src === dest && throw(ArgumentError("src must not alias dest"))

    first_n = 1

    if !isempty(dest)
        old = copy(dest)
        resize!(dest, 0)

        @inbounds for x in old
            pos, inserted = _sorted_push_unify!(dest, x; lt=lt, unify=unify, first_n=first_n)
            first_n = inserted ? pos : max(1, pos - 1)
        end

        # Reset search start before processing `src` since they may contain elements
        # smaller than the last reinserted element from `dest`.
        first_n = 1
    end

    @inbounds for y in src
        pos, inserted = _sorted_push_unify!(dest, y; lt=lt, unify=unify, first_n=first_n)
        first_n = inserted ? pos : max(1, pos - 1)
    end

    return dest
end

"""
    sorted_append_unify(dest::Vector{T}, src::AbstractVector{<:T}, unify::Function;
                        lt=isless) -> Vector{T}

Pure version of [`sorted_append_unify!`](@ref).

Returns a new sorted vector that is the unified merge of `dest` and `src`.
"""
function sorted_append_unify(dest::Vector{T}, src::AbstractVector{<:T}, unify::Function;
                             lt=isless) where {T}
    out = copy(dest)
    sorted_append_unify!(out, src, unify; lt=lt)
    return out
end
