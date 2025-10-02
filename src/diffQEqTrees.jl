module EqTrees

using Base: OneTo, LinearIndices, CartesianIndices
using ..QExpressions: diffQEq, QAbstract, diffQEqOrdered, QCumulantOrdered, OrderbyOperator, QAtomProduct, QExpr, QTerm, contains_which_t_indexes, max_order_of_terms, which_abstracts
using ..QSpaces: QSpace, SubSpace
using ..OffsetArrays: OffsetArray

const OpIndex = Vector{Int}
const OpBlock = Vector{OpIndex}
const OpPath = Vector{OpBlock}

const StoredTypes = Union{diffQEqOrdered, QCumulantOrdered}
# =====================> Tree Start <===========================================================
abstract type AbstractOpIndexNode{P} end  # P = payload type in leaves
struct OpIndexNode{P} <: AbstractOpIndexNode{P}
    children::OffsetArray{AbstractOpIndexNode{P}}
end
mutable struct OpIndexLeaf{P} <: AbstractOpIndexNode{P}
    payload::P
    how_many::Int
end
struct OpIndexRoot{P} <: AbstractOpIndexNode{P}
    where_ensembles::Vector{Int}
    subspaces::Vector{SubSpace}
    children::OffsetArray{AbstractOpIndexNode{P}}
end

struct diffQEqTree{N}
    qspace::QSpace
    diff_qeq::diffQEq
    abstract::QAbstract
    max_order::Int
    where_ensembles::Vector{Int}
    diff_trees::OpIndexRoot{diffQEqOrdered}
    cumulant_trees::OpIndexRoot{QCumulantOrdered}
    function diffQEqTree(qspace::QSpace, diff_qeq::diffQEq, abstract::QAbstract, max_order::Int)
        where_ensembles::Vector{Int} = qspace.subspace_info.where_ensembles
        how_many_ensembles::Int = length(where_ensembles)
        z = zeros(Int, how_many_ensembles)
        diff_trees = OpIndexRoot{diffQEqOrdered}(where_ensembles, qspace.subspaces, OffsetArray{OpIndexNode{diffQEqOrdered}}(z, z))
        cumulant_trees = OpIndexRoot{QCumulantOrdered}(where_ensembles, qspace.subspaces, OffsetArray{OpIndexNode{QCumulantOrdered}}(z, z))
        return diffQEqTree{N}(qspace, diff_qeq, abstract, max_order, where_ensembles, diff_trees, cumulant_trees)
    end
end

@inline function Base.getindex(root::OpIndexRoot{P}, pos::OpPath) where {P <: StoredTypes}
    # first step: pick the correct entry from the root by ensemble lengths
    start_idxs = Tuple(length(pos[i]) for i in root.where_ensembles)
    curr = root.children[start_idxs...]
    curr === nothing && return nothing

    # then descend through blocks
    @inbounds for block in pos
        for op_ind in block
            curr isa OpIndexNode{P} || return nothing
            curr = curr.children[op_ind...]
            curr === nothing && return nothing
        end
    end
    # final node must be a leaf
    @assert isa(curr, OpIndexLeaf{P}) "The final element should have been a OpIndexLeaf, but got a $(typeof(curr)) instead. "
    return curr 
end

# Child grid sized by the *current* subspace's bounds.
@inline function _grid_for(::Type{P}, ss::SubSpace) where {P <: StoredTypes}
    m, M = ss.min_ints, ss.max_ints
    #length(m) == length(M) || error("SubSpace.min_ints / max_ints length mismatch")
    return OffsetArray{AbstractOpIndexNode{P}}(m, M)  # -1 in max_ints -> 0 length init (your rule)
end
# Promote/ensure a node at a position; its children are sized by `next_ss`.
@inline function _ensure_node!(kids::OffsetArray{AbstractOpIndexNode{P}},idxs::AbstractVector{Int},ss::SubSpace,) where {P <: StoredTypes}
    ex = kids[idxs...]  # may be nothing / leaf / node
    if ex isa OpIndexNode{P}
        return ex
    else
        nd = OpIndexNode{P}(_grid_for(P, ss))
        kids[idxs...] = nd
        return nd
    end
end
# Place/update a leaf at a position.
@inline function _place_leaf!(kids::OffsetArray{AbstractOpIndexNode{P}}, idxs::AbstractVector{Int}, payload::P, how_many::Int) where {P <: StoredTypes}
    ex = kids[idxs...]
    if isnothing(ex) || ex isa OpIndexLeaf{P}
        kids[idxs...] = OpIndexLeaf{P}(payload, how_many)
    else#if ex isa OpIndexNode{P}
        error("Expected OpIndexLeaf or Nothing at end of path, found $(typeof(ex)).")
    end
    return payload
end

function Base.setindex!(root::OpIndexRoot{P}, payload::P, how_many::Int,::P, pos::OpPath) where {P <: StoredTypes}
    start_idxs = Tuple(length(pos[i]) for i in root.where_ensembles)
    curr = root.children[start_idxs...]

    @inbounds for (block, ss) in zip(pos[1:end-1], root.subspaces[1:end-1])
        for op_ind in block 
            kids = (curr::OpIndexNode{P}).children
            curr = _ensure_node!(kids, op_ind, ss)
        end
    end
    block, ss = pos[end], root.subspaces[end]
    for op_ind in block[1:end-1]
        kids = (curr::OpIndexNode{P}).children
        curr = _ensure_node!(kids, op_ind, ss)
    end            
    # final_step
    _place_leaf!(kids, op_ind, payload, how_many)
    return payload
end

# =======================================> Generate EQ Sets <==================================================

function gen_diffQEqSet(diff_qeq::diffQEq, initial_operators::Vector{QExpr}, max_order::Int=2)::diffQEqTree
    extracted_initial_operators = QAtomProduct[]
    for init in initial_operators
        @assert length(init) == 1 "initial_operators need a single term, got $init instead."
        term = init.terms[1]
        @assert term isa QAtomProduct "initial_operators must contain a single QAtomProduct, got $(typeof(term)) instead."
        push!(extracted_initial_operators, term)
    end
    return gen_diffQEqSet(diff_qeq, extracted_initial_operators, max_order)
end

function gen_diffQEqSet(diff_qeq::diffQEq, initial_operators::Vector{QAtomProduct}, max_order::Int=2)::diffQEqTree
    for expr in initial_operators
        @assert length(expr.expr) == 1 && expr.expr[1] isa QTerm "initial_operators must contain a single QTerm inside each QAtomProduct; got $(expr.expr) instead."
    end

    qspace = diff_qeq.qspace
    time_indexes = contains_which_t_indexes(diff_qeq)
    @assert sum(time_indexes) <= 1 "Differential equation contains too many time indexes, found time indexes: $(findall(time_indexes) .- 1)."
    first_time = findfirst(time_indexes)
    curr_time_index = isnothing(first_time) ? -1 : first_time - 1

    isempty(initial_operators) || @assert maximum(max_order_of_terms.(initial_operators)) <= max_order "initial_operators contain QTerms that exceed max_order"

    defined_abstracts = which_abstracts(diff_qeq.left_hand_side)
    lengths_abstract = length.(defined_abstracts)
    @assert sum(lengths_abstract) == 1 "diff_qeq requires exactly one QAbstract."
    key_index = findfirst(>(0), lengths_abstract)
    sub_index = defined_abstracts[key_index][1]
    operator_type = qspace.operatortypes[key_index]
    abstract = QAbstract(operator_type, key_index, sub_index, 1, false, curr_time_index)

    ordered_diff = OrderbyOperator(diff_qeq)
    tree = diffQEqTree(qspace; source_equation=ordered_diff, root_abstract=abstract, max_order=max_order)

    # now for each initial_operator, substitute abstract in diff_qeq with it, then add it to the tree 
    # (but store it as the transformed diffQEqOrdered as well as the number of terms that we have for it), 
    # then for each qterm in it also use it to substitute into diff_qeq and then add the result to the tree until we reach max_order, 
    # terms that surpass max_order we add to a second tree in which we store cumulants
end
