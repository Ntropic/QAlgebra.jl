module EqTrees

struct OpIndexNode
    op_index::Int
    children::Vector{OpIndexNode}
    diff_QEq::Union{Nothing, diff_QEq}  # only valid at leaves
    how_many::Int
end

OpIndexNode(op_index::Int) = OpIndexNode(op_index, OpIndexNode[], nothing)

mutable struct OpIndexTree
    root::OpIndexNode
end

OpIndexTree() = OpIndexTree(OpIndexNode(-1))  # sentinel root

mutable struct diff_QEqTree
    op_order::Vector{Int}
    trees::Vector{OpIndexTree}
end

diff_QEqTree(op_order::Vector{Int}) = diff_QEqTree(op_order, [OpIndexTree() for _ in 1:length(op_order)])

"Traverse the tree with a Vector{Vector{Int}} path"
function navigate(tree::OpIndexTree, path::Vector{Vector{Int}}; create=false)
    node = tree.root
    for segment in path
        isempty(segment) && continue   # skip empty vectors

        # Step 1: go to child by segment length
        len = length(segment)
        child = findfirst(c -> c.op_index == len, node.children)
        if child === nothing
            if create
                newchild = OpIndexNode(len)
                push!(node.children, newchild)
                node = newchild
            else
                return nothing
            end
        else
            node = node.children[child]
        end

        # Step 2: traverse inner vector contents
        for idx in segment
            child = findfirst(c -> c.op_index == idx, node.children)
            if child === nothing
                if create
                    newchild = OpIndexNode(idx)
                    push!(node.children, newchild)
                    node = newchild
                else
                    return nothing
                end
            else
                node = node.children[child]
            end
        end
    end
    return node
end

function Base.setindex!(dqt::diff_QEqTree, eq::String, order::Int, path::Vector{Vector{Int}})
    idx = findfirst(==(order), dqt.op_order)
    idx === nothing && error("Order $order not found in op_order")

    node = navigate(dqt.trees[idx], path; create=true)
    if !isempty(node.children)
        error("Cannot assign diff_QEq to non-leaf node")
    end
    node.diff_QEq = eq
end

function Base.getindex(dqt::diff_QEqTree, order::Int, path::Vector{Vector{Int}})
    idx = findfirst(==(order), dqt.op_order)
    idx === nothing && return nothing

    node = navigate(dqt.trees[idx], path; create=false)
    return node === nothing ? nothing : node.diff_QEq
end

function Base.haskey(dqt::diff_QEqTree, key::Tuple{Int,Vector{Vector{Int}}})
    order, path = key
    idx = findfirst(==(order), dqt.op_order)
    idx === nothing && return false
    node = navigate(dqt.trees[idx], path; create=false)
    return node !== nothing && node.diff_QEq !== nothing
end
