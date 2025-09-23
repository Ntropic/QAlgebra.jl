module Cumulants

import Base: show, string

function signed_int_string(i::Int)::String
    if i < -1
        return " -" * string(-i) * " "
    elseif i == -1
        return " - "
    elseif i == 0
        return " +0"
    elseif i == 1
        return " + "
    else
        return " +" * string(i) * " "
    end
end

intvec_to_do_braket(is::Vector{Int})::String = "<" * join(is, ",") * ">"

struct IndexedProduct
    coeff::Int
    indices::Vector{Vector{Int}}
end

function string(ind_prod::IndexedProduct)::String
    return signed_int_string(ind_prod.coeff) * join(intvec_to_do_braket.(ind_prod.indices), "")
end
function show(io::IO, ind_prod::IndexedProduct)
    print(io, "IndexedProduct: ", string(ind_prod))
end

struct IndexedCumulant
    n::Int
    cumulant::Vector{IndexedProduct}
end

function string(ind_cum::IndexedCumulant)::String
    return join(string.(ind_cum.cumulant), "")
end
function show(io::IO, ind_cum::IndexedCumulant)
    print(io, "IndexedCumulant: ", string(ind_cum))
end

import Base: -
(-(a::IndexedProduct)) = IndexedProduct(-a.coeff, a.indices)

function FirstIndexedCumulant(i::Int=1)
    return IndexedCumulant(1, [IndexedProduct(1, [[i]])])
end

function raising_operator(new_op::Int, prev_product::IndexedProduct)::Vector{IndexedProduct}
    n = length(prev_product.indices)
    prev_coeff = prev_product.coeff

    new_product = IndexedProduct[]
    for i in 1:length(prev_product.indices)
        new_c = [copy(v) for v in prev_product.indices]
        push!(new_c[i], new_op)
        push!(new_product, IndexedProduct(prev_coeff, new_c))
    end

    c = [copy(v) for v in prev_product.indices]
    push!(c, [new_op])
    push!(new_product, IndexedProduct(-n * prev_coeff, c))
    return new_product
end

function raising_operator(new_op::Int, prev_cumulant::IndexedCumulant)
    cumulant = prev_cumulant.cumulant
    new_cumulant = raising_operator(new_op, cumulant[1])
    for i in 2:length(cumulant)
        append!(new_cumulant, raising_operator(new_op, cumulant[i]))
    end
    return IndexedCumulant(prev_cumulant.n + 1, new_cumulant)
end

function IndexedCumulant(order::Int)
    order <= 0 && error("Order must be greater than 0")
    cum = FirstIndexedCumulant()
    for i in 2:order
        cum = raising_operator(i, cum)
    end
    return cum
end

struct ReducedIndexedCumulant
    operator::IndexedProduct
    approximation::Vector{IndexedProduct}
end

function ReducedIndexedCumulant(order::Int)
    full_cumulant = IndexedCumulant(order)
    operator = full_cumulant.cumulant[1]
    approximation = .-full_cumulant.cumulant[2:end]
    return ReducedIndexedCumulant(operator, approximation)
end

function string(cumulant_list::ReducedIndexedCumulant)
    return string(cumulant_list.operator) * " ≈ " * join(string.(cumulant_list.approximation), "")
end
function show(io::IO, cumulant_list::ReducedIndexedCumulant)
    print(io, "ReducedIndexedCumulant: ", string(cumulant_list))
end

mutable struct ReducedCumulantList
    curr_order::Int
    max_order::Int
    last_cumulant::IndexedCumulant
    reduced_cumulants::Vector{ReducedIndexedCumulant}
end

function string(cumulant_list::ReducedCumulantList)
    return join(string.(cumulant_list.reduced_cumulants), "\n")
end
function show(io::IO, cumulant_list::ReducedCumulantList)
    print(io, "ReducedCumulantList with curr_order = ", cumulant_list.curr_order, "\n", string(cumulant_list))
end

function ReducedCumulantList(curr_order::Int=1; max_order::Int=10^12)
    curr_order <= 0 && error("Current order must be greater than 0")
    curr_order > max_order && error("Current order must be less than or equal to max order")

    reduced_cumulants = ReducedIndexedCumulant[]
    full_cumulant = FirstIndexedCumulant()
    operator = full_cumulant.cumulant[1]
    approximation = .-full_cumulant.cumulant[2:end]
    push!(reduced_cumulants, ReducedIndexedCumulant(operator, approximation))

    for i in 2:curr_order
        full_cumulant = raising_operator(i, full_cumulant)
        operator = full_cumulant.cumulant[1]
        approximation = .-full_cumulant.cumulant[2:end]
        push!(reduced_cumulants, ReducedIndexedCumulant(operator, approximation))
    end

    return ReducedCumulantList(curr_order, max_order, full_cumulant, reduced_cumulants)
end

function (cumulant_list::ReducedCumulantList)(order::Int)
    if order > cumulant_list.curr_order
        expand_cumulant_list!(cumulant_list, order)
    end
    return cumulant_list.reduced_cumulants[order]
end

function expand_cumulant_list!(cumulant_list::ReducedCumulantList, order::Int)::Nothing
    order > cumulant_list.max_order && error("Order must be less than or equal to max_order")
    order <= cumulant_list.curr_order && (@warn "Order can only be increased if target order is greater than current order"; return)

    full_cumulant = cumulant_list.last_cumulant
    for i in cumulant_list.curr_order + 1:order
        full_cumulant = raising_operator(i, full_cumulant)
        operator = full_cumulant.cumulant[1]
        approximation = .-full_cumulant.cumulant[2:end]
        push!(cumulant_list.reduced_cumulants, ReducedIndexedCumulant(operator, approximation))
    end
    cumulant_list.curr_order = order
    cumulant_list.last_cumulant = full_cumulant
    return
end

export IndexedProduct, IndexedCumulant, ReducedIndexedCumulant, ReducedCumulantList,
       FirstIndexedCumulant, raising_operator, expand_cumulant_list!

end # module Cumulants
