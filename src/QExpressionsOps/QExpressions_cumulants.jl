#### Preparations ---> Index Cumulants first, then QExpressions further down 
import Base: show, string 
function signed_int_string(i::Int)::String
    if i < -1
        return " -" * string(-i)* " "
    elseif i == -1
        return " - "  
    elseif i == 0
        return " +0" 
    elseif i == 1
        return " + " 
    else
        return " +" * string(i)* " "
    end
end
function intvec_to_braket(is::Vector{Int})::String
    return string("<", join(is, ","), ">")
end

# Define necessary Structs
struct IndexedProduct
    coeff::Int 
    indices::Vector{Vector{Int}}
end

function string(ind_prod::IndexedProduct)::String
    return signed_int_string(ind_prod.coeff)*join(intvec_to_braket.(ind_prod.indices),"")
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

import Base:-
function -(a::IndexedProduct)::IndexedProduct 
    return IndexedProduct(-a.coeff, a.indices)
end

# Define Operations to Construct and Raise the Level of the Structs!
function FirstIndexedCumulant(i::Int=1)
    return IndexedCumulant(1, [IndexedProduct(1, [[i]])])
end
function raising_operator(new_op::Int, prev_product::IndexedProduct)::Vector{IndexedProduct}
    # Raising operator is distributve with respect to products
    # new_op | [i1, i2, i3] = [new_op, i1, i2, i3] - [i1, i2, i3] [new_op]
    # new_op | [[A],[B],...] = [new_op | [A], new_op | [B], ...]
    # The sum of applying itto each Vector{Int} in the product
    n = length(prev_product.indices)
    prev_coeff = prev_product.coeff

    new_product::Vector{IndexedProduct} = []
    for i in 1:length(prev_product.indices)
        new_c = [ copy(v) for v in prev_product.indices ]
        push!(new_c[i], new_op)
        push!(new_product, IndexedProduct(prev_coeff, new_c))
    end
    c = [ copy(v) for v in prev_product.indices ]
    push!(c, [new_op])
    push!(new_product, IndexedProduct(-n*prev_coeff, c))
    return new_product 
end

function raising_operator(new_op::Int, prev_cumulant::IndexedCumulant)
    # Raising operator is distributve with respect to products  
    # and linear with respect to sum   
    cumulant = prev_cumulant.cumulant
    new_cumulant = raising_operator(new_op, cumulant[1]) 
    for i in 2:length(cumulant)
        append!(new_cumulant, raising_operator(new_op, cumulant[i]))
    end
    return IndexedCumulant(prev_cumulant.n+1, new_cumulant)
end
function IndexedCumulant(order::Int)
    if order == 0 
        error("Order must be greater than 0")
    end
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

# Create Cumulant list 
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

function ReducedCumulantList(curr_order::Int=1;max_order::Int=10^12)
    if curr_order > max_order
        error("Current order must be less than or equal to max order")
    end
    if curr_order > 0 
        reduced_cumulants::Vector{ReducedIndexedCumulant} = []
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
    else
        error("Current order must be greater than 0")
    end
end
# Access as a function i.e. reduced_cumulant_list(3) gives the 3rd order cumulant and if necessary constructs it up to that order
function (cumulant_list::ReducedCumulantList)(order::Int)
    if order > cumulant_list.curr_order
        expand_cumulant_list!(cumulant_list, order)
    end
    return cumulant_list.reduced_cumulants[order]
end

function expand_cumulant_list!(cumulant_list::ReducedCumulantList, order::Int)::Nothing
    if order > cumulant_list.max_order
        error("Order must be less than or equal to max_order")
    end
    if order <= cumulant_list.curr_order
        warning("Order can only be increased if target order is greater than current order")
        return
    end
    full_cumulant = cumulant_list.last_cumulant
    for i in cumulant_list.curr_order+1:order
        full_cumulant = raising_operator(i, full_cumulant)
        operator = full_cumulant.cumulant[1]
        approximation = .-full_cumulant.cumulant[2:end]
        push!(cumulant_list.reduced_cumulants, ReducedIndexedCumulant(operator, approximation))
    end
    cumulant_list.curr_order = order
    cumulant_list.last_cumulant = full_cumulant
    return
end

function where_acting_index(q::QTerm, statespace::StateSpace)::Vector{Int}
    return [i for (i, op) in enumerate(q.op_indices) if op!=statespace.neutral_op[i]]
end
# Order of QTerm 
function order(q::QTerm, statespace::StateSpace)::Int
    # Determines the order of a QTerm operator 
    return length(where_acting_index(q, statespace))
end

struct QCumulant <:QComposite  # Must be QComposite to be in qExpr's
    statespace::StateSpace
    atom::QAtom
    expr_::qExpr
    order::Int
    where_acting::Vector{Int}
end

#Is = Union{Int,Vector{Int}}
function replace_indexes(neutral_op::Vector{Is}, curr_op_indexes::Vector{Is}, indexes::Vector{Int}) 
    new_op = copy(neutral_op)
    for ind in indexes
        new_op[ind] = copy(curr_op_indexes[ind])
    end
    return new_op
end

# Add conspiracies later!
function Cumulant(qprod::QAtomProduct)
    if length(qprod.terms) != 1
        error("Cumulant only defined for single QAtom in QAtomProduct")
    end
    return Cumulant(qprod.coeff_fun, qprod.terms[1].expr[1], qprod.statespace)
end
function Cumulant(q::QAbstract, statespace::StateSpace, kwargs...)::QCumulant
    error("Cumulant not defined for qAbstract!")
end
function Cumulant(coeff_fun::FFunction, atom::QTerm, statespace::StateSpace)::QCumulant
    where_acting = where_acting_index(atom, statespace)
    order = length(where_acting)
    red_cum = ReducedIndexedCumulant(order)
    return _Cumulant(coeff_fun, atom, statespace, order, where_acting, red_cum)
end

function Cumulant(qprod::QAtomProduct, red_cum_list::ReducedCumulantList)::QCumulant
    if length(qprod.terms) != 1
        error("Cumulant only defined for single QAtom in QAtomProduct")
    end
    return Cumulant(qprod.coeff_fun, qprod.terms[1].expr[1], qprod.statespace, red_cum_list)
end
function Cumulant(coeff_fun::FFunction, atom::QTerm, statespace::StateSpace, red_cum_list::ReducedCumulantList)::QCumulant
    where_acting = where_acting_index(atom, statespace)
    order = length(where_acting)
    return _Cumulant(coeff_fun, atom, statespace, order, where_acting, red_cum_list(order))
end

@inline function _Cumulant(coeff_fun::FFunction, atom::QTerm, statespace::StateSpace, order::Int, where_acting::Vector{Int}, red_cum::ReducedIndexedCumulant)::QCumulant
    curr_op_indexes = atom.op_indices
    neutral_op = statespace.neutral_op
    # create indexed cumulant to this order 
    # create qExpr from this 
    qexpr::Vector{qComposite} = Vector{qComposite}(undef, length(red_cum.approximation))
    for (j, app) in enumerate(red_cum.approximation)
        # create qTerm from app, where acting and atom 
        coeff = app.coeff
        indexes = app.indices
        curr_atoms = Vector{QTerm}(undef, length(indexes))
        # replace indexes from neutral with tthose in qterm 
        for (i, ind) in enumerate(indexes)
            curr_atoms[i] = QTerm(statespace, replace_indexes(neutral_op, curr_op_indexes, ind))
        end
        qexpr[j] = QAtomProduct(statespace, coeff_fun*coeff, curr_atoms)
    end
    return QCumulant(statespace, copy(atom), qexpr, order, where_acting)
end

