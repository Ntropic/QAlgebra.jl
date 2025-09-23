"""
### Baker-Campbell-Hausdorff (BCH) Expansion
We use the following combinatorial representation of the BCH expansion to derive the terms algorithmically.
$$
Z=\log\left(e^X e^Y\right) = \sum_{n=1}^\infty \frac{(-1)^{n-1}}{n} \sum_{\substack{r_1 + s_1 > 0 \\ r_2 + s_2 > 0 \\ \vdots \\ r_n + s_n > 0}} \frac{[X^{r_1} Y^{s_1} X^{r_2} Y^{s_2} \dotsm X^{r_n} Y^{s_n}]}{\left(\sum_{j=1}^n (r_j + s_j)\right) \cdot \prod_{i=1}^n r_i! s_i!}
$$

as an infinite series involving nested commutators of $X$ and $Y$:

$$
Z = X + Y + \frac{1}{2}[X, Y] + \frac{1}{12}[X, [X, Y]] - \frac{1}{12}[Y, [X, Y]] - \frac{1}{24}[Y, [X, [X, Y]]] + \dots
$$

where we used the nested commutator notation $[ABC\dots] = [A, [B, [C, \dots]]]$, where we substitutesd repeated operators with their powers, e.g., $X^2 = X X$.
Each term in the series captures the non-commutative nature of the operators and involves increasingly nested commutators.

### Algorithm for Generating the BCH Terms

To generate the BCH terms up to a specified order $ N $, we perform the following steps:

1. **Loop over $ n $ from 1 to $ N $**:
   - For each $ n $, generate all sequences $ \{(r_i, s_i)\}_{i=1}^n $ such that $ r_i + s_i > 0 $ and $ r_i, s_i \geq 0 $.
2. **For each sequence**:
   - Calculate the total degree $ D = \sum_{j=1}^n (r_j + s_j) $.
   - Compute the coefficient:
     $$
     c = \frac{(-1)^{n-1}}{n \cdot D \cdot \prod_{i=1}^n r_i! s_i!}
     $$
   - Construct the nested commutator corresponding to the sequence $ X^{r_1} Y^{s_1} \dotsm X^{r_n} Y^{s_n} $.
3. **Collect terms of the same nested commutator**:
   - Combine coefficients for identical nested commutators.

This algorithm allows us to systematically generate all terms up to order $ N $ without predefining them.
"""

include("Sum.jl")


mutable struct Commutator
    subterms::Vector{Union{Int, Commutator}}
    max_ind::Int
    order::Int
    function Commutator(subterms::Vector)
        comm = new(Vector{Union{Int, Commutator}}(), 0, 0)
        if length(subterms) > 2
            error("Commutator must have at most 2 terms")
        end 
        for term in subterms
            if isa(term, Int)
                comm.max_ind = max(comm.max_ind, term)
                comm.order += 1
                push!(comm.subterms, term)
            elseif isa(term, Vector) # Make Commutator from Vector
                new_term = Commutator(term)
                push!(comm.subterms, new_term)
                comm.max_ind = max(comm.max_ind, new_term.max_ind)
                comm.order += new_term.order
            elseif isa(term, Commutator)
                push!(comm.subterms, term)
                comm.max_ind = max(comm.max_ind, term.max_ind)
                comm.order += term.order
            else
                error("Invalid term type, each term must be Int, Commutator, or Vector. ")
            end
        end
        return comm
    end
    function Commutator(A::Union{Int, Vector, Commutator}, B::Union{Int, Vector, Commutator})
        return Commutator([A, B])
    end
    function Commutator(subterms::Vector, max_ind::Int, order::Int)
        return new(subterms, max_ind, order)
    end
end
function Base.getindex(comm::Commutator, i::Int)
    return comm.subterms[i]
end
function Base.setindex!(comm::Commutator, value, i::Int)
    comm.subterms[i] = value
end
function Base.length(comm::Commutator)
    return length(comm.subterms)
end
# Define a custom show method for the Commutator type
function Base.show(io::IO, comm::Commutator)
    print(io, "Commutator: ", format(comm))
end
# Helper function to format the commutator recursively
function format(comm::Commutator)::String
    if length(comm) == 1
        return format(comm[1])
    else 
        return "[" * format(comm[1]) * "," * format(comm[2]) * "]"
    end
end
function format(element::Int)::String
    return string(element)
end
function Base.string(comm::Commutator)
    return format(comm)
end
# Test 
Commutator([1,[2,3]])


function Base.:+(A::Commutator, B::Commutator)
    return Sum([1//1, 1//1], [A, B])
end
function Base.:-(A::Commutator, B::Commutator)
    return Sum([1//1, -1//1], [A, B])
end
function max_ind(comm::Union{Int, Commutator})
    if isa(comm, Int)
        return comm
    else
        return comm.max_ind
    end
end
function order(comm::Union{Int, Commutator})
    if isa(comm, Int)
        return 1
    else
        return comm.order
    end
end
function max_ind_order(comm::Union{Int, Commutator})
    if isa(comm, Int)
        return comm, 0
    else
        return comm.max_ind, comm.order
    end
end
function Base.:<(A::Commutator, B::Int)
    return false
end
function Base.:<(A::Int, B::Commutator)
    return true
end
function Base.:<(A::Commutator, B::Commutator)
    max_ind_A, max_order_A = A.max_ind, A.order
    max_ind_B, max_order_B = B.max_ind, B.order
    if max_order_A < max_order_B
        return true
    elseif max_order_A > max_order_B
        return false
    else
        if max_ind_A < max_ind_B
            return true
        elseif max_ind_A > max_ind_B
            return false
        else
            # recursively check both subterms
            if A[1] < B[1]
                return true
            elseif A[1] > B[1]
                return false
            elseif A[2] < B[2]
                return true
            else
                return false
            end
        end
    end
end
function Base.isless(A::Commutator, B::Commutator)
    return A < B 
end
function Base.:>(A::Union{Int, Commutator}, B::Union{Int, Commutator})
    return B < A
end
function Base.:(==)(A::Commutator, B::Commutator)
    if A.max_ind != B.max_ind || A.order != B.order
        return false
    end
    return A.subterms == B.subterms
end
# Test
A = Commutator([1,[2,3]])
B = Commutator([1,[2,3]])
A == B


function flip_AB(C::Commutator)
    # Flip the order of the elements in the commutator
    if length(C) == 1
        return C
    else
        return Commutator([C[2], C[1]], C.max_ind, C.order)
    end
end
function A_smaller_B(C::Commutator)
    # Check if A < B
    if length(C) == 1
        error("Commutator must have 2 terms")
    else
        return C[1] < C[2]
    end
end
function A_equals_B(C::Commutator)
    # Check if A = B
    if length(C) == 1
        error("Commutator must have 2 terms")
    else
        return C[1] == C[2]
    end
end
# sign gets flipped when flipping A and B
function repartition_Commutator(comm_in::Commutator; sign::Int=1)::Tuple{Commutator, Int}
    # recursively reorganize to create nesting to the right (until every level is A<B within the stack)
    comm = deepcopy(comm_in)
    inner_sign = sign
    if comm.order == 0
        return comm, sign
    elseif comm.order == 1
        # deepest level reached, 
        if A_smaller_B(comm)
            return comm, sign
        else
            return flip_AB(comm), -inner_sign
        end
    else
        is_smaller = A_smaller_B(comm)
        if isa(comm[1], Commutator)
            comm[1], new_sign = repartition_Commutator(comm[1]; sign=sign)
            inner_sign *= new_sign
        end
        if isa(comm[2], Commutator)
            comm[2], new_sign = repartition_Commutator(comm[2]; sign=sign)
            inner_sign *= new_sign
        end
        if is_smaller
            return comm, inner_sign
        else
            return flip_AB(comm), -inner_sign
        end
    end
end

# Test 
comm = Commutator([[[4,2],[3,6]],[1,5]])
repartition_Commutator(comm)


# Function to check if A and B in a Commutator are equal
function is_zero(comm::Commutator)::Bool
    # Check if commutator comm is zero 
    # checks if any commutator in the stack is equal on both sides 
    if length(comm) == 1
        return false
    end
    if isa(comm[1], Commutator)
        # check inside of it for is_zero 
        if is_zero(comm[1])
            return true
        end
    end
    if isa(comm[2], Commutator)
        # check inside of it for is_zero 
        if is_zero(comm[2])
            return true
        end
    end
    if isa(comm[1], Int) && isa(comm[2], Int)
        return comm[1] == comm[2]
    elseif isa(comm[1], Commutator) && isa(comm[2], Commutator)
        return A_equals_B(comm)
    end
    return false
end


mutable struct OperatorProduct
    operators::Vector{Int}
    function OperatorProduct(operators::Vector{Int})
        return new(operators)
    end
end
function format(op::OperatorProduct)::String
    return join(string.(op.operators), "*")
end
function Base.show(io::IO, op::OperatorProduct)
    print(io, "OperatorProduct: ", format(op))
end
function Base.string(op::OperatorProduct)
    return format(op)
end
# Test 
op = OperatorProduct([1,2,3])
op


function Base.:*(A::OperatorProduct, B::OperatorProduct)
    return OperatorProduct([A.operators..., B.operators...])
end
function Base.:+(A::OperatorProduct, B::OperatorProduct)
    return Sum([1//1, 1//1], [A, B])
end
function Base.:-(A::OperatorProduct, B::OperatorProduct)
    return Sum([1//1, -1//1], [A, B])
end
function Base.isless(A::OperatorProduct, B::OperatorProduct)
    if length(A.operators) < length(B.operators)
        return true
    elseif length(A.operators) > length(B.operators)
        return false
    else
        return A.operators < B.operators
    end
end
function Base.:<(A::OperatorProduct, B::OperatorProduct)
    return isless(A, B)
end
function Base.:>(A::OperatorProduct, B::OperatorProduct)
    return isless(B, A)
end
function Base.:(==)(A::OperatorProduct, B::OperatorProduct)
    return A.operators == B.operators
end


function expand(comm::Commutator)::Sum{OperatorProduct}
    # Expand the commutator into an OperatorProduct
    if length(comm) == 1
        return Sum([OperatorProduct([comm[1]])])
    else # Recursively expand the commutator 
        sum_A = expand(comm[1])
        sum_B = expand(comm[2])
        return sum_A * sum_B - sum_B * sum_A
    end
end
function expand(s::Sum{Commutator})::Sum{OperatorProduct}
    # Expand the sum of commutators into a sum of OperatorProducts
    sumy = expand(s.elements[1]) * s.coeffs[1]
    for i in 2:length(s)
        sumy += expand(s.elements[i]) * s.coeffs[i]
    end
    return sumy
end
function expand(comm::Int)::Sum{OperatorProduct}
    return Sum([OperatorProduct([comm])])
end
# Test 
comm = Commutator([1,[2,3]])
sumy = expand(comm)


A = Commutator([1,[2,3]])
B = Commutator([2,[3,1]])
sumy = simplify(expand(A+B))

function vector_to_nested(v::Vector{Int})
    # Convert a vector of integers to a nested commutator
    if length(v) == 1
        return Commutator(v)
    else
        return Commutator([v[1], vector_to_nested(v[2:end])])
    end
end

function rs_to_indexes(r_seq, s_seq)
    # Convert r and s sequences to indexes
    indexes = []
    if length(r_seq) != length(s_seq)
        throw(ArgumentError("r_seq and s_seq must have the same length"))
    end
    for i in 1:length(r_seq)
        for j in 1:r_seq[i]
            push!(indexes, 1)
        end
        for j in 1:s_seq[i]
            push!(indexes, 2)
        end
    end
    return indexes
end
# Test 
rs_to_indexes([1,2,1], [1,1,2])
function generate_rs_sequences(m::Int; curr_rs::Vector{Int} = Int[], curr_ss::Vector{Int} = Int[], curr_sum::Int = 0)
    rs::Vector{Vector{Int}} = []
    ss::Vector{Vector{Int}} = []
    remaining = m - curr_sum 
    for curr_increase in 1:remaining-1
        # all combinations of two numbers summing up to curr_increase
        for r in 0:curr_increase
            s = curr_increase - r
            new_rs = [curr_rs..., r]
            new_ss = [curr_ss..., s]
            # recursively call generate_rs_sequences with the new values
            all_new_rs, all_new_ss = generate_rs_sequences(m; curr_rs=new_rs, curr_ss=new_ss, curr_sum=curr_sum+curr_increase)
            for (new_r, new_s) in zip(all_new_rs, all_new_ss)
                push!(rs, new_r)
                push!(ss, new_s)
            end
        end
    end
    for r in 0:remaining
        s = remaining - r
        new_rs = [curr_rs..., r]
        new_ss = [curr_ss..., s]
        push!(rs, new_rs)
        push!(ss, new_ss)
    end
    return rs, ss
end
function factorial_array(n::Int)::Vector{Int}
    fact = 1
    fact_array::Vector{Int} = [1]
    for i in 1:n
        fact *= i
        push!(fact_array, fact)
    end
    return fact_array
end

# Generate BCH expansion for a given order N 
function construct_BCH_terms(N::Int)::Tuple{Vector{Vector{Vector{Int}}}, Vector{Vector{Rational{Int}}}}
    fact::Vector{Int} = factorial_array(N)
    prefactors::Vector{Rational{Int}} = [(-1)^(n-1)//n for n in 1:N]
    indexes::Vector{Vector{Vector{Int}}} = []   # sorted by order
    coeffs::Vector{Vector{Rational{Int}}} = []  # sorted by order
    for m in 1:N
        rs, ss = generate_rs_sequences(m)
        curr_indexes::Vector{Vector{Int}} = []
        curr_coeffs::Vector{Rational{Int}} = []
        for (r_seq, s_seq) in zip(rs, ss)
            inds = rs_to_indexes(r_seq, s_seq)
            factorials_rs = prod([fact[r+1] for r in r_seq])*prod([fact[s+1] for s in s_seq])
            push!(curr_indexes, inds)
            push!(curr_coeffs, prefactors[length(r_seq)]//factorials_rs//m )
        end
        push!(indexes, curr_indexes)
        push!(coeffs, curr_coeffs)
    end
    return indexes, coeffs
end
function simplify_BCH_terms(indexes::Vector{Vector{Vector{Int}}}, coeffs::Vector{Vector{Rational{Int}}})::Tuple{Vector{Vector{Int}}, Vector{Rational{Int}}}
    # first sort terms and coeffs accordingly #
    reduced_inds::Vector{Vector{Int}} = []
    reduced_coeffs::Vector{Rational{Int}} = []
    for (inds, c) in zip(indexes, coeffs)
        if length(inds) == 0
            continue
        end
        order = sortperm(inds)   # sorted by operators
        sorted_inds = inds[order]
        sorted_coeffs = c[order]
        # check iteratively if term is repeated 
        new_inds::Vector{Vector{Int}} = []
        new_coeffs::Vector{Rational{Int}} = []
        len_2 = length(sorted_inds[1]) > 1

        curr_inds = sorted_inds[1]
        curr_coeff = sorted_coeffs[1]
        if len_2 
            if curr_inds[end-1] > curr_inds[end] # check if last operators are sorted 
                curr_inds[end-1], curr_inds[end] = curr_inds[end], curr_inds[end-1]
                sorted_coeffs[1] = -sorted_coeffs[1]
            end
        end
        for i in 2:length(sorted_inds)
            inds = sorted_inds[i]
            coeff = sorted_coeffs[i]
            if len_2 
                if inds[end-1] > inds[end] # check if last operators are sorted 
                    inds[end-1], inds[end] = inds[end], inds[end-1]
                    coeff = -coeff
                end
            end
            if inds == curr_inds
                curr_coeff += coeff
            else
                if len_2
                    if curr_inds[end-1] == curr_inds[end]
                        curr_coeff = 0
                    end
                end
                if curr_coeff != 0
                    push!(new_inds, curr_inds)
                    push!(new_coeffs, curr_coeff)
                end
                curr_inds = sorted_inds[i]
                curr_coeff = sorted_coeffs[i]
            end
        end
        append!(reduced_inds, new_inds)
        append!(reduced_coeffs, new_coeffs)
    end
    return reduced_inds, reduced_coeffs
end
function BCH(N::Int)
    inds, coeffs = construct_BCH_terms(N)
    reduced_inds, reduced_coeffs = simplify_BCH_terms(inds, coeffs)
    # Construct a sum of commutators 
    # use vector_to_nested to convert the indexes to a nested commutator
    return Sum(reduced_coeffs, Commutator[vector_to_nested(inds) for inds in reduced_inds])
end
# Test 
bch = BCH(5)