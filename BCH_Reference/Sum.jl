using LaTeXStrings

mutable struct Sum{T, V}
    coeffs::Vector{V}
    elements::Vector{T}
    function Sum(coeffs::Vector{V}, elements::Vector{T}) where {T, V}
        return new{T, V}(coeffs, elements)
    end

    function Sum(elements::Vector{T}) where {T}
        return new{T, Rational{Int}}(ones(Rational{Int}, length(elements)), elements)
    end
end

Base.length(cs::Sum) = length(cs.elements)

function Base.push!(s::Sum{T, V}, coeff::V, elem::T) where {T, V}
    push!(s.elements, elem)
    push!(s.coeffs, coeff)
end

function Base.deleteat!(s::Sum{T, V}, i::Int) where {T, V}
    deleteat!(s.elements, i)
    deleteat!(s.coeffs, i)
end

function Base.deleteat!(s::Sum{T, V}, is::Vector{Int}) where {T, V}
    deleteat!(s.elements, is)
    deleteat!(s.coeffs, is)
end

function Base.getindex(s::Sum{T, V}, i::Int) where {T, V}
    return s.coeffs[i], s.elements[i]
end

function Base.iterate(s::Sum{T, V}, state::Int = 1) where {T, V}
    if state > length(s.elements)
        return nothing  # End of iteration
    end
    return (s.coeffs[state], s.elements[state]), state + 1
end
function rational2str(r::Rational{Int}; latex::Bool=false)::String
    if r.den == 1
        return string(r.num)
    else
        if latex
            curr_str = "\\frac{" * string(abs(r.num)) * "}{" * string(r.den) * "}"
        else
            curr_str = string(abs(r.num)) * "/" * string(r.den)
        end
        if r.num < 0
            curr_str = "-" * curr_str
        end
        return curr_str
    end
end
function Base.string(rat::Rational{Int}; latex::Bool=false)::String
    return rational2str(rat)
end

function do_string(obj; latex::Bool=false)::String
    if !latex 
        return string(obj)
    else 
        try 
            return string(obj, latex=true)
        catch 
            return string(obj)
        end
    end
end
function Base.string(s::Sum{T, V}; latex::Bool=false)::String where {T, V}
    stringer = ""
    if latex 
        stringer *= "\$"
    end
    if length(s) == 0
        stringer *= "0"
    end
    ind = 0
    for (coeff, elem) in s
        # if hasmethod(string, Tuple{V, Bool})
        string_elem = do_string(elem, latex=latex)
        string_coeff = do_string(coeff, latex=latex)
        ind += 1
        if occursin("+", string_coeff[2:end]) || occursin("-", string_coeff[2:end]) 
            if ind > 1
                stringer *=" +"
            end
            stringer *= "(" * strip(string_coeff) * ")"
        else
            if string_coeff[1] == '+' 
                if ind > 1
                    stringer *=" +"
                end
                string_coeff = string_coeff[2:end] 
            elseif string_coeff[1] == '-' 
                stringer *=" -"
                string_coeff = string_coeff[2:end] 
            else
                if ind > 1
                    stringer *= " +"
                end
                if string_coeff[1] == ' '
                    string_coeff = string_coeff[2:end]
                end
            end
            if !(string_coeff=="1" && length(string_coeff) > 0)
                stringer *= string_coeff
                if latex 
                    stringer *= "\\cdot"
                else
                    stringer *= " "
                end
            end
        end
        stringer *= string_elem
    end
    if latex 
        stringer *= " \$"
    end
    return stringer 
end
function Base.show(io::IO, ::MIME"text/plain", s::Sum{T, V}) where {T, V}
    print(io, string(s))
end
function Base.show(io::IO, ::MIME"text/latex", s::Sum{T, V}) where {T, V}
    print(io, latexstring(string(s, latex=true)))
end


# Now the important math functions
function Base.:*(A::Sum{T, V}, B::Sum{T, V}) where {T, V}
    coeffs = V[]
    elements = T[]
    for (coeff_A, elem_A) in A
        for (coeff_B, elem_B) in B
            push!(coeffs, coeff_A * coeff_B)
            push!(elements, elem_A * elem_B)
        end
    end
    return Sum(coeffs, elements)
end
function Base.:*(coeff::V, A::Sum{T, V}) where {T, V}
    return Sum(coeff * A.coeffs, A.elements)
end
function Base.:*(A::Sum{T, V}, coeff::V) where {T, V}
    return Sum([coeff * A_coeff for A_coeff in A.coeffs], A.elements)
end

function Base.:+(A::Sum{T, V}, B::Sum{T, V}) where {T, V}
    return simplify(Sum([A.coeffs..., B.coeffs...], [A.elements..., B.elements...]))
end
function Base.:+(A::Sum{T, V}, term::T)::Sum{T,V} where {T, V}
    return simplify(Sum([A.coeffs..., V(1)], [A.elements..., term]))
end
function Base.:+(term::T, B::Sum{T, V})::Sum{T,V} where {T,V}
    return simplify(Sum([V(1), B.coeffs...], [term, B.elements...]))
end
function Base.:-(A::Sum{T, V}, B::Sum{T, V}) where {T, V}
    new_coeffs = deepcopy(A.coeffs)
    for coeff in B.coeffs
        push!(new_coeffs, -coeff)
    end
    return simplify(Sum(new_coeffs, [A.elements..., B.elements...]))
end

function Base.sort(S::Sum{T, V}) where {T, V}
    sorted = sortperm(S.elements)
    return Sum(S.coeffs[sorted], S.elements[sorted])
end

function Base.sort!(S::Sum{T, V}) where {T, V}
    sorted = sortperm(S.elements)
    S.coeffs = S.coeffs[sorted]
    S.elements = S.elements[sorted]
end

function simplify(S::Sum{T, V}) where {T, V}
    sort!(S)
    i = 1
    while i < length(S.elements)
        if S.coeffs[i] == 0
            deleteat!(S, i)
        elseif S.elements[i] == S.elements[i + 1]
            S.coeffs[i] += S.coeffs[i + 1]
            deleteat!(S, i + 1)
            if S.coeffs[i] == 0
                deleteat!(S, i)
            end
        else
            i += 1
        end
    end
    return S
end

