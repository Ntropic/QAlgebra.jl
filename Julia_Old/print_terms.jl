# A set of functions for printing terms in the operator basis
include("operator_terms.jl")
using LaTeXStrings
using PrettyTables
using Printf


#### Fractions ###############################################################
function closest_fraction_complex(x, tol::Float64=1e-12)Tuple{Int,Int,Int,Int}
    # find the closest fraction to a potentially complex number
    # x: number to be approximated
    # tol: tolerance for the approximation
    # if has complex part then return two Fractions
    # else return one Fraction
    if x isa Int
        x = Complex(x)
    end
    num_real::Int, den_real::Int = 0, 0
    num_imag::Int, den_imag::Int = 0, 0
    if abs(real(x)) > tol
        frac = rationalize(real(x), tol=tol)
        num_real = numerator(frac)
        den_real = denominator(frac)
    end
    if abs(imag(x)) > tol
        frac = rationalize(imag(x), tol=tol)
        num_imag = numerator(frac)
        den_imag = denominator(frac)
    end
    return num_real, den_real, num_imag, den_imag
end
#closest_fraction_complex(1)

function print_int(num::Int; do_imag::Bool=false, var_str::String="")::String
    sign_str::String = ""
    if num < 0
        sign_str = "-"
        num = -num
    end
    imag_str::String = ""
    if do_imag
        imag_str = "i"
    end
    num_str::String = ""
    if num == 1
        num_str = ""
        if length(var_str) == length(imag_str) == 0
            num_str = "1"
        end
    elseif num == 0
        num_str = ""
        if length(var_str) == length(imag_str) == 0
            num_str = "0"
        end
    else
        num_str = string(num)
    end
    return sign_str * num_str * imag_str * var_str
end
## Test 
#vals = [-1, 0, 1, 2]
#var_strs = ["", "x"]
#do_imags = [false, true]
#for val in vals
#    for var_str in var_strs
#        for do_imag in do_imags
#            print(print_int(val, do_imag=do_imag, var_str=var_str))
#            print("  |  ")
#        end
#    end
#end

function frac_to_str(num::Int, den::Int; do_latex::Bool=false, do_imag::Bool=false, var_str::String="")::String
    # uses print_int to transform numerator and denominator to strings (num, den)
    # do_latex: if true then return LaTeX string \frac{num}{den} else return string (num/den)
    # var_str is part of num, sign goes before \frac
    if den == 0
        return "0"
    elseif den == 1 # no frac
        return print_int(num, do_imag=do_imag, var_str=var_str)
    else # do frac
        sign_str::String = ""
        if num < 0
            sign_str = "-"
            num = -num
        end
        num_str::String = ""
        if do_latex
            num_str = print_int(num, do_imag=do_imag, var_str=var_str)
            return sign_str * "\\frac{" * num_str * "}{" * string(den) * "}"
        else
            imag_str::String = ""
            num_str = print_int(num, do_imag=do_imag, var_str=var_str)
            return sign_str * num_str * "/" * string(den)
        end
    end
end
## Test
#vals = [1/2, 1, 2, 3/2]
#var_strs = ["", "x"]
#do_imags = [false, true]
#for val in vals
#    for var_str in var_strs
#        for do_imag in do_imags
#            num, den = closest_fraction_complex(val)
#            display(latexstring(frac_to_str(num, den, do_latex=true, do_imag=do_imag, var_str=var_str)))
#            display(frac_to_str(num, den, do_latex=false, do_imag=do_imag, var_str=var_str))
#        end
#    end
#end

function complex_to_fraction_to_str(value; do_latex::Bool=false, var_str::String="")::String
    if !(value isa ComplexF64)
        value = ComplexF64(value)
    end
    num_real, den_real, num_imag, den_imag = closest_fraction_complex(value)
    is_real::Bool = false
    is_imag::Bool = false
    is_both::Bool = false
    if num_real != 0
        is_real = true
    end
    if num_imag != 0
        is_imag = true
    end
    if is_real && is_imag
        is_both = true
    end
    frac_str::String = ""
    if is_both
        frac_str *= frac_to_str(num_real, den_real, do_latex=do_latex, do_imag=false, var_str="") * " + "
        frac_str *= frac_to_str(num_imag, den_imag, do_latex=do_latex, do_imag=true, var_str="")
        if do_latex
            if length(var_str) > 0
                frac_str *= "\\left(" * frac_str * "\\right)" * var_str
            end
        else
            if length(var_str) > 0
                frac_str *= "(" * frac_str * ")" * var_str
            end
        end
    elseif is_real
        frac_str *= frac_to_str(num_real, den_real, do_latex=do_latex, do_imag=false, var_str=var_str)
    elseif is_imag
        frac_str *= frac_to_str(num_imag, den_imag, do_latex=do_latex, do_imag=true, var_str=var_str)
    else
        frac_str = "0"
    end
    return frac_str
end

function varstring_from_var_array(vars::Vector{String}; do_latex::Bool=true)::String
    # convert a vector of variables to a string
    # vars: vector of variables
    # do_latex: if true then return a latex string
    if !do_latex  # substitutions i.e. "\beta" -> "β"
        var_substitutions = Dict(raw"\alpha" => "α", raw"^{*}" => "*", raw"\beta" => "β", raw"\gamma" => "γ", raw"\delta" => "δ", raw"\kappa" => "κ", raw"\lambda" => "λ",
            raw"\mu" => "μ", raw"\nu" => "ν", raw"\rho" => "ρ", raw"\sigma" => "σ", "\tau" => "τ", raw"\phi" => "φ", raw"\chi" => "χ", raw"\psi" => "ψ", raw"\omega" => "ω",
            raw"\epsilon" => "ε", raw"\zeta" => "ζ", raw"\eta" => "η", raw"\theta" => "θ", raw"\iota" => "ι", raw"\xi" => "ξ", raw"\pi" => "π", raw"\upsilon" => "υ", raw"\omega" => "ω",
            raw"\Gamma" => "Γ", raw"\Delta" => "Δ", raw"\Theta" => "Θ", raw"\Lambda" => "Λ", raw"\Xi" => "Ξ", raw"\Pi" => "Π", raw"\Sigma" => "Σ", raw"\Upsilon" => "Υ", raw"\Phi" => "Φ",
            raw"\Psi" => "Ψ", raw"\Omega" => "Ω", raw"\sqrt{\kappa}" => "√̅κ")
    end
    var_str::String = ""
    for i in 1:length(vars)
        curr_var = vars[i]
        if do_latex
            var_str *= curr_var * " "
        else
            # Substitute certain variables with substitution table (for ascii) -> substitute for every element in the substitution Dict
            for (key, value) in var_substitutions
                curr_var = replace(curr_var, key => value)
            end
            var_str *= curr_var
        end
    end
    return var_str #*" "
end

function term2str(term::Term; do_sigma::Bool=false, do_latex::Bool=true, do_braket::Bool=false, do_vars::Bool=true)::String
    # create a string (for latex or printing) of a term
    spin_types::String = term.spin_types
    spin_indices::String = term.spin_indices
    bosons::String = term.bosons
    exponents::Vector{Int} = term.exponents
    coeff::ComplexF64 = term.coeff
    vars::Vector{String} = term.vars
    is_pauli::Bool = term.is_pauli

    boson_substitutions::Dict{Char,String} = Dict()
    spin_substitutions::Dict{Char,String} = Dict()
    if do_latex
        boson_substitutions = Dict('+' => raw"\hat{a}^\dagger", '-' => raw"\hat{a}")
        if do_sigma
            if is_pauli
                spin_substitutions = Dict('x' => raw"\hat{\sigma}_x", 'y' => raw"\hat{\sigma}_y", 'z' => raw"\hat{\sigma}_z")
            else
                spin_substitutions = Dict('m' => raw"\hat{\sigma}_{-}", 'p' => raw"\hat{\sigma}_{+}", 'z' => raw"\hat{\sigma}_z")
            end
        else
            if is_pauli
                spin_substitutions = Dict('x' => raw"\hat{x}", 'y' => raw"\hat{y}", 'z' => raw"\hat{z}")
            else
                spin_substitutions = Dict('m' => raw"\hat{m}", 'p' => raw"\hat{p}", 'z' => raw"\hat{z}")
            end
        end
    else # substitutions i.e. "\beta" -> "β"
        boson_substitutions = Dict('+' => raw"a†", '-' => raw"a")
        superscript_indexes = Dict('a' => "ᵃ", 'b' => "ᵇ", 'c' => "ᶜ", 'd' => "ᵈ", 'e' => "ᵉ", 'f' => "ᶠ",
            'g' => "ᵍ", 'h' => "ʰ", 'i' => "ⁱ", 'j' => "ʲ", 'k' => "ᵏ", 'l' => "ˡ", 'm' => "ᵐ", 'n' => "ⁿ",
            'o' => "ᵒ", 'p' => "ᵖ", 'q' => "ᵠ", 'r' => "ʳ", 's' => "ˢ", 't' => "ᵗ", 'u' => "ᵘ", 'v' => "ᵛ",
            'w' => "ʷ", 'x' => "ˣ", 'y' => "ʸ", 'z' => "ᶻ")
        subscript_indexes = Dict('h' => "ₕ", 'i' => "ᵢ", 'j' => "ⱼ", 'k' => "ₖ", 'l' => "ₗ", 'm' => "ₘ", 'n' => "ₙ", 'o' => "ₒ", 'p' => "ₚ", '1' => "₁", '2' => "₂", '3' => "₃", '4' => "₄", '5' => "₅", '6' => "₆", '7' => "₇", '8' => "₈", '9' => "₉")
        superscript_numbers = Dict(1 => "¹", 2 => "²", 3 => "³", 4 => "⁴", 5 => "⁵", 6 => "⁶", 7 => "⁷", 8 => "⁸", 9 => "⁹", 10 => "¹⁰", 11 => "¹¹", 12 => "¹²", 13 => "¹³", 14 => "¹⁴", 15 => "¹⁵", 16 => "¹⁶", 17 => "¹⁷", 18 => "¹⁸", 19 => "¹⁹", 20 => "²⁰")
    end
    # Do in parts, 1st the coeff part, 2nd the vars part, 3rd the boson part and 4th the spin part
    # Coeff part
    coeff_str::String = ""
    # Vars part
    vars_str::String = ""
    if do_vars
        vars_str = varstring_from_var_array(vars, do_latex=do_latex)
        coeff_str = complex_to_fraction_to_str(coeff, do_latex=do_latex, var_str=vars_str)
    else
        coeff_str = complex_to_fraction_to_str(coeff, do_latex=do_latex, var_str="")
    end
    # Boson Part (cavity part)
    boson_str::String = ""
    # iterate over bosons and exponents
    for (boson, exponent) in zip(bosons, exponents)
        if exponent == 1
            boson_str *= boson_substitutions[boson]
        elseif exponent > 1
            if do_latex
                if boson == '-'
                    boson_str *= boson_substitutions[boson] * "^{$exponent} "
                else
                    boson_str *= "\\hat{a}^{\\dagger $exponent} "
                end
            else
                boson_str *= boson_substitutions[boson] * superscript_numbers[exponent]
            end
        end
    end
    # Spin Part
    spin_str::String = ""
    # iterate over spin types and spin indices
    for i in 1:length(spin_types)
        spin_type = spin_types[i]
        spin_index = spin_indices[i]
        if do_latex
            spin_str *= spin_substitutions[spin_type]
            if !(spin_index == ' ')
                if do_sigma
                    spin_str *= "^{(" * spin_index * ")} "
                else
                    spin_str *= "_{" * spin_index * "} "
                end
            end
        else
            spin_str *= string(spin_type)
            if !(spin_index == ' ')
                spin_str *= subscript_indexes[spin_index]    # Old way: "⁽"*superscript_indexes[spin_index]*"⁾"
            end
        end
    end

    # Put it all together	
    if do_latex
        if length(boson_str) == length(spin_str) == 0
            return coeff_str
        else
            if coeff_str == "1"
                coeff_str = ""
            elseif coeff_str == "-1"
                coeff_str = "-"
            end
            if do_braket
                spin_str = coeff_str * " \\braket{" * boson_str * " " * spin_str * "}"
            else
                spin_str = coeff_str * " \\phantom{,}" * boson_str * " " * spin_str * "\\phantom{,}"
            end
        end
    else
        if length(boson_str) == length(spin_str) == 0
            spin_str = coeff_str
        else
            if coeff_str == "1"
                coeff_str = ""
            elseif coeff_str == "-1"
                coeff_str = "-"
            end
            if do_braket
                #return coeff_str * "〈" * boson_str * spin_str * "〉"
                spin_str = coeff_str * "⟨" * boson_str * spin_str * "⟩"
            else
                spin_str = coeff_str * " " * boson_str * spin_str * " "
            end
        end
    end
    if term.conjugate
        if do_latex
            spin_str *= "^{*}"
        else
            spin_str *= "*"
        end
    end
    return spin_str
end
## Test
#term = make_term("++--+-xiyjzk", 2.0, ["\\kappa", "\\alpha"])
#clean_up_bosons!(term)
#display(latexstring(term2str(term, do_latex=true, do_braket=true)))
#display(term2str(term, do_latex=false))
#
#term = make_term("++-x", .5im, ["A", "B"])
#clean_up_bosons!(term)
#display(latexstring(term2str(term, do_latex=true, do_braket=true)))
#display(term2str(term, do_latex=false))
#
#term = make_term("", 1.0, ["A", "B"])
#display(latexstring(term2str(term, do_latex=true, do_braket=true)))
#display(term2str(term, do_latex=false, do_braket=true))
#
#term = make_term("", 1.0)
#display(latexstring(term2str(term, do_latex=true, do_braket=true)))
#display(term2str(term, do_latex=false, do_braket=true))


############################################################################################################
##### Term Group Printing ##################################################################################
############################################################################################################

function same_vars(termA::Term, termB::Term)::Bool
    # Checks if two terms have the same variables
    # termA: first term
    # termB: second term
    # return: true if they have the same variables, false otherwise
    varsA::Vector{String} = termA.vars
    varsB::Vector{String} = termB.vars
    if length(varsA) != length(varsB)
        return false
    else
        # First sort the terms alphabetically
        sort!(varsA)
        sort!(varsB)
        for i in 1:length(varsA)
            if varsA[i] != varsB[i]
                return false
            end
        end
        return true
    end
end
## Test
#termA = make_term("++--+-xiyjzk", 2.0, ["A", "B"])
#termB = make_term("++-xiyj", 2.0, ["B", "A"])
#println(same_vars(termA, termB))    # should be true
#termC = make_term("++--+-xiyjzk", 2.0, ["A"])
#println(same_vars(termA, termC))    # should be false

function group_terms_by_vars(term_list::Vector{Term})::Dict{String,Vector{Term}}
    # group terms by their variables
    # term_list: list of terms
    # return: Dict of terms grouped by their variables
    # check for every term, if other terms have the same variables
    # if so, add them to the same group
    # if not, create a new group
    term_groups::Dict{String,Vector{Term}} = Dict()
    index_array::Vector{Int} = collect(1:length(term_list)) # remove indexes that are already in a group
    curr_term::Term = Term()
    other_term::Term = Term()
    curr_vars::Vector{String} = String[]
    curr_group::Vector{Term} = Term[]
    curr_indexes::Vector{Int} = Int[]
    while length(index_array) > 0
        curr_term = term_list[index_array[1]]
        curr_vars = curr_term.vars
        sort!(curr_vars)
        curr_group = [curr_term]
        curr_indexes = [1]
        for i in 2:length(index_array)
            other_term = term_list[index_array[i]]
            if same_vars(curr_term, other_term)
                push!(curr_group, other_term)
                push!(curr_indexes, i)
            end
        end
        # remove the indexes from index_array via curr_indexes
        index_array = deleteat!(index_array, curr_indexes)
        # add the group to the Dict
        curr_str_vars = join(curr_vars, "_")
        term_groups[curr_str_vars] = curr_group
    end
    return term_groups
end
## Test
#term_list = [make_term("+-x", 1, ["A", "B"]), make_term("+-y", 2, ["A"]), make_term("--z", 3, ["A", "B"])]
#term_groups = group_terms_by_vars(term_list)

function common_coeffs(coeff_list::Vector; tol::Float64=1e-12)::Tuple{Bool,ComplexF64}
    # find a possible common coeff of all terms
    # coeff_list: list of coefficients
    # return: do common denominator, common denominator

    # if any coeff_list element is close to zero -> return 1
    if any(abs.(coeff_list) .< tol)
        return false, ComplexF64(1)
    end

    # make sure all elements of coeff_list are of same datatype ComplexF64
    for i in 1:length(coeff_list)
        if !(coeff_list[i] isa ComplexF64)
            coeff_list[i] = ComplexF64(coeff_list[i])
        end
    end
    # All real or all complex?
    any_both::Bool = false
    all_imag::Bool = true
    new_coeff_list::Vector{Float64} = Float64[]
    for i in 1:length(coeff_list)
        if abs(imag(coeff_list[i])) > tol
            push!(new_coeff_list, imag(coeff_list[i]))
        else
            push!(new_coeff_list, real(coeff_list[i]))
            all_imag = false
        end
        if abs(imag(coeff_list[i])) > tol && abs(real(coeff_list[i])) > tol   # both 
            return false, ComplexF64(1)
        end
    end

    # are all multiples of each other?
    # smallest element
    smallest::Float64 = minimum(abs.(new_coeff_list))
    # check for all elements if they are multiples of smallest
    factor::Float64 = 0
    for i in 1:length(new_coeff_list)
        factor = new_coeff_list[i] / smallest
        if abs(factor - round(factor)) > tol
            return false, ComplexF64(1)
        end
    end
    # if all negative, then make smallest positive
    if maximum(new_coeff_list) < 0
        smallest = -smallest
    end
    # if all multiples of each other, then return true and smallest
    if all_imag
        return true, ComplexF64(0, smallest)
    else
        return true, ComplexF64(smallest)
    end
end
## Term
#coeff_list = [2, 2, 4.0]
#common_coeffs(coeff_list)

function common_coeff_terms(term_list::Vector{Term}; tol::Float64=1e-12)::Tuple{Bool,ComplexF64,Vector{Term}}
    # find a possible common coeff of all terms
    # term_list: list of terms
    # return: do_coeffs, coefficient, list of terms with common coeff
    # get the coefficients
    coeff_list::Vector{ComplexF64} = [term.coeff for term in term_list]
    # find the common coeff
    do_common, common_coeff = common_coeffs(coeff_list, tol=tol)
    # if no common coeff, then return the original list
    new_term_list::Vector{Term} = deepcopy(term_list)
    if !do_common
        return false, ComplexF64(1), new_term_list
    end
    # if common coeff, then divide all terms by common coeff as it is outside the sum
    for i in 1:length(new_term_list)
        new_term_list[i].coeff /= common_coeff
    end
    return true, common_coeff, new_term_list
end
## Test
#term_list = [make_term("+-x", 2, ["A", "B"]), make_term("+-y", 2, ["A"]), make_term("--z", 4, ["A", "B"])]
#do_coeff, coeff, terms = common_coeff_terms(term_list)
#println(coeff)
#for term in terms
#    println(term2str(term, do_latex=false))
#end

function terms2str(term_array::Vector{Term}; group_terms::Bool=true, do_sigma::Bool=false, do_latex::Bool=true, do_braket::Bool=false)::String
    # Creates a string from a list of terms
    # term_array: Vector of terms
    # group_terms: if true, group terms if they have the same variables (term.vars)
    # do_sigma (default false): if true, use sigma_{x,y,z} operators instead of x,y,z
    # do_latex (default true): if true, use latex symbols formatting, else use ascii
    # do_braket (default false): if true, use \braket notation, else use angle brackets (only matters for latex)
    # return: string of terms
    all_terms_str::String = ""
    if group_terms
        term_groups::Dict{String,Vector{Term}} = group_terms_by_vars(term_array)
        term_str::String = ""
        i = 0
        for (key, value) in term_groups
            i += 1
            if length(value) == 1
                term_str = term2str(value[1], do_sigma=do_sigma, do_latex=do_latex, do_braket=do_braket)
            else
                # find a possible common coeff of all terms
                term_str = ""
                do_common, common_coeff, new_term_list = common_coeff_terms(value)
                curr_vars = value[1].vars
                var_str = varstring_from_var_array(curr_vars, do_latex=do_latex)
                if do_common
                    term_str *= complex_to_fraction_to_str(common_coeff, do_latex=do_latex, var_str=var_str)
                    if term_str == "1"
                        term_str = ""
                    elseif term_str == "-1"
                        term_str = "-"
                    end
                else
                    term_str *= var_str
                end
                if do_latex
                    term_str *= "\\left("
                else
                    term_str *= "("
                end
                for j in 1:length(new_term_list)
                    term = new_term_list[j]
                    curr_term_str = term2str(term, do_sigma=do_sigma, do_latex=do_latex, do_braket=do_braket, do_vars=false)

                    curr_term_str = strip(curr_term_str)
                    if curr_term_str[1] == '-'
                        term_str *= " - " * curr_term_str[2:end]
                    elseif j > 1
                        term_str *= " + " * curr_term_str
                    else
                        term_str *= curr_term_str
                    end
                end
                if do_latex
                    term_str *= "\\right)"
                else
                    term_str *= ")"
                end
            end
            term_str = strip(term_str)
            if term_str[1] == '-'
                all_terms_str *= " - " * term_str[2:end]
            elseif i > 1
                all_terms_str *= " + " * term_str
            else
                all_terms_str *= term_str
            end
        end
    else # Don't group terms
        curr_term::Term = term_array[1]
        all_terms_str = term2str(curr_term, do_sigma=do_sigma, do_latex=do_latex, do_braket=do_braket)
        for i in 2:length(term_array)
            curr_term = term_array[i]
            curr_str = term2str(curr_term, do_sigma=do_sigma, do_latex=do_latex, do_braket=do_braket)
            curr_str = strip(curr_str)
            if curr_str[1] == '-'
                all_terms_str *= " - " * curr_str[2:end]
            else
                all_terms_str *= " + " * curr_str
            end
        end
    end
    return all_terms_str
end
## Test
#term_list = [make_term("+-x", 2.0, ["A", "B"]), make_term("+-y", 4.0, ["A", "B"]), make_term("--z", 2.0, ["C"])]
#display(latexstring(terms2str(term_list, do_latex=true, do_braket=true)))


## Print Terms
# For console display
function Base.show(io::IO, ::MIME"text/plain", x::Term)
    str = term2str(x, do_sigma=false, do_latex=false, do_braket=false)
    print(io, str)
end
# For Jupyter Notebook display as LaTeX
function Base.show(io::IO, ::MIME"text/latex", x::Term)
    latex_code = latexstring(term2str(x, do_sigma=false, do_latex=true, do_braket=false))
    print(io, latex_code)
end
## Test
#Op = make_term("z")
#display(Op)


function sumstring(ind::String="i", not_ind::Union{String,Vector{String}}=""; do_latex::Bool=true)::String
    # create a sum string
    # ind: index of sum
    # not_ind: what is not summed over
    # do_latex: if true then return a latex string
    if isa(not_ind, String)
        if length(not_ind) == 0
            not_ind = []
        else
            not_ind = [not_ind]
        end
    end
    do_not_ind::Bool = false
    if length(not_ind) > 0
        do_not_ind = true
    end
    if do_latex
        if do_not_ind
            not_ind_str::String = ""
            for i in 1:length(not_ind)
                not_ind_str *= not_ind[i]
            end
            return "\\sum_{" * ind * " \\neq " * not_ind_str * "}"
        else
            return "\\sum_{" * ind * "}"
        end
    else
        sum_str::String = ""
        neq_symbol::String = "̡̡₌"  # "̷̷₌"
        subscript_indexes::Dict{String,String} = Dict("h" => "ₕ", "i" => "ᵢ", "j" => "ⱼ", "k" => "ₖ", "l" => "ₗ", "m" => "ₘ", "n" => "ₙ", "o" => "ₒ", "p" => "ₚ", "s" => "ₛ", "t" => "ₜ")
        sum_str = "∑"
        if length(ind) == 1
            sum_str *= subscript_indexes[ind]
        end
        if do_not_ind
            #"∑["*ind*" ≠ "*not_ind*"] "
            sum_str *= neq_symbol
            for i in 1:length(not_ind)
                sum_str *= subscript_indexes[not_ind[i]]
            end
        end
        return sum_str
    end
end
## Test
#sumstring("i", ["j"], do_latex=false)

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

mutable struct DE_Term              # Differential Equation Term (for a Single Qubit)
    exp_op::Term                        # Expectation value, for which we derived the time derivative
    which_ind::Vector{String}           # Indexes contained in exp_op
    terms::Vector{Term}                 # Terms
    #cumulant_terms::Vector{Term}        # Terms that are cumulants (rethink this solution)
    # Inner constructor with default values
    function DE_Term(exp_op::Term, which_ind::Vector{String}; terms::Vector{Term}=Term[])#, cumulant_terms::Vector{Term}=Term[])
        new(exp_op, which_ind, terms)#, cumulant_terms)
    end
end

function terms2str(x::DE_Term; group_terms::Bool=true, do_sigma::Bool=false, do_latex::Bool=true, do_braket::Bool=true)::String
    exp_op::Term = x.exp_op
    terms::Vector{Term} = x.terms
    # create the string
    eq_str::String = ""
    exp_op_str::String = term2str(exp_op, do_sigma=do_sigma, do_latex=do_latex, do_braket=do_braket)
    terms_str::String = terms2str(terms, group_terms=group_terms, do_sigma=do_sigma, do_latex=do_latex, do_braket=do_braket)
    if do_latex
        eq_str = "\\frac{\\mathrm{d}}{\\mathrm{d}t}" * exp_op_str * "= " * terms_str
    else
        eq_str = "d/dt" * exp_op_str * " = " * terms_str
    end
    return eq_str
end

function Base.show(io::IO, ::MIME"text/plain", x::DE_Term)
    latex_string = terms2str(x, do_latex=false, do_braket=true)
    print(io, latex_string)
end
function Base.show(io::IO, ::MIME"text/latex", x::DE_Term)
    latex_string = terms2str(x, do_latex=true, do_braket=true)
    print(io, latexstring(latex_string))
end


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

mutable struct DE_Term_Multi          # Differential Equation Term (for Multi Qubit)
    exp_op::Term                        # Expectation value, for which we derived the time derivative
    which_ind::Vector{String}           # Indexes contained in exp_op
    sum_ind::String                     # Indexes to sum over
    sum_i_terms::Vector{Term}           # Terms to sum over (in py terms[0])
    sum_i_neq_j_terms::Vector{Term}     # Terms to sum over (in py terms[1])
    non_sum_terms::Vector{Term}         # Terms that do not contain a sum (in py terms[2] - i=ind[0] + ... + i=ind[-1] + term_kappa+drive_terms)

    # Inner constructor with default values
    function DE_Term_Multi(exp_op::Term, which_ind::Vector{String}, sum_ind::String; sum_i_terms::Vector{Term}=Term[], sum_i_neq_j_terms::Vector{Term}=Term[], non_sum_terms::Vector{Term}=Term[])
        new(exp_op, which_ind, sum_ind, sum_i_terms, sum_i_neq_j_terms, non_sum_terms)
    end
end

function terms2str(multi_term::DE_Term_Multi; group_terms::Bool=true, do_sigma::Bool=false, do_latex::Bool=true, do_braket::Bool=true)::String
    # Extract data from a DE_Term_Multi and create a string from it
    exp_op::Term = multi_term.exp_op
    which_ind = multi_term.which_ind
    sum_ind::String = multi_term.sum_ind
    sum_i::Vector{Term} = multi_term.sum_i_terms
    sum_i_neq_ind::Vector{Term} = multi_term.sum_i_neq_j_terms
    rest::Vector{Term} = multi_term.non_sum_terms
    # create the string
    eq_str::String = ""
    exp_op_str::String = term2str(exp_op, do_sigma=do_sigma, do_latex=do_latex, do_braket=do_braket)
    sum_i_str::String = ""
    if length(sum_i) > 0
        sum_i_str = sumstring(sum_ind, do_latex=do_latex)

        if length(sum_ind) > 1
            if do_latex
                sum_i_neq_ind_str *= "\\left("
            else
                sum_i_neq_ind_str *= "("
            end
        end
        sum_i_str *= terms2str(sum_i, group_terms=group_terms, do_sigma=do_sigma, do_latex=do_latex, do_braket=do_braket)
        if length(sum_ind) > 1
            if do_latex
                sum_i_neq_ind_str *= "\\right)"
            else
                sum_i_neq_ind_str *= ")"
            end
        end
    end
    sum_i_neq_ind_str::String = ""
    if length(sum_i_neq_ind) > 0
        sum_i_neq_ind_str = sumstring(sum_ind, which_ind, do_latex=do_latex)
        if length(sum_i_neq_ind) > 1
            if do_latex
                sum_i_neq_ind_str *= "\\left("
            else
                sum_i_neq_ind_str *= "("
            end
        end
        sum_i_neq_ind_str *= terms2str(sum_i_neq_ind, group_terms=group_terms, do_sigma=do_sigma, do_latex=do_latex, do_braket=do_braket)
        if length(sum_i_neq_ind) > 1
            if do_latex
                sum_i_neq_ind_str *= "\\right)"
            else
                sum_i_neq_ind_str *= ")"
            end
        end
    end
    rest_str = ""
    if length(rest) > 0
        rest_str = terms2str(rest, group_terms=group_terms, do_sigma=do_sigma, do_latex=do_latex, do_braket=do_braket)
    end
    # put it all together
    if do_latex
        eq_str = "\\frac{\\mathrm{d} }{\\mathrm{d}t} " * exp_op_str * " = "
    else
        eq_str = "d/dt" * exp_op_str * " = " #* sum_i_str * " + " * sum_i_neq_ind_str
    end
    if length(sum_i_str) > 0
        eq_str *= sum_i_str
    end
    if length(sum_i_neq_ind_str) > 0
        if length(sum_i_str) > 0
            eq_str *= " + "
        end
        eq_str *= sum_i_neq_ind_str
    end
    if length(rest_str) > 0
        if !(rest_str[1] == "-")
            if length(sum_i_str) > 0 || length(sum_i_neq_ind_str) > 0
                eq_str *= " + "
            end
            eq_str *= rest_str
        else
            eq_str *= " - " * rest_str[2:end]
        end
    end
    return eq_str
end

function Base.show(io::IO, ::MIME"text/plain", x::DE_Term_Multi)
    latex_string = terms2str(x, do_latex=false)
    print(io, latex_string)
end
function Base.show(io::IO, ::MIME"text/latex", x::DE_Term_Multi)
    latex_string = terms2str(x, do_latex=true)
    print(io, latexstring(latex_string))
end


function pretty_table_from_tuple_vector(tuple_vector, names::Vector{String})
    # create a pretty table from a vector of tuples
    # names is a vector of strings
    # Example:
    # tuple_vector = [(1, 2, 3), (4, 5, 6)]
    # names = ["a", "b", "c"]
    # pretty_table_from_tuple_vector(tuple_vector, names)
    # a | b | c
    # 1 | 2 | 3
    # 4 | 5 | 6
    arrays = permutedims(hcat([[string(a) for a in A] for A in tuple_vector]...), (2, 1))
    pretty_table(arrays, header=names)
end
## Test
#x = all_eqs_indexed[1].linear_terms
#names = ["i(Op.).", "Conj.", "i(Var.)", "Coeff."]
#pretty_table_from_tuple_vector(x, names)


function pretty_table_from_indexed_terms(indexed_term_vector::Vector{Indexed_Term}, names::Vector{String})
    # create a pretty table from a vector of tuples
    # names is a vector of strings
    # Example:
    # tuple_vector = [(1, 2, 3), (4, 5, 6)]
    # names = ["a", "b", "c"]
    # pretty_table_from_tuple_vector(tuple_vector, names)
    # a | b | c
    # 1 | 2 | 3
    # 4 | 5 | 6
    # indexed_terms to tuple_vector
    function indexed_term2tuple(indexed_term::Indexed_Term)
        return (indexed_term.operator_index, indexed_term.conjugate, indexed_term.variables_index, indexed_term.coefficient)
    end
    tuple_vector = [indexed_term2tuple(indexed_term) for indexed_term in indexed_term_vector]
    if length(tuple_vector) > 0
        arrays = permutedims(hcat([[string(a) for a in A] for A in tuple_vector]...), (2, 1))
        pretty_table(arrays, header=names)
    end
end

function pretty_table_from_constant_terms(indexed_term_vector::Vector{Constant_Term}, names::Vector{String})
    # create a pretty table from a vector of tuples
    # names is a vector of strings
    # Example:
    # tuple_vector = [(1, 2, 3), (4, 5, 6)]
    # names = ["a", "b", "c"]
    # pretty_table_from_tuple_vector(tuple_vector, names)
    # a | b | c
    # 1 | 2 | 3
    # 4 | 5 | 6
    # indexed_terms to tuple_vector
    function constant_term2tuple(constant_term::Constant_Term)
        return (constant_term.variables_index, constant_term.coefficient)
    end
    tuple_vector = [constant_term2tuple(indexed_term) for indexed_term in indexed_term_vector]
    if length(tuple_vector) > 0
        arrays = permutedims(hcat([[string(a) for a in A] for A in tuple_vector]...), (2, 1))
        pretty_table(arrays, header=names)
    end
end
## Test
#x = all_eqs_indexed[1].linear_terms
#names = ["i(Op.).", "Conj.", "i(Var.)", "Coeff."]
#pretty_table_from_tuple_vector(x, names)



#### String formatting ####
function exp2str(float_value::Float64, precision::Int=2)
    exponent = Int(floor(log10(abs(float_value))))
    if exponent < 0
        if exponent < -precision
            f = Printf.Format("%0." * string(precision) * "e")
        else
            f = Printf.Format("%0." * string(precision) * "f")
        end
    else
        if exponent > precision
            f = Printf.Format("%0." * string(precision) * "e")
        else
            f = Printf.Format("%0." * string(max(precision - exponent, 1)) * "f")
        end
    end
    return Printf.format(f, float_value)
end
function fprintln(args...; digits::Int=2)
    # transform all args, if arg in args is a Float64, then transform it to a string with digits digits (default 2)
    new_args = []
    for arg in args
        if typeof(arg) == Float64
            push!(new_args, exp2str(arg, digits))
        else
            push!(new_args, arg)
        end
    end
    println(new_args...)
end


#### Time printing

function time2str(time_in_seconds::Float64)::String
    if time_in_seconds < 10^12
        hours::Int = floor(time_in_seconds / 3600)
        remaining::Float64 = time_in_seconds - hours * 3600
        minutes::Int = floor(remaining / 60)
        seconds::Int = floor(remaining - minutes * 60)
        if hours < 10
            hour_str = "0" * string(hours)
        else
            hour_str = string(hours)
        end
        if minutes < 10
            min_str = "0" * string(minutes)
        else
            min_str = string(minutes)
        end
        if seconds < 10
            sec_str = "0" * string(seconds)
        else
            sec_str = string(seconds)
        end
        time_str::String = ""
        if hours > 0
            time_str = hour_str * ":" * min_str * ":" * sec_str
        else
            time_str = min_str * ":" * sec_str
        end
        return time_str
    else
        return "∞"
    end
end

function elapsed_remaining_time_str(t_start, status, total)::Tuple{String,String}
    curr_elapsed::Float64 = time() - t_start
    ratio::Float64 = status / total
    curr_remaining::Float64 = curr_elapsed / ratio - curr_elapsed
    # separate into minutes and seconds
    elapsed_str::String = time2str(curr_elapsed)
    remaining_str::String = time2str(curr_remaining)
    return elapsed_str, remaining_str
end

function dotted_str(num::Int)
    str = string(num)
    # add dots for thousand million... numbers
    new_str = ""
    for i in 1:length(str)
        if i > 1 && (length(str) - i + 1) % 3 == 0
            new_str *= "."
        end
        new_str *= str[i]
    end
    return new_str
end

function progress(iterator, how_many=-1; text::String="", printing::Bool=true)
    # what is the length of the iterator
    #iterator = collect(it)
    len = length(iterator)
    if how_many < 0
        how_many = len
    end
    # when to print the progress bar so that len is optimally split into how_many parts
    how_many_ints::Int = 0
    t0 = time()
    dot_len = dotted_str(len)
    # Print underscore line ____ above to show width of progress bar
    chnl = Channel() do channel
        for (counter, elem) in enumerate(iterator)
            curr_floor = floor(Int, counter / len * how_many)
            if curr_floor > how_many_ints
                if printing
                    elapsed_str, remaining_str = elapsed_remaining_time_str(t0, counter, len)
                    println("$(text) ", lpad(dotted_str(counter), length(dot_len)), "/", dot_len, " ( Elapsed: ", elapsed_str, " - Remaining: ", remaining_str, " )")
                end
                how_many_ints = curr_floor
            end
            put!(channel, elem)
        end
    end
    return chnl
end
function progress(iterator, how_many=-1; text::String="", printing::Bool=true, alternative_counter::Vector{Int}=Int[], min_for_counting::Int=-1)
    # what is the length of the iterator
    # check if length(iterator) is equal to length(alternative_counter)
    len_it = length(iterator)
    len = length(alternative_counter)
    total_sum = sum(alternative_counter)
    if len == 0
        alternative_counter = ones(Int, length(iterator))
    end
    if len_it != len
        error("Length of iterator and alternative_counter must be the same.")
    end
    if how_many < 0
        how_many = len
    end
    if min_for_counting > 0
        if total_sum < min_for_counting
            printing = false
        end
    end
    # when to print the progress bar so that len is optimally split into how_many parts
    how_many_ints::Int = 0
    counter::Int = 0
    t0 = time()
    dot_len = dotted_str(total_sum)
    index::Int = 0
    # Print underscore line ____ above to show width of progress bar
    chnl = Channel() do channel
        for (curr_counter, next1) in zip(alternative_counter, iterator)
            counter += curr_counter
            curr_floor = floor(Int, counter / total_sum * how_many)
            #println(next1, " ", index, " ", curr_floor, " ", how_many_ints)
            if curr_floor > how_many_ints
                if printing
                    elapsed_str, remaining_str = elapsed_remaining_time_str(t0, counter, total_sum)
                    println("$(text) ", lpad(dotted_str(counter), length(dot_len)), "/", dot_len, " ( Elapsed: ", elapsed_str, " - Remaining: ", remaining_str, " )")
                end
                how_many_ints = curr_floor
            end
            put!(channel, next1)
        end
    end
    return chnl
end