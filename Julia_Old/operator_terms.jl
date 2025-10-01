using Base: deepcopy
using LinearAlgebra: conj

#### CONSTRUCT TERMS FROM STRING #################################################################

mutable struct Term
    spin_types::String
    spin_indices::String
    bosons::String
    exponents::Vector{Int}
    coeff::ComplexF64
    vars::Vector{String}
    conjugate::Bool # only use after computations are done (not taken into account by multiplication or commutators...)
    is_pauli::Bool

    # Inner constructor with default values
    function Term(; spin_types::String="", spin_indices::String="", bosons::String="", exponents::Vector{Int}=Int[], coeff::ComplexF64=0.0 + 0.0im, vars::Vector{String}=String[], conjugate::Bool=false, is_pauli::Bool=true)
        new(spin_types, spin_indices, bosons, copy(exponents), coeff, vars, conjugate, is_pauli)
    end
end
## Test
#term = Term(spin_types="xy", spin_indices="ij")
##term.spin_types = "xy"
##term.spin_indices = "ij"
#term.bosons = "+-"
#term.exponents = [1, 2]
#term.coeff = 1.0+0.0im
#print(term)

function str2int_array(s::String)::Vector{Int}
    # transform string into integer array with an integer for each character
    return Int.(collect(s))
end

function int_array2str(a::Vector{Int})::String
    # transform integer array into string
    return join(Char.(a))
end

function sort_spins_by_index!(term::Term)
    indices::String = term.spin_indices
    types::String = term.spin_types
    int_indices = str2int_array(indices)
    ind = sortperm(int_indices)
    int_indices = int_indices[ind]
    term.spin_types = types[ind]
    term.spin_indices = int_array2str(int_indices)
end
## test
#term = Term(spin_types="xyz", spin_indices="kji")
#display(term)
#sort_spins_by_index!(term)
#display(term)

function type_chars2index(s::String; is_pauli::Bool=true)::Vector{Int}
    # transforms xyz+- into 12345
    indexes = str2int_array(s)
    ind_vec::Vector{Int} = []
    if is_pauli
        ind_vec = [120, 121, 122, 43, 45] # x, y, z, +, -
    else
        ind_vec = [112, 109, 122, 73, 43, 45] # p, m, z, I, +, -
    end
    # replace indexes element wise 
    for (ind, i) in enumerate(ind_vec)
        indexes[indexes.==i] .= ind
    end
    if any(indexes .< 1) || any(indexes .> length(ind_vec))
        error("type_chars2index: invalid index")
    end
    return indexes
end

function index2type_chars(indexes::Vector{Int}; is_pauli::Bool=true)::String
    # transforms 12345 into xyz+-
    if is_pauli
        type_chars = ['x', 'y', 'z', '+', '-']
    else
        type_chars = ['p', 'm', 'z', 'I', '+', '-']
    end
    str = ""
    for i in indexes
        if i < 1 || i > length(type_chars)
            error("index2type_chars: invalid index")
        end
        str *= type_chars[i]
    end
    return str
end
## test
#indexes = type_chars2index("xyz+-")
#display(indexes)
#string = index2type_chars(indexes)
#display(string)

function get_spin_coeff()::Tuple{Matrix{Int},Matrix{ComplexF64}}
    # Generates a cayley table for the spin product coefficients, separated into operator index and coefficient
    # Returns tables for spin product coefficients x,y,z,I are on indexes i= 1,2,3,4 respectively
    ij_k = Int[4 3 2 1;
        3 4 1 2;
        2 1 4 3;
        1 2 3 4]
    ij_coeff = ComplexF64[1 im -im 1;
        -im 1 im 1;
        im -im 1 1;
        1 1 1 1]
    return ij_k, ij_coeff
end
## test
#ij_k, ij_coeff = get_spin_coeff()
#display(ij_k)
#display(ij_coeff)
function get_spin_coeff_pm()::Tuple{Matrix{Int},Matrix{Int},Matrix{ComplexF64},Matrix{ComplexF64}}
    # 1,2,3,4 are m,p,z,I (0 means not used! and if used anyways would lead to an indexing error)
    #           p,m,z,I 
    ij_k1 = Int[4 4 1 1; # p 
        4 4 2 2; # m
        1 2 4 3; # z
        1 2 3 4] # I
    #           p,m,z,I
    ij_k2 = Int[0 3 0 0; # p 
        3 0 0 0; # m
        0 0 0 0; # z
        0 0 0 0] # I

    ij_coeff1 = ComplexF64[0.0 0.5 -1.0 1.0;
        0.5 0.0 1.0 1.0;
        1.0 -1.0 1.0 1.0;
        1.0 1.0 1.0 1.0]
    ij_coeff2 = ComplexF64[0 0.5 0 0;
        -0.5 0 0 0;
        0 0 0 0;
        0 0 0 0]
    return ij_k1, ij_k2, ij_coeff1, ij_coeff2
end
function get_spin_coeff_pm_how_many()::Matrix{Int}
    how_many = Int[1 2 1 1; # p 
        2 1 1 1; # m
        1 1 1 1; # z
        1 1 1 1] # I
    return how_many
end

function pauli_simplify_terms(old_spin_type_ind::Vector{Int}, old_spin_indices::String, old_coeff::ComplexF64)::Tuple{Vector{Int},String,ComplexF64}
    ij_k, ij_coeff = get_spin_coeff()
    i = 1

    spin_type_ind = deepcopy(old_spin_type_ind)
    spin_indices = deepcopy(old_spin_indices)
    coeff = deepcopy(old_coeff)
    n = length(spin_type_ind)
    while i <= n - 1
        if spin_indices[i] == spin_indices[i+1]
            j = spin_type_ind[i]
            k = spin_type_ind[i+1]
            new_type = ij_k[j, k]
            new_coeff = ij_coeff[j, k]
            coeff *= new_coeff

            if new_type == 4
                deleteat!(spin_type_ind, [i, i + 1])
                spin_indices = spin_indices[1:i-1] * spin_indices[i+2:end]
                n -= 2
            else
                spin_type_ind[i] = new_type
                deleteat!(spin_type_ind, i + 1)
                spin_indices = spin_indices[1:i] * spin_indices[i+2:end]
                n -= 1
            end
        else
            i += 1
        end
    end
    return spin_type_ind, spin_indices, coeff
end

function remove_four!(spin_type_ind::Vector{Int}, spin_indices::String)
    # remove indeices that are 4 from spin_type_ind and spin_indices
    i::Int = 1
    while i <= length(spin_type_ind)
        if spin_type_ind[i] == 4
            deleteat!(spin_type_ind, i)
            spin_indices = spin_indices[1:i-1] * spin_indices[i+1:end]
        else
            i += 1
        end
    end
end
function remove_four!(spin_type_inds::Vector{Vector{Int}}, spin_indices::Vector{String})
    for x in 1:length(spin_type_inds)
        i::Int = 1
        while i <= length(spin_type_inds[x])
            if spin_type_inds[x][i] == 4
                deleteat!(spin_type_inds[x], i)
                spin_indices[x] = spin_indices[x][1:i-1] * spin_indices[x][i+1:end]
            else
                i += 1
            end
        end
    end
end
function pm_simplify_terms(spin_type_ind::Vector{Int}, spin_indices::String, coeff::ComplexF64)
    ij_k1, ij_k2, ij_coeff1, ij_coeff2 = get_spin_coeff_pm()
    how_many = get_spin_coeff_pm_how_many()
    i::Int = 1
    n::Int = length(spin_type_ind)
    how_long::Vector{Int} = []
    while i <= n
        if i <= n - 1 && spin_indices[i] == spin_indices[i+1]
            j = spin_type_ind[i]
            k = spin_type_ind[i+1]
            push!(how_long, how_many[j, k])
            i += 2
        else
            push!(how_long, 1)
            i += 1
        end
    end
    how_many_total::Int = 2 .^ (sum(how_long .- 1))
    #println(spin_type_ind, " ", spin_indices, " ", how_long, " ", how_many_total)
    # Now we create empty arrays for all combinations 
    len::Int = length(how_long)
    new_spin_type_ind::Vector{Vector{Int}} = [zeros(Int, len) for i in 1:how_many_total]
    new_spin_indices::String = ""
    new_coeffs::Vector{ComplexF64} = [coeff for i in 1:how_many_total]
    how_many_prev::Int = 1 # how often the schema is repeated in the list before element 1 gets replaced by element 2
    i = 1
    curr_i = 1
    while i <= n
        if i <= n - 1 && spin_indices[i] == spin_indices[i+1]
            j = spin_type_ind[i]
            k = spin_type_ind[i+1]
            type1 = ij_k1[j, k]
            coeff1 = ij_coeff1[j, k]
            curr_how_many = how_many[j, k]
            new_spin_indices *= spin_indices[i]
            if curr_how_many > 1
                type2 = ij_k2[j, k]
                coeff2 = ij_coeff2[j, k]
                ji = 1
                # do how_many_prev repetition of adding type_1, then equally many type_2 and repeat this until how_many_total is reached to fill up curr_i elements of the new_spin_type_ind
                while ji <= how_many_total
                    for curr_j in ji:ji+how_many_prev-1
                        new_spin_type_ind[curr_j][curr_i] = type1
                        new_coeffs[curr_j] *= coeff1
                    end
                    ji += how_many_prev
                    for curr_j in ji:ji+how_many_prev-1
                        new_spin_type_ind[curr_j][curr_i] = type2
                        new_coeffs[curr_j] *= coeff2
                    end
                    ji += how_many_prev
                end
                how_many_prev *= curr_how_many
            else
                for j in 1:how_many_total
                    new_spin_type_ind[j][curr_i] = type1
                    new_coeffs[j] *= coeff1
                end
            end
            i += 2
        else
            j = spin_type_ind[i]
            new_spin_indices *= spin_indices[i]
            for k in 1:how_many_total
                new_spin_type_ind[k][curr_i] = j
            end
            i += 1
        end
        curr_i += 1
    end

    # clean up, remove spin_type_ind == 4 -> I is not needed 
    repeated_new_spin_indices::Vector{String} = [new_spin_indices for i in 1:how_many_total]
    remove_four!(new_spin_type_ind, repeated_new_spin_indices)
    # check if any new_spin_indices are still double?
    # if there are still double indices, we use pm_simplify_terms again on every element 
    # of the new_spin_type_ind and new_spin_indices
    any_doubles::Bool = false
    if length(unique(new_spin_indices)) != length(new_spin_indices)
        any_doubles = true
    end
    if any_doubles
        super_new_spin_indices::Vector{String} = []
        super_new_spin_type_ind::Vector{Vector{Int}} = []
        super_new_coeffs::Vector{ComplexF64} = []
        for i in 1:how_many_total
            my_new_spin_type_ind, curr_new_spin_indices, my_new_coeffs = pm_simplify_terms(new_spin_type_ind[i], repeated_new_spin_indices[i], new_coeffs[i])
            if abs(my_new_coeffs) > 10^-14
                remove_four!(my_new_spin_type_ind, curr_new_spin_indices)
                append!(super_new_spin_type_ind, my_new_spin_type_ind)
                append!(super_new_spin_indices, curr_new_spin_indices)
                append!(super_new_coeffs, my_new_coeffs)
            end
        end
        return super_new_spin_type_ind, super_new_spin_indices, super_new_coeffs
    end
    return new_spin_type_ind, repeated_new_spin_indices, new_coeffs
end

function simplify_spin_terms(term::Term)
    # Simplifies spin terms by combining operators acting on the same subspace (same index) -> sorts by index first
    sort_spins_by_index!(term)
    is_pauli = term.is_pauli
    spin_types = term.spin_types
    spin_indices = term.spin_indices
    coeff = term.coeff
    #display(term)
    n = length(spin_types)
    if length(spin_types) != length(spin_indices)
        throw(ArgumentError("Invalid spin term, spin_types and spin_indices must be the same length"))
    end
    if n >= 2
        spin_type_ind = type_chars2index(spin_types, is_pauli=is_pauli)
        if is_pauli
            new_spin_type_ind, new_spin_indices, new_coeff = pauli_simplify_terms(spin_type_ind, spin_indices, coeff)
            new_spin_types = join(index2type_chars(new_spin_type_ind, is_pauli=is_pauli))
            return Term(spin_types=new_spin_types, spin_indices=new_spin_indices, bosons=term.bosons, exponents=term.exponents, coeff=new_coeff, vars=term.vars, is_pauli=is_pauli)
        else
            new_spin_type_ind, new_spin_indices, new_coeff = pm_simplify_terms(spin_type_ind, spin_indices, coeff)
            new_terms::Vector{Term} = []
            for i in 1:length(new_spin_type_ind)
                new_spin_types = join(index2type_chars(new_spin_type_ind[i], is_pauli=is_pauli))
                push!(new_terms, Term(spin_types=new_spin_types, spin_indices=new_spin_indices[i], bosons=term.bosons, exponents=term.exponents, coeff=new_coeff[i], vars=term.vars, is_pauli=is_pauli))
            end
            if length(new_terms) == 1
                return new_terms[1]
            end
            return new_terms
        end
    end
    return term
end
## test
#spin_types = "xyzyyzzx"
#spin_indices = "iikiijjj"
#term = Term(spin_types=spin_types, spin_indices=spin_indices)
#term.coeff = 1.0+0.0im
#display(term)
#term = simplify_spin_terms(term)
#display(term)

function clean_up_bosons!(term::Term)
    # Check for double -- and ++, reduce them to - and +, respectively and add up their exponents
    bosons::String = term.bosons
    exponents::Vector{Int} = term.exponents
    i = 1
    if length(bosons) != length(exponents)
        throw(ArgumentError("Invalid boson types and exponents, not same length"))
    end

    while i < length(bosons)
        if bosons[i] == bosons[i+1]
            # Same type, add exponents and remove one
            exponents[i] += exponents[i+1]
            bosons = bosons[1:i] * bosons[i+2:end]
            deleteat!(exponents, i + 1)
        else
            i += 1
        end
    end
    term.bosons = bosons
    term.exponents = exponents
end
function clean_up_bosons!(terms::Vector{Term})
    for term in terms
        clean_up_bosons!(term)
    end
end
## test
#term.bosons = "++--"
#term.exponents = [1, 1, 1, 1]
#display(term)
#clean_up_bosons!(term)
#display(term)

function remove_borders!(term::Term)
    bosons = term.bosons
    exponents = term.exponents
    bosons, exponents = clean_up_bosons!(bosons, exponents)
    if !isempty(bosons)
        if exponents[1] == 0
            bosons = bosons[2:end]
            exponents = exponents[2:end]
        end
        if exponents[end] == 0
            bosons = bosons[1:end-1]
            exponents = exponents[1:end-1]
        end
    end
    term.bosons = bosons
    term.exponents = exponents
end

function add_empty_borders!(term::Term)
    bosons = term.bosons
    exponents = term.exponents
    if !isempty(bosons)
        if bosons[1] != '+'
            bosons = "+" * bosons
            exponents = vcat(0, exponents)
        end
        if bosons[end] != '-'
            bosons = bosons * "-"
            exponents = vcat(exponents, 0)
        end
    else
        bosons = "+-"
        exponents = [0, 0]
    end
    term.bosons = bosons
    term.exponents = exponents
end

function remove_empty_elements!(term::Term, borders::Bool=false)
    # Remove empty elements from bosons and exponents and add empty borders if needed
    clean_up_bosons!(term)
    bosons = term.bosons
    exponents = term.exponents
    border_int::Int = borders ? 0 : 1
    i = border_int + 1

    while i <= length(exponents) - border_int
        if exponents[i] == 0
            bosons = bosons[1:i-1] * bosons[i+1:end]
            exponents = vcat(exponents[1:i-1], exponents[i+1:end])
        else
            i += 1
        end
    end
    term.bosons = bosons
    term.exponents = exponents
    clean_up_bosons!(term)
end
## test
#term.bosons = "++++-"
#term.exponents = [0, 1, 0, 1, 0]
#display(term)
#remove_empty_elements!(term)
#display(term)
#display("-"^50)
#term.bosons = "+-"
#term.exponents = [0, 0]
#display(term)
#remove_empty_elements!(term, true)
#display(term)
function term_str_to_term_vec(str_term::String)::Vector{String}
    op_types = ['x', 'y', 'z', '+', '-']
    str_term = replace(str_term, "_" => "")
    for s in op_types
        str_term = replace(str_term, s => "*" * s)
    end
    string_list::Vector{String} = split(str_term, "*")
    # remove empty strings in op_list
    string_list = filter(p -> p != "", string_list)
    return string_list
end
## Test
#println(term_str_to_term_vec("++x_iy_j"))

function make_term(str_term::String, coeff=1.0 + 0im, vars::Union{String,Vector{String}}=String[], conjugate::Bool=false; optim::Bool=true, is_pauli::Bool=true)::Union{Term,Vector{Term}}
    # * is needed to separate the operators
    if !isa(coeff, ComplexF64)
        coeff = ComplexF64(coeff)
    end
    if isa(vars, String)
        vars = String[vars]
    end
    # Create Term struct from string
    term = Term()
    indices::String = ""
    types::String = ""

    if is_pauli
        pauli_strs = ['x', 'y', 'z']
    else
        pauli_strs = ['m', 'p', 'z']
    end
    cavity_str::String = ""
    new_str_term::String = ""
    for s in str_term
        if s in ['+', '-']
            cavity_str *= s
        else
            new_str_term *= s
        end
    end
    term.bosons = cavity_str
    str_term = new_str_term

    str_term = replace(str_term, "_" => "")
    string_list = split(str_term, "*")
    # remove empty strings in op_list
    string_list = filter(p -> p != "", string_list)
    for (i, p) in enumerate(string_list)
        if p[1] in pauli_strs
            types *= p[1]
            if length(p) == 1
                indices *= " "
            else
                indices *= p[2:end]
            end
        else
            error("Operator not recognized")
        end
    end

    term.spin_types = types
    term.spin_indices = indices
    term.coeff = coeff
    term.vars = vars
    term.exponents = ones(Int, length(term.bosons))
    term.conjugate = conjugate
    term.is_pauli = is_pauli
    if optim
        ## Add post processing
        # 1st sort spins by sort_spins_by_index
        # 2nd simplify spin terms
        term = simplify_spin_terms(term)
        if isa(term, Vector{Term}) && length(term) == 1
            term = term[1]
        end
        # 3rd clean up bosons
        clean_up_bosons!(term)
    end
    return term
end
function make_terms(str_terms::Vector{String}, coeff=1.0 + 0im, vars::Union{String,Vector{String}}=String[], conjugate::Bool=false; optim::Bool=true, is_pauli::Bool=true)::Vector{Term}
    terms::Vector{Term} = []
    for str_term in str_terms
        term = make_term(str_term, coeff, vars, conjugate, optim=optim, is_pauli=is_pauli)
        push!(terms, term)
    end
    return terms
end

function string2term(str_term::String, coeff=1.0 + 0im, vars::Union{String,Vector{String}}=String[])::Term
    return make_term(str_term, coeff, vars)
end
## test the function
#term = make_term("+-x_j*yi")
#display(term)
#term = make_term("+-xjyi", 1.0+0.0im, ["a", "b"])
#display(term)
#term = make_term("xy", 1.0+0.0im, ["a", "b"])
#display(term)

function terms_combinable(A::Term, B::Term)
    if !(A.spin_types == B.spin_types)
        return false
    elseif !(A.spin_indices == B.spin_indices)
        return false
    elseif !(A.bosons == B.bosons)
        return false
    elseif !(A.exponents == B.exponents)
        return false
    elseif !(A.vars == B.vars)
        return false
    else
        return true
    end
end

function combine_terms(terms::Vector{Term}, eps::Float64=1e-14)::Vector{Term}
    term_list = deepcopy(terms)
    i = 1
    while i <= length(term_list)
        j = i + 1
        while j <= length(term_list)
            if terms_combinable(term_list[i], term_list[j])
                term_list[i].coeff += term_list[j].coeff
                deleteat!(term_list, j)
            else
                j += 1
            end
        end
        i += 1
    end
    # remove vanishingly small terms
    term_list = filter(term -> abs(term.coeff) > eps, term_list)
    return term_list
end
## test
#terms = [make_term("+y_j*x_i"), make_term("+x_i*y_j"), make_term("+x_i*y_j")]
#term_list = combine_terms(terms)
#display(term_list[1])

function normal_order(term::Term)::Vector{Term}
    # Create normal ordering of a term, returns a Vector of Terms)
    function flip_first_plus_element!(term::Term)::Term
        remove_empty_elements!(term, false)
        term_copy = deepcopy(term)
        bosons = term.bosons
        exponents = term.exponents
        if length(term.bosons) > 2 && term.exponents[3] > 0
            if term.bosons[2] != '-'
                error("Invalid boson types, bosons[2] is not -. (", term.bosons, ", ", term.exponents, ")")
            end
            if term.bosons[3] != '+'
                error("Invalid boson types, bosons[3] is not +. (", term.bosons, ", ", term.exponents, ")")
            end
            multiplier = term.exponents[2]
            term_copy.exponents[2] -= 1
            term_copy.exponents[3] -= 1
            term.exponents[1] += 1
            term.exponents[3] -= 1
            term_copy.coeff *= multiplier

            remove_empty_elements!(term, false)
            remove_empty_elements!(term_copy, false)
        else
            term_copy.coeff = 0.0 + 0.0im
        end
        return term_copy
    end
    term2 = deepcopy(term)
    remove_empty_elements!(term2, true)
    add_empty_borders!(term2)

    term_copy = flip_first_plus_element!(term2)
    term_array::Vector{Term} = Term[]
    if term_copy.coeff == 0
        term_array = [term]
    else
        term_array = normal_order(term2)
        append!(term_array, normal_order(term_copy))
    end
    # remove duplicates and vanishingly small terms
    for i in 1:length(term_array)
        remove_empty_elements!(term_array[i], true)
    end
    term_array = combine_terms(term_array)
    return term_array
end
function normal_order(terms::Vector{Term})::Vector{Term}
    term_array::Vector{Term} = Term[]
    for term in terms
        term_array = vcat(term_array, normal_order(term))
    end
    return term_array
end
## Test
#term.bosons = "+-+"
#term.exponents = [0,2,2]
#term.coeff = 1.0+0.0im
#display(term)
#term_array = normal_order(term)
#for term in term_array
#    bos = term.bosons
#    expo = term.exponents
#    coeff = term.coeff
#    println(bos, " ", expo, " ", coeff)
#end

function check_is_pauli(term::Term)::Tuple{Bool,Bool} # is_pauli, matters (length > 0)
    return term.is_pauli, true
end
function check_is_pauli(term::Vector{Term})::Tuple{Bool,Bool} # is_pauli, matters (length > 0)
    if length(term) == 0
        return true, false
    end
    all_same_pauli::Bool = true
    is_pauli::Bool = term[1].is_pauli
    for t in term
        if t.is_pauli != is_pauli
            all_same_pauli = false
            break
        end
    end
    if !all_same_pauli
        # raise error 
        error("All terms in the vector should be of the same type")
    end
    return is_pauli, true
end
function compare_is_pauli(termA::Union{Term,Vector{Term}}, termB::Union{Term,Vector{Term}})::Tuple{Bool,Bool}
    is_pauli_A, matters_A = check_is_pauli(termA)
    is_pauli_B, matters_B = check_is_pauli(termB)

    # If neither term matters, return default values
    if !matters_A && !matters_B
        return true, false
    end
    # If only one of the terms matters, return its values
    if matters_A && !matters_B
        return is_pauli_A, true
    elseif !matters_A && matters_B
        return is_pauli_B, true
    end

    # If both terms matter, compare their properties
    if !(is_pauli_A == is_pauli_B)
        error("Both terms should be of the same type")
    end
    return is_pauli_A, true
end

function multiply_terms(A::Union{Term,Vector{Term}}, B::Union{Term,Vector{Term}}; normalize::Bool=true)::Union{Term,Vector{Term}}
    is_pauli, matters = compare_is_pauli(A, B)
    if isa(A, Vector{Term})
        res = []
        for term in A
            curr_res = multiply_terms(term, B, normalize=normalize)
            for c in curr_res
                if !isa(c, Vector{Term})
                    push!(res, c)
                else
                    append!(res, c)
                end
            end
        end
        return combine_terms(res)
    elseif isa(B, Vector{Term})
        res::Vector{Term} = []
        for term in B
            curr_res = multiply_terms(A, term, normalize=normalize)
            for c in curr_res
                if !isa(c, Vector{Term})
                    push!(res, c)
                else
                    append!(res, c)
                end
            end
        end
        return combine_terms(res)
    else
        # Add implementation for `get_spin_coeff`
        # Make sure to define and import all used methods and fields here
        spin_types::String = A.spin_types * B.spin_types
        spin_indices::String = A.spin_indices * B.spin_indices
        bosons::String = A.bosons * B.bosons
        exponents::Vector{Int} = vcat(A.exponents, B.exponents)
        coeff::ComplexF64 = A.coeff * B.coeff
        vars::Vector{String} = filter(x -> !isempty(x), vcat(A.vars, B.vars))
        term = Term(spin_types=spin_types, spin_indices=spin_indices, bosons=bosons, exponents=exponents, coeff=coeff, vars=vars, is_pauli=is_pauli)
        sort_spins_by_index!(term)
        term = simplify_spin_terms(term)
        clean_up_bosons!(term)
        if normalize
            return normal_order(term)
        end
        return term
    end
end
## Test
## Generate two terms that can be multiplied
#A = make_term("+--x_i*y_j")
#B = make_term("+y_i*x_j")
#display(A)
#display(B)
#C = multiply_terms(A, B)
#display(C)

function add_terms(termA::Union{Term,Vector{Term}}, termB::Union{Term,Vector{Term}}, args::Union{Term,Vector{Term}}...)::Vector{Term}
    # Adds two terms together
    is_pauli, matters = compare_is_pauli(termA, termB) # parameters don't matter, just the error check

    term_list::Vector{Term} = Term[]
    if isa(termA, Term)
        termA = Term[termA]
    end
    if isa(termB, Term)
        termB = Term[termB]
    end
    append!(term_list, termA)
    append!(term_list, termB)
    if length(args) > 0
        for term in args
            if isa(term, Term)
                term = Term[term]
            end
            append!(term_list, term)
        end
    end
    term_list = combine_terms(term_list)
    return term_list
end
## Test
#termA = Term[make_term("+x_i*y_j")]
#termB = make_term("+-x_i*y_j")
#term_list = add_terms(termA, termB)
#display(term_list)
function flip_m_p(curr_spin_terms)
    temp_char = '#'  # Make sure this character does not appear in curr_spin_terms
    # Replace 'm' with temporary character
    step1 = replace(curr_spin_terms, 'm' => temp_char)
    # Replace 'p' with 'm'
    step2 = replace(step1, 'p' => 'm')
    # Replace temporary character with 'p'
    result = replace(step2, temp_char => 'p')

    return result
end
function dagger_term(term::Union{Term,Vector{Term}}; normalize::Bool=true, eps::Float64=1e-14)::Union{Term,Vector{Term}}
    # Does not conjugate vars!
    if isa(term, Vector)
        res::Vector{Term} = []
        for t in term
            curr_res = dagger_term(t, normalize=normalize, eps=eps)
            if !isa(curr_res, Vector{Term})
                push!(res, curr_res)
            else
                append!(res, curr_res)
            end
        end
        return res
    else
        dag_term = deepcopy(term)
        # Change the signs of bosons and reverse them 
        # replace p and m terms and for every p and m -> complex conjugate the coefficient
        curr_spin_terms = term.spin_types
        curr_spin_terms = flip_m_p(curr_spin_terms)
        dag_term.spin_types = curr_spin_terms
        dag_term.bosons = reverse(term.bosons)
        # replace + with - and - with +
        dag_term.bosons = map(x -> x == '+' ? '-' : '+', dag_term.bosons)
        dag_term.exponents = reverse(dag_term.exponents) # reverse exponents order
        dag_term.coeff = conj(dag_term.coeff) # conjugate coeff
        if normalize
            return normal_order(dag_term)
        else
            return dagger_term
        end
    end
end

function scale_term(term::Union{Term,Vector{Term}}, scale)::Union{Term,Vector{Term}}
    if isa(term, Vector{Term})
        return [scale_term(t, scale) for t in term]
    else
        new_term = deepcopy(term)
        new_term.coeff *= scale
        return new_term
    end
end

function scale_term!(term::Union{Term,Vector{Term}}, scale)
    if isa(term, Vector{Term})
        for i in 1:length(term)
            term[i].coeff *= scale
        end
    else
        term.coeff *= scale
    end
end

function commutator_terms(A::Union{Term,Vector{Term}}, B::Union{Term,Vector{Term}}, scale=1)::Union{Term,Vector{Term}}
    terms_1 = multiply_terms(A, B)
    terms_2 = scale_term(multiply_terms(B, A), -1)
    terms::Vector{Term} = vcat(terms_1, terms_2)
    if terms[1].is_pauli
        terms = combine_terms(terms)
    end
    if scale != 1
        return scale_term(terms, scale)
    else
        return terms
    end
end

function sigma_plus(index=""; is_pauli::Bool=true)::Vector{Term}
    if is_pauli
        return [make_term("x" * index, 0.5, is_pauli=true), make_term("y" * index, 0.5im, is_pauli=true)]
    else
        return [make_term("p" * index, is_pauli=false)]
    end
end

function sigma_minus(index=""; is_pauli::Bool=true)::Vector{Term}
    if is_pauli
        return [make_term("x" * index, 0.5, is_pauli=true), make_term("y" * index, -0.5im, is_pauli=true)]
    else
        return [make_term("m" * index, is_pauli=false)]
    end
end

function sigma_x(index::String=""; is_pauli::Bool=true)::Vector{Term}
    if is_pauli
        return [make_term("x" * index, is_pauli=true)]
    else
        return [make_term("p" * index, is_pauli=false), make_term("m" * index, is_pauli=false)]
    end
end
function sigma_y(index::String=""; is_pauli::Bool=true)::Vector{Term}
    if is_pauli
        return [make_term("y" * index, is_pauli=true)]
    else
        return [make_term("p" * index, -1.0im, is_pauli=false), make_term("m" * index, 1.0im, is_pauli=false)]
    end
end
function sigma_z(index::String=""; is_pauli::Bool=true)::Vector{Term}
    return [make_term("z" * index, is_pauli=is_pauli)]
end

function add_vars(term::Union{Term,Vector{Term}}, vars::Union{String,Vector{String}})::Union{Term,Vector{Term}}
    if isa(vars, String)
        vars = String[vars]
    elseif !isa(vars, Vector{String})
        error("vars should be a string or a vector of strings")
    end

    if isa(term, Vector{Term})
        return [add_vars(t, vars) for t in term]
    else
        new_term = deepcopy(term)
        append!(new_term.vars, vars)
        return new_term
    end
end
## test
#A = sigma_plus("i") 
#B = dagger_term(sigma_minus("i"))
#A = add_vars(A, "t")
#B = add_vars(B, "t")
#C = Commutator_terms(A, B)
#C = scale_term(C, 0.5)
#display(C)

function add_vars!(term::Union{Term,Vector{Term}}, vars::Union{String,Vector{String}})
    if isa(vars, String)
        vars = String[vars]
    elseif !isa(vars, Vector{String})
        error("vars should be a string or a vector of strings")
    end

    if isa(term, Vector{Term})
        for i in 1:length(term)
            append!(term[i].vars, vars)
        end
    else
        append!(term.vars, vars)
    end
end


# Get the number of operators in a term
function how_many_operators_in_term(term)::Tuple{Int,Int,Int}
    # Returns for a term:
    # total_op: total number of operators
    # ferm_op: number of fermionic operators
    # bos_op: number of bosonic operators
    ferm_op = length(term.spin_types)
    bos_op = 0
    exponents = term.exponents
    for i in 1:length(exponents)
        bos_op += exponents[i]
    end
    return ferm_op + bos_op, ferm_op, bos_op
end
## Test
#term = make_term("++-xi*yj")
#how_many_operators_in_term(term)


#exp_op = all_eqs[1][1].exp_op
#println(term2str_identifier(exp_op)) 
function term2str_identifier(term::Term, do_indexes::Bool=false)::Tuple{String,Int,Int}
    # Create string identifier for term
    # Example: " "+-xiyj"
    # if an operator has more creation (a^\dagger) than annihilation (a) operators, it is conjugated, and the conjugate flag is set to true, 
    spin_types::String = term.spin_types
    spin_indices::String = term.spin_indices
    bosons::String = term.bosons
    exponents::Vector{Int} = term.exponents
    str::String = ""
    for i in 1:length(exponents)
        str *= bosons[i]^exponents[i]
    end
    # add spin operators
    if do_indexes
        for i in 1:length(spin_indices)
            if spin_indices[i] == ' '
                str *= spin_types[i]
            else
                str *= spin_types[i] * spin_indices[i]
            end
            if i < length(spin_indices)
                str *= "*"
            end
        end
    else
        for i in 1:length(spin_types)
            str *= spin_types[i]
        end
    end
    spin_order::Int = length(spin_types)
    order::Int = sum(exponents) + spin_order
    return str, order, spin_order
end
## Test
#exp_op = make_term("++z")
#println(term2str_identifier(exp_op)) 

function string2order(term_str::String)::Int
    # transform a term string to it's order
    # Example: " "+-xi*yj" -> 4 (2 spins, 2 bosons)
    # count +, -, x, y, z, m, p
    # separate into cavity and spin operators
    # separate spin operators by * and count them
    order::Int = 0
    new_str_term::String = ""
    for s in term_str
        if s in ['+', '-']
            order += 1
        else
            new_str_term *= s
        end
    end
    string_list = split(new_str_term, "*")
    # remove empty strings in op_list
    string_list = filter(p -> p != "", string_list)
    for (i, p) in enumerate(string_list)
        if p[1] in ['x', 'y', 'z', 'm', 'p']
            order += 1
        else
            error("Operator not recognized (", p, ") in term ", term_str, " at position ", i, " of ", string_list, ".")
        end
    end
    return order
end
## Test 
#println(string2order("+-xi*yj"))

function string2order_spinorder(term_str::String)::Tuple{Int,Int,Int,Int}
    # transform a term string to it's order
    # Example: " "+-xiyj" -> 4 (2 spins, 2 bosons)
    # count +, -, x, y, z
    creation_order::Int = 0
    annihilation_order::Int = 0
    spin_order::Int = 0
    new_str_term::String = ""
    for s in term_str
        if s == '+'
            creation_order += 1
        elseif s == '-'
            annihilation_order += 1
        else
            new_str_term *= s
        end
    end
    string_list = split(new_str_term, "*")
    # remove empty strings in op_list
    string_list = filter(p -> p != "", string_list)
    for (i, p) in enumerate(string_list)
        if p[1] in ['x', 'y', 'z', 'm', 'p']
            spin_order += 1
        else
            error("Operator not recognized (", p, ") in term ", term_str, " at position ", i, " of ", string_list, ".")
        end
    end
    order::Int = creation_order + annihilation_order + spin_order
    return order, spin_order, creation_order, annihilation_order
end
## Test
#println(string2order_spinorder("+-xiyj"))

function count_occurances(str::String, char::Char)::Int
    count::Int = 0
    for i in 1:length(str)
        if str[i] == char
            count += 1
        end
    end
    return count
end
function count_occurances(str::String, chars::Vector{Char})::Vector{Int}
    counts::Vector{Int} = zeros(Int, length(chars))
    for i in 1:length(chars)
        counts[i] = count_occurances(str, chars[i])
    end
    return counts
end
function count_occurances(str::String, chars::String)::Vector{Int}
    counts::Vector{Int} = zeros(Int, length(chars))
    for i in 1:length(chars)
        counts[i] = count_occurances(str, chars[i])
    end
    return counts
end
function term2reducedstr(term::Term)::Tuple{String,Vector{Int},Vector{Int}}
    if term.is_pauli
        operators = ['x', 'y', 'z']
    else
        operators = ['p', 'm', 'z']
    end
    occurances::Vector{Int} = count_occurances(term.spin_types, operators)
    operator_string::String = ""
    cavity_occurances::Vector{Int} = zeros(Int, 2)
    for (b, e) in zip(term.bosons, term.exponents)
        operator_string *= b^e
        if b == '+'
            cavity_occurances[1] += e
        elseif b == '-'
            cavity_occurances[2] += e
        end
    end
    for (o, n) in zip(operators, occurances)
        operator_string *= o^n
    end
    return operator_string, occurances, cavity_occurances
end
## Test 
#term2reducedstr(make_term("+-mi*pj", is_pauli=false))
function str(term::Term)::String
    return term2reducedstr(term)[1]
end



#### Indexification of Operators #################################################################
mutable struct Constant_Term
    variables_index::Int
    coefficient::ComplexF64
    function Constant_Term(variables_index::Int, coefficient::ComplexF64)
        new(variables_index, coefficient)
    end
end

mutable struct Abstract_Constant_Term
    variables_index::Int
    coefficient::ComplexF64
    sample_vars::Vector{Tuple{String,Int}}
    summed::Int # 0 -> not summed, 1 -> summed, 2 -> neq_summed
    function Abstract_Constant_Term(variables_index::Int, coefficient::ComplexF64, sample_vars::Vector{Tuple{String,Int}}, summed::Int)
        new(variables_index, coefficient, sample_vars, summed)
    end
end

mutable struct Indexed_Term # Term for within DE_Term_indexed and DE_Term_Multi_indexed, all strings replaced by indexes to apropriate vectors
    # Tuple{Int,Bool,Int,ComplexF64} - Tuple consists of (operator index, conjugate, variables index, coefficient)
    operator_index::Int
    conjugate::Bool
    variables_index::Int
    coefficient::ComplexF64
    function Indexed_Term(operator_index::Int, conjugate::Bool, variables_index::Int, coefficient::ComplexF64)
        new(operator_index, conjugate, variables_index, coefficient)
    end
end

mutable struct Abstract_Term # Abstract reduced representation of a term, 
    # Abstract terms are separated into x y and z operators to allow the easy calculation of indexing into a vector of operator expectation values
    # keeps references to variables in Term form 
    term_str::String
    spin_orders::Vector{Int}
    spin_indices::Vector{Vector{Int}} # hijk -> 0123 but separated into x y and z groups 
    coeff::ComplexF64
    vars::Int
    sample_vars::Vector{Tuple{String,Int}} # (var, index) 
    conjugate::Bool
    summed::Int # 0 -> not summed, 1 -> summed, 2 -> neq_summed
    function Abstract_Term(term_str::String, spin_orders::Vector{Int}, spin_indices::Vector{Vector{Int}}, coeff::ComplexF64, vars::Int, sample_vars::Vector{Tuple{String,Int}}, conjugate::Bool=false, summed::Int=0)
        new(term_str, spin_orders, spin_indices, coeff, vars, sample_vars, conjugate, summed)
    end
end

# Function that transforms a Term into an Abstract_Term separating the x, y and z operators into different groups
function term_to_abstract_term(term::Term, vars_dict::Dict{String,Int}, summed::Int, min_ind::Char='h')::Abstract_Term # Modifies vars_dict!!! 
    # term is a Term from operator_terms.jl
    # returns a vector of how many of each operator is in the term counting the number of +, -, x, y and z operators
    # e.g. for the term +x*y*z*x*y*z the vector [0, 0, 2, 2, 2] is returned
    spin_orders::Vector{Int} = zeros(Int, 3)
    indices = term.spin_indices  # transform into numbers 
    inds = [Int(i) for i in indices]

    min_ind_ind = Int(min_ind)
    inds = inds .- min_ind_ind # min ind is zero
    spin_indices::Vector{Vector{Int}} = Vector{Vector{Int}}(undef, 3)
    for i in 1:3
        spin_indices[i] = []
    end
    if term.is_pauli
        operators = ['x', 'y', 'z']
    else
        operators = ['p', 'm', 'z']
    end
    for (s, i) in zip(term.spin_types, inds)
        if s == operators[1]
            spin_orders[1] += 1
            push!(spin_indices[1], i)
        elseif s == operators[2]
            spin_orders[2] += 1
            push!(spin_indices[2], i)
        elseif s == operators[3]
            spin_orders[3] += 1
            push!(spin_indices[3], i)
        else
            if term.is_pauli
                error("spin type must be x, y or z, got $s")
            else
                error("spin type must be p, m or z, got $s")
            end
        end
    end
    vars::Vector{String} = sort(term.vars)
    # separate into sample_vars (with _i or _h...) and other_vars
    other_vars::Vector{String} = []
    sample_vars::Vector{Tuple{String,Int}} = []
    for v in vars
        if '_' in v # sample variable
            prefix::String, suffix::String = split(v, "_")
            if length(suffix) > 1
                error("term_to_abstract_term: sample variable has more than one index")
            end
            suffix_ind = Int(suffix[1]) - min_ind_ind
            push!(sample_vars, (prefix, suffix_ind))
        else
            push!(other_vars, v)
        end
    end
    var_str::String = ""
    if length(other_vars) > 0
        var_str = other_vars[1]
    end
    for i in 2:length(other_vars)
        var_str *= "?" * other_vars[i]   # separate at ? (special symbol)
    end
    # if var_str is not in vars_dict, add it
    if !haskey(vars_dict, var_str)
        term_var_index = length(vars_dict) + 1
        vars_dict[var_str] = term_var_index
    else
        term_var_index = vars_dict[var_str]
    end
    # Is var in vars_dict? if not add it
    # transform operator into a string 
    term_str, _, _ = term2reducedstr(term)
    return Abstract_Term(term_str, spin_orders, spin_indices, term.coeff, term_var_index, sample_vars, term.conjugate, summed)
end

function term_to_abstract_constant_term(term::Term, vars_dict::Dict{String,Int}, summed::Int, min_ind::Char='h')::Abstract_Constant_Term
    min_ind_ind = Int(min_ind)
    vars::Vector{String} = sort(term.vars)
    # separate into sample_vars (with _i or _h...) and other_vars
    other_vars::Vector{String} = []
    sample_vars::Vector{Tuple{String,Int}} = []
    for v in vars
        if '_' in v # sample variable
            prefix::String, suffix::String = split(v, "_")
            if length(suffix) > 1
                error("term_to_abstract_term: sample variable has more than one index")
            end
            suffix_ind = Int(suffix[1]) - min_ind_ind
            push!(sample_vars, (prefix, suffix_ind))
        else
            push!(other_vars, v)
        end
    end
    var_str::String = ""
    if length(other_vars) > 0
        var_str = other_vars[1]
    end
    for i in 2:length(other_vars)
        var_str *= "?" * other_vars[i]   # separate at ? (special symbol)
    end
    # if var_str is not in vars_dict, add it
    if !haskey(vars_dict, var_str)
        term_var_index = length(vars_dict) + 1
        vars_dict[var_str] = term_var_index
    else
        term_var_index = vars_dict[var_str]
    end
    return Abstract_Constant_Term(term_var_index, term.coeff, sample_vars, summed)
end

#### Check Operator Types for op_str #################################################################
function check_consistent_operator_type(op_str::String, is_pauli::Bool)
    if !is_pauli 
        # make sure there is no x or y in op_str 
        if 'x' in op_str || 'y' in op_str
            error("Operator string $(op_str) contains x or y, which is not allowed for non-pauli systems")
        end
    else
        # make sure there is no p or m in op_str
        if 'p' in op_str || 'm' in op_str
            error("Operator string $(op_str) contains p or m, which is not allowed for pauli systems")
        end
    end
end
function check_operator_ordering(op_str::String)
    # operator ordering order should be +, -, x, p, y, m, z
    order::Vector{Char} = ['+', '-', 'x', 'p', 'y', 'm', 'z']
    last_idx::Int = 1
    for i in 1:length(op_str)
        if !(op_str[i] in order)
            error("Operator string $(op_str) contains $(op_str[i]) which is not allowed")
        end
        # get index in order list of current char 
        idx::Int = findfirst(x -> x == op_str[i], order)
        # check if the previous index is <= than the current index
        if last_idx > idx
            error("Operator string $(op_str) is not ordered correctly. Order needs to be +, -, x, p, y, m, z.")
        end
        last_idx = idx
    end
end
function only_know_conjugate(op_str::String)
    do_conjugate, conjugate_operator = conjugate_to_basis(op_str)
    if do_conjugate
        error("Operator string $(op_str) is not in the basis, it's conjugate might be in the basis $(conjugate_operator)")
    end
end