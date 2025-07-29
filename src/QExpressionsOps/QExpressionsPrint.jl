using Printf
using LaTeXStrings
import Base: string

export string, latex_string, QExpr2string


function operators2string(op_indices::Vector{Is}, statespace::StateSpace; do_latex::Bool=false)::String
    op_str::String = ""
    subspaces = statespace.subspaces
    not_neutral = false
    i = 0
    for subspace in subspaces
        op_set = subspace.op_set
        for (ind, key) in zip(subspace.statespace_inds, subspace.keys)
            i += 1
            curr_op_ind = op_indices[i]
            if op_set.neutral_element != curr_op_ind
                not_neutral = true
                if do_latex
                    op_str *= op_set.op2latex(curr_op_ind, key)
                else
                    op_str *= op_set.op2str(curr_op_ind, key)
                end
            end
        end
    end
    return op_str 
end
function qAtom2string(qatom::QTerm, statespace::StateSpace; do_latex::Bool=false)::String
    return operators2string(qatom.op_indices, statespace; do_latex=do_latex)
end
function qAtom2string(q::QAbstract, statespace::StateSpace; do_latex::Bool=false)::String
    type = q.operator_type
    name = type.name 
    if do_latex 
        curr_string = raw"\hat{" * name * "}"
        if q.sub_index != -1 
            curr_string *= "_"*string(q.sub_index)
        end
        if q.dag
            curr_string *= raw"^{\dagger}"
            if q.exponent != 1 
                curr_string = raw"\left("*curr_string* raw"\right)^{"*string(q.exponent)*"}"
            end
        elseif q.exponent != 1 
            curr_string *= "^{"*string(q.exponent)*"}"
        end
    else
        curr_string = name 
        if q.sub_index != -1 
            curr_string *= str2sub(string(q.sub_index))
        end
        if q.dag 
            curr_string *= "'"
        end
        if q.exponent != 1 
            curr_string *= str2sup(string(q.exponent))
        end
    end
    return curr_string
end 

function variable_str_vec(q::T; do_latex::Bool=true)::Vector{String} where T<:QComposite
    if do_latex 
        return String[t.var_latex for t in q.statespace.vars]
    else
        return String[t.var_str for t in q.statespace.vars]
    end
end
function variable_str_vec(statespace::StateSpace; do_latex::Bool=true)::Vector{String}
    if do_latex 
        return String[t.var_latex for t in statespace.vars]
    else
        return String[t.var_str for t in statespace.vars]
    end
end

function sum_symbol_str(s::QSum; do_latex::Bool=false)
    if do_latex
        s_index_str = join(s.indexes, ",")
    else
        s_index_str = join(s.indexes, "")
    end
    n = length(s.indexes)
    equal_sign = "="
    if !s.neq
        if do_latex
            return "\\sum_{$s_index_str}^{=}"
        else
            return "∑" * str2sub(s_index_str) * "⁼"
        end
    else
        if n == 1
            if do_latex
                return "\\sum_{$s_index_str}^{\\neq}"
            else
                return "∑" * str2sub(s_index_str)
            end
        else
            if do_latex
                return "\\sum_{($s_index_str)}^{\\neq}"  # \\in \\mathcal{C}_N^{$n}
            else
                return "∑" * str2sub("(" * s_index_str * ")")
            end
        end
    end
end

# braced not used here as an argument use, it in other QComposites that contain QExpr to determine groupings!
function QComposite2string(q::QAtomProduct; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    if is_numeric(q)
        curr_sign, curr_str = to_stringer(q.coeff_fun, variable_str_vec(q, do_latex=do_latex), braced=false, do_frac=do_frac)
        if return_if_braced
            return curr_sign, curr_str, false
        else
            return curr_sign, curr_str
        end
    else
        curr_sign, curr_str = to_stringer(q.coeff_fun, variable_str_vec(q, do_latex=do_latex), braced=true, do_frac=do_frac, has_op=true)
        operator_str = join([qAtom2string(t, q.statespace, do_latex=do_latex) for t in q.expr], "")
        connector =do_latex ? raw" " : ""
        if return_if_braced
            return curr_sign, curr_str * connector * operator_str, false
        else
            return curr_sign, curr_str * connector * operator_str
        end
    end
end
function QComposite2string(term::QSum; do_latex::Bool=false, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    sum_str = sum_symbol_str(term, do_latex=do_latex)
    first_sign, total_string, single_group = QExpr2string(term.expr; do_latex=do_latex, braced=braced, do_frac=do_frac, return_grouping=true)
    if return_if_braced
        if single_group
            return first_sign, sum_str * total_string, false
        else
            return false, sum_str * brace(total_string, do_latex=do_latex), true
        end
    else
        if single_group
            return first_sign, sum_str * total_string
        else
            return false, sum_str * brace(total_string, do_latex=do_latex)
        end
    end
end
function QComposite2string(q::QCompositeProduct; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    total_sign = false
    all_strings::Vector{String} = []
    for term in q.expr
        new_sign, new_str, is_grouped = QExpr2string(term, do_latex=do_latex, braced=braced, do_frac=do_frac, return_grouping=true)
        total_sign = xor(total_sign, new_sign)
        if !is_grouped
            new_str = brace(new_str, do_latex=do_latex)
        end
        push!(all_strings, new_str)
    end
    # concatenate strings 
    connector = do_latex ? raw"\cdot" : "⋅"
    total_string = join(all_strings, connector) 
    if return_if_braced 
        return total_sign, total_string, true 
    else 
        return total_sign, total_string 
    end 
end

function do_return_braced_false(return_argument1::Bool, return_argument2::String, return_if_braced::Bool)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    if return_if_braced 
        return return_argument1, return_argument2, false
    else 
        return return_argument1, return_argument2
    end 
end

function QComposite2string(q::QExp; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    first_sign, total_str = QExpr2string(q.expr, do_latex=do_latex, braced=braced, do_frac=do_frac) 
    total_str = first_sign ? "-"*total_str : total_str 
    prefix = do_latex ? raw"\exp" : "exp"
    total_str = prefix * brace(total_str, do_latex=do_latex)
    return do_return_braced_false(false, total_str, return_if_braced) 
end
function QComposite2string(q::QLog; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    first_sign, total_str = QExpr2string(q.expr, do_latex=do_latex, braced=braced, do_frac=do_frac) 
    total_str = first_sign ? "-"*total_str : total_str 
    prefix = do_latex ? raw"\log" : "log"
    total_str = prefix * brace(total_str, do_latex=do_latex)
    return do_return_braced_false(false, total_str, return_if_braced) 
end
function QComposite2string(q::QPower; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    first_sign, total_str = QExpr2string(q.expr, do_latex=do_latex, braced=braced, do_frac=do_frac) 
    total_str = first_sign ? "-"*total_str : total_str 
    total_str = brace(total_str, do_latex=do_latex)
    if do_latex 
        return do_return_braced_false(false, total_str * "^{"*string(q.n)*"}", return_if_braced)
    else
        return do_return_braced_false(false, total_str * str2sup(string(q.n)), return_if_braced)
    end
end
function QComposite2string(q::QRoot; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    first_sign, total_str = QExpr2string(q.expr, do_latex=do_latex, braced=braced, do_frac=do_frac) 
    total_str = first_sign ? "-"*total_str : total_str 
    
    if do_latex 
        if q.n == 2
            return do_return_braced_false(false, raw"\sqrt{" * total_str * "}", return_if_braced)
        else
            return do_return_braced_false(false, raw"\sqrt["*string(q.n)*"]{" * total_str * "}", return_if_braced)
        end
    else
        total_str = brace(total_str, do_latex=do_latex)
        return do_return_braced_false(false, total_str * str2sup("1="*string(q.n)), return_if_braced)
    end
end
function QComposite2string(q::QCommutator; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    my_strings::Vector{String} = []
    for expr in q.expr
        first_sign, total_str = QExpr2string(expr, do_latex=do_latex, braced=braced, do_frac=do_frac)
        total_str = first_sign ? "-"*total_str : total_str 
        push!(my_strings, total_str)
    end
    if do_latex 
        return do_return_braced_false(false, raw"\left["*join(my_strings, raw",\,") * raw"\right]", return_if_braced)
    else
        return do_return_braced_false(false, "["*join(my_strings, ", ") * "]", return_if_braced)
    end
end

function QComposites2string(terms::AbstractVector{<: QComposite}; do_latex::Bool=false, braced::Bool=true, do_frac::Bool=true, separate_sign::Bool=false)::String
    substrings::Vector{Tuple{Bool, String}} = [QComposite2string(t, do_latex=do_latex, braced=braced, do_frac=do_frac) for t in terms]
    # connect substrings  
    if separate_sign
        string = substrings[1][2]
    else
        string = substrings[1][1] ? "-" * substrings[1][2] : substrings[1][2]
    end
    for (curr_sign, curr_str) in substrings[2:end]
        string *= curr_sign ? "-" * curr_str : "+" * curr_str
    end
    if separate_sign 
        return substrings[1][1], string
    end
    return string
end

import ..FFunctions: how_to_combine_Fs
function group_qAtomProducts(qs::Vector{QAtomProduct})::Vector{Union{QAtomProduct, Tuple{Union{FAtom, FSum}, Vector{QAtomProduct}}}}
    coeffs_funs::Vector{Union{FAtom, FSum}} = [q.coeff_fun for q in qs]
    coeff_groups, indexes = how_to_combine_Fs(coeffs_funs)
    new_qs = []
    for (coeffs, indexes) in zip(coeff_groups, indexes)
        if !( coeffs isa Tuple )
            push!(new_qs, copy(qs[indexes[1]]))
        else
            pre_F = coeffs[1]
            post_F = coeffs[2]
            post_q = QAtomProduct[]
            for (F, i) in zip(post_F, indexes) 
                new_p = copy(qs[i]) 
                new_p.coeff_fun = F 
                push!(post_q, new_p) 
            end
            push!(new_qs, (pre_F, post_q)) 
        end
    end 
    return new_qs 
end 

function qAtomProduct_group2string(qs::QAtomProduct; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true)::Tuple{Bool, String, String}
    curr_sign, operator_str = QComposite2string(qs, do_latex=do_latex, braced=braced, do_frac=do_frac)
    return curr_sign, "", operator_str
end
function qAtomProduct_group2string(qs::Tuple{Union{FAtom, FSum}, Vector{QAtomProduct}}; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true)::Tuple{Bool, String, String}
    # assume the qs can be simple grouped (see the functions: simple_combinable_Fs, group_Fs)
    F = qs[1]
    qs = qs[2]
    has_op = any([!is_numeric(s) for s in qs])
    f_sign, f_str = to_stringer(F, variable_str_vec(qs[1], do_latex=do_latex); do_latex=do_latex, braced=braced, do_frac=do_frac, has_op=has_op)
    q_str = QComposites2string(qs; do_latex=do_latex, braced=braced, do_frac=do_frac, separate_sign=false)  # don't worry about internal signs, this has already been taken care off by the sign handling of the grouping 
    return f_sign, f_str, q_str
end 

function allnegative(x::Tuple{Bool, String})::Bool 
    return x[1]
end
function allnegative(x::Vector{Tuple{Bool, String}})::Bool
    return all(allnegative, x)
end
function brace(x::String; do_latex::Bool=true)::String
    if do_latex 
        return raw"\left(" * x * raw"\right)"
    else
        return "(" * x * ")"
    end
end
function QExpr2string(q::QExpr; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_grouping::Bool=false)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    # outputs sign, string, {optional return_grouping:} single_group::Bool   => return grouping implies that the expression will be braced if it isn't already! , hence the outputted sign is handled differently 
    # first sort terms 
    # q_sorted = q   # activate if debugging for checking if sorting is the problem 
    # q_sorted = sort(q)
    q_sorted = simplify(q)
    if !braced # outside (not inside of a QComposite)
        if !return_grouping
            return QComposites2string(q_sorted.terms, do_latex=do_latex, braced=braced, do_frac=do_frac, separate_sign=true)
        else
            first_sign, total_string = QComposites2string(q_sorted.terms, do_latex=do_latex, braced=braced, do_frac=do_frac, separate_sign=true)
            return first_sign, total_string, false
        end
    else  # inside another QComposite
        # separate into QAtomProduct and other QComposite terms -> sorting puts qAtomProducts first 
        # then group qAtomProducts by their factors 
        first_non_qAtomProduct = findfirst(x -> !isa(x, QAtomProduct), q_sorted.terms)
        if first_non_qAtomProduct === nothing
            first_non_qAtomProduct = length(q_sorted.terms) + 1
        end
        qAtomProduct_terms::Vector{QAtomProduct} = q_sorted.terms[1:first_non_qAtomProduct-1]
        other_terms = q_sorted.terms[first_non_qAtomProduct:end]
        groups = group_qAtomProducts(qAtomProduct_terms)
        # create strings for each element 
        all_strings::Vector{Tuple{Bool, String}} = []
        for group in groups
            #if isonelike(group[1])  # => remove grouping with brace 
            curr_sign, first, second = qAtomProduct_group2string(group, do_latex=do_latex, braced=braced, do_frac=do_frac)
            if isa(group, Tuple)
                push!(all_strings, (curr_sign, first*brace(second, do_latex=do_latex)))
            else
                push!(all_strings, (curr_sign, first*second))
            end
        end
        for term in other_terms
            push!(all_strings, QComposite2string(term, do_latex=do_latex, braced=braced, do_frac=do_frac))
        end
        
        # make QComposites2string better, by allowing a third output in case of single group. to traverse multiple layers of Composites within composites. 
        if return_grouping 
            # do we switch the sign? 
            if allnegative(all_strings) || ( get_default(:FLIP_IF_FIRST_TERM_NEGATIVE) && all_strings[1][1] )
                # switch signs 
                first_sign = true 
                all_strings = [(!sign, s) for (sign, s) in all_strings]
            else
                first_sign = false
            end
            if length(groups) == 1 && length(all_strings) == 1
                total_string = all_strings[1][2]
            else
                total_string = all_strings[1][1] ? "-" : ""
                total_string *= all_strings[1][2]
            end
        else
            first_sign = all_strings[1][1]
            total_string = all_strings[1][2]
        end
        for (sign, s) in all_strings[2:end]
            if sign 
                total_string *= "-" * s
            else 
                total_string *= "+" * s
            end
        end
        if return_grouping 
            single_group::Bool = length(all_strings) == 1   # length(groups) == 1 || 
            return first_sign, total_string, single_group
        else
            return first_sign, total_string
        end
    end
end

function diff_qEQ2string(eq::diff_QEq; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true)::String
    curr_sign, curr_string = QExpr2string(eq.expr, do_latex=do_latex, braced=braced, do_frac=do_frac)
    right_hand_side = curr_sign ? "-" * curr_string : curr_string

    curr_sign, curr_string = QComposite2string(eq.left_hand_side, do_latex=do_latex, braced=braced, do_frac=do_frac)
    left_hand_side_op_str = curr_sign ? "-" * curr_string : curr_string
    left_hand_side_op_str = lstrip(left_hand_side_op_str, '+')

    if do_latex
        #left_hand_side = raw"\frac{\text{d} \phantom{t}}{\text{d}t}\!" * left_hand_side_op_str
        left_hand_side = raw"\frac{\text{d} "*left_hand_side_op_str*raw"}{\text{d}t}"
    else
        left_hand_side = "d(" * left_hand_side_op_str * ") / dt"
    end
    return left_hand_side * " = " * right_hand_side
end

#### String ##########################################################################################################################
""" 
    string(eq::QExpr) -> String
    string(eq::QAtomProduct) -> String
    string(eq::diff_QEq) -> String

Returns a string representation of the QExpr, QAtomProduct or diff_QEq object. The string is formatted in a human-readable way, but without LaTeX formatting.
"""
function string(eq::QExpr)::String
    # add default variables for do_Frac, braced and so on. take care of this by writing a single function called by every string and latex string function 
    curr_sign, curr_string = QExpr2string(eq, do_latex=false, braced=get_default(:DO_BRACED))
    total_string = curr_sign ? "-" * curr_string : curr_string
    return total_string
end
function string(eq::QAtomProduct)::String
    # add default variables for do_Frac,, braced and so on. take care of this by writing a single function called by every string and latex string function 
    sign, total_string = QComposite2string(eq; do_latex=false, braced=get_default(:DO_BRACED))
    total_string = sign ? "-" * total_string : total_string
    return total_string
end
function string(eq::diff_QEq)::String
    return diff_qEQ2string(eq; do_latex=false, braced=get_default(:DO_BRACED))
end

#### LaTeX-String ##########################################################################################################################
""" 
    latex_string(eq::QExpr) -> String
    latex_string(eq::QAtomProduct) -> String
    latex_string(eq::diff_QEq) -> String

Returns a LaTeX string representation of the QExpr, QAtomProduct or diff_QEq object. 
"""
function latex_string(eq::QExpr)::String
    curr_sign, curr_string = QExpr2string(eq, do_latex=true, braced=get_default(:DO_BRACED))
    total_string = curr_sign ? "-" * curr_string : curr_string
    return total_string
end
function latex_string(eq::QAtomProduct)::String
    # add default variables for do_Frac, braced and so on. take care of this by writing a single function called by every string and latex string function 
    sign, total_string = QComposite2string(eq; do_latex=true, braced=get_default(:DO_BRACED))
    total_string = sign ? "-" * total_string : total_string
    return total_string
end
function latex_string(eq::diff_QEq)::String
    return diff_qEQ2string(eq; do_latex=true, braced=get_default(:DO_BRACED))
end

#### Show off #########################################################################################################################
function show(io::IO, x::QExpr)
    print(io, string(x))
end
function show(io::IO, ::MIME"text/latex", x::QExpr)
    print(io, latexstring(latex_string(x)))
end

function show(io::IO, x::QAtomProduct)
    print(io, string(x))
end
function show(io::IO, ::MIME"text/latex", x::QAtomProduct)
    print(io, latexstring(latex_string(x)))
end

function show(io::IO, q::diff_QEq)
    print(io, string(q))
end
function show(io::IO, ::MIME"text/latex", q::diff_QEq)
    print(io, latexstring(latex_string(q)))
end

