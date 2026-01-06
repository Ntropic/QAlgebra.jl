using Printf
using LaTeXStrings
import Base: string
using ..QSpaces: AbstractIndex2string

export string, latex_string, cumulant_string

const EnsembleBits = Vector{BitSet}

include("QExpressions_print_sum.jl")


function _compact_cumulant_string(q::QCumulant; do_latex::Bool)
    coeff_sign, coeff_str = to_stringer(q.coeff_fun; do_latex=do_latex, braced=DO_BRACED, do_frac=true, has_op=true)
    atom_product = QAtomProduct(q.qspace, q.atom, false)
    atom_sign, atom_body = QComposite2string(atom_product, which_ensemble_acting(atom_product); do_latex=do_latex, braced=true, do_frac=true, do_braket=true)
    total_sign = xor(coeff_sign, atom_sign)

    body = ""
    if !isempty(coeff_str)
        body *= coeff_str
        if !isempty(atom_body)
            body *= (do_latex ? raw" " : " ") * atom_body
        end
    else
        body = atom_body
    end
    body = isempty(body) ? (do_latex ? "0" : "0") : body
    suffix = do_latex ? "_{C}" : str2sup("c")
    body *= suffix
    return total_sign, body
end

function _expanded_cumulant_string(q::QCumulant; do_latex::Bool)
    lhs_product = QAtomProduct(q.qspace, q.coeff_fun, [q.atom], false)
    lhs_sign, lhs_body = QComposite2string(lhs_product, which_ensemble_acting(lhs_product, do_abstract=true); do_latex=do_latex, braced=DO_BRACED, do_frac=true, do_braket=true)
    lhs_display = lhs_sign ? "-" * lhs_body : lhs_body
    lhs_display = isempty(lhs_display) ? (do_latex ? "0" : "0") : lhs_display

    rhs_sign, rhs_body = QExpr2string(q.expr, which_ensemble_acting(q.expr, do_abstract=true); do_latex=do_latex, braced=DO_BRACED, do_braket=true)
    rhs_display = rhs_sign ? "-" * rhs_body : rhs_body
    rhs_display = isempty(rhs_display) ? (do_latex ? "0" : "0") : rhs_display

    approx_symbol = do_latex ? raw" \approx " : " ≈ "
    return lhs_display * approx_symbol * rhs_display
end


function qAtom2string(q::QTerm, qspace::QSpace; do_latex::Bool=false)::String
    op_indices = q.op_indices
    op_str::String = ""
    subspaces = qspace.subspaces
    for particle in op_indices
        idx = particle.index
        subspace = subspaces[idx.subspace]
        op_set = subspace.op_set
        particle.operator == op_set.neutral_element && continue
        add_index = idx.ensemble != 0 || PRINT_NON_ENSEMBLE_INDEXES
        if add_index
            if idx.ensemble == 0
                label = subspace.key
            elseif do_latex
                label = AbstractIndex2string(qspace.subspace_info, idx; do_latex=true, as_index=false)
            else
                ensemble = subspace.ensemble
                index_symbol = idx.summation ? ensemble.sum_string : ensemble.non_sum_string
                subindex = idx.slot == 0 ? "" : string(idx.slot)
                label = index_symbol * subindex
            end
        else
            label = ""
        end
        if do_latex
            op_str *= op_set.op2latex(particle.operator, label; add_index=add_index)
        else
            op_str *= op_set.op2str(particle.operator, label; add_index=add_index)
        end
    end
    if q.time_index.order != -1
        op_str *= "("*t_suffix(q.time_index.order, do_latex=do_latex)*")"
    end
    return op_str 
end

function qAtom2string(q::QAbstract, qspace::QSpace; do_latex::Bool=false)::String
    type = q.operator_type
    name = type.name 
    if do_latex 
        op_str = raw"\hat{" * name * "}"
        if q.sub_index != -1 
            op_str *= "_"*string(q.sub_index)
        end
        if q.dag
            if q.exponent == 1 
                op_str *= raw"^{\dagger}"
            else
                op_str *= raw"^{\dagger "*string(q.exponent)*"}"
            end
        elseif q.exponent != 1 
            op_str *= "^{"*string(q.exponent)*"}"
        end
    else
        op_str = name 
        if q.sub_index != -1 
            op_str *= str2sub(string(q.sub_index))
        end
        if q.dag 
            op_str *= "'"
        end
        if q.exponent != 1 
            op_str *= str2sup(string(q.exponent))
        end
    end
    if q.time_index != -1
        op_str *= "("*t_suffix(q.time_index, do_latex=do_latex)*")"
    end
    return op_str
end 


# braced not used here as an argument use, it in other QComposites that contain QExpr to determine groupings!
function QComposite2string(q::QAtomProduct, ::EnsembleBits; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false, do_braket::Bool=false)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    if isnumeric(q)
        curr_sign, curr_str = to_stringer(q.coeff_fun, braced=false, do_frac=do_frac, do_latex=do_latex)
        if return_if_braced
            return curr_sign, curr_str, length(q.coeff_fun) == 1
        else
            return curr_sign, curr_str
        end
    else
        curr_sign, curr_str = to_stringer(q.coeff_fun, braced=true, do_frac=do_frac, has_op=true, do_latex=do_latex)
        show_braket = do_braket || q.braket
        if !show_braket
                operator_str = join([qAtom2string(t, q.qspace, do_latex=do_latex) for t in q.expr], "")
        else
            if q.separate_expectation_values
                operator_str = join([braket(qAtom2string(t, q.qspace, do_latex=do_latex), do_latex=do_latex) for t in q.expr], "")
            else
                operator_str = braket(join([qAtom2string(t, q.qspace, do_latex=do_latex) for t in q.expr], ""), do_latex=do_latex)
            end
        end
        connector =do_latex ? raw" " : ""
        if return_if_braced
            return curr_sign, curr_str * connector * operator_str, false
        else
            return curr_sign, curr_str * connector * operator_str
        end
    end
end

function QComposite2string(q::QAtomIndexed, where_acting::EnsembleBits; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false, do_braket::Bool=false)
    coeff_sign, coeff_str = to_stringer(q.coeff_fun; braced=true, do_frac=do_frac, has_op=true, do_latex=do_latex)
    full_indices = recompose_op_indices(q.op_indices, q.ensemble_indices, q.qspace)
    term = QTerm(full_indices, q.time_index)
    operator_str = qAtom2string(term, q.qspace; do_latex=do_latex)

    flat_indices = [string(idx) for ensemble in q.concrete_indices.indices for idx in ensemble]
    if !isempty(flat_indices)
        operator_str *= indices2str(flat_indices; do_latex=do_latex)
    end

    connector = do_latex ? raw" " : ""
    total = coeff_str * connector * operator_str
    if return_if_braced
        return coeff_sign, total, false
    else
        return coeff_sign, total
    end
end
function QComposite2string(term::AbstractQSum, where_acting::EnsembleBits; do_latex::Bool=false, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false, do_braket::Bool=false)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    sum_str, measure = sum_symbol_str(term, where_acting; do_latex=do_latex)
    first_sign, total_string, single_group = QExpr2string(term.expr, where_acting; do_latex=do_latex, braced=braced, do_frac=do_frac, return_grouping=true, do_braket=do_braket)
    measure_str = measure
    if return_if_braced
        if single_group
            return first_sign, sum_str * total_string * measure_str, false
        else
            return false, sum_str * brace(total_string, do_latex=do_latex) * measure_str, true
        end
    else
        if single_group
            return first_sign, sum_str * total_string * measure_str
        else
            return false, sum_str * brace(total_string, do_latex=do_latex) * measure_str
        end
    end
end
function QComposite2string(q::QCompositeProduct, where_acting::EnsembleBits; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false, do_braket::Bool=false)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    parts = String[]
    total_sign, coeff_str = to_stringer(q.coeff_fun; braced=true, do_frac=do_frac, has_op=true, do_latex=do_latex)
    if !isempty(coeff_str)
        push!(parts, coeff_str)
    end
    for term in q.expr
        term_sign, term_str, _ = QComposite2string(term, where_acting; do_latex=do_latex, braced=braced, do_frac=do_frac, return_if_braced=true, do_braket=do_braket)
        total_sign = xor(total_sign, term_sign)
        if !isempty(term_str)
            push!(parts, term_str)
        end
    end

    connector = do_latex ? raw" " : " "
    total_string = join(parts, connector)
    if return_if_braced
        return total_sign, total_string, true
    else
        return total_sign, total_string
    end
end
function QComposite2string(q::QCommutator, where_acting::EnsembleBits; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false, do_braket::Bool=false)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    all_strings::Vector{String} = []
    total_sign, coeff_str = to_stringer(q.coeff_fun; braced=true, do_frac=do_frac, has_op=true, do_latex=do_latex)
    for expr in q.expr
        first_sign, total_str = QExpr2string(expr, where_acting; do_latex=do_latex, braced=braced, do_frac=do_frac, do_braket=do_braket)
        total_str = first_sign ? "-"*total_str : total_str 
        push!(all_strings, total_str)
    end
    if do_latex 
        return do_return_braced_true(total_sign, coeff_str*raw"\left["*join(all_strings, raw",\,") * raw"\right]", return_if_braced)
    else
        return do_return_braced_true(total_sign, coeff_str*"["*join(all_strings, ", ") * "]", return_if_braced)
    end
end

function QComposite2string(q::QCumulant, ::EnsembleBits; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false, do_braket::Bool=false)
    if EXPAND_CUMULANTS
        body = _expanded_cumulant_string(q; do_latex=do_latex)
        if return_if_braced
            return false, body, false
        else
            return false, body
        end
    else
        sign, body = _compact_cumulant_string(q; do_latex=do_latex)
        if return_if_braced
            return sign, body, false
        else
            return sign, body
        end
    end
end


wrap(q::QExp, s::String; do_latex) = (do_latex ? raw"\exp" : "exp") * brace(s; do_latex=do_latex)
wrap(q::QLog, s::String; do_latex) = (do_latex ? raw"\log" : "log") * brace(s; do_latex=do_latex)
wrap(q::QPower, s; do_latex) = begin
    base = brace(s; do_latex=do_latex)
    do_latex ? base * "^{" * string(q.n) * "}" : base * str2sup(string(q.n))
end
function wrap(q::QRoot, s; do_latex)
    if do_latex
        q.n == 2 ? (raw"\sqrt{" * s * "}") : (raw"\sqrt[" * string(q.n) * "]{" * s * "}")
    else
        brace(s; do_latex=do_latex) * str2sup("1/" * string(q.n))
    end
end

@inline _ret(sign::Bool, s::String, return_if_braced::Bool; grouped::Bool) = return_if_braced ? (sign, s, grouped) : (sign, s)


function QComposite2string(q::T, where_acting::EnsembleBits; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_if_braced::Bool=false, do_braket::Bool=false)::Union{Tuple{Bool,String},Tuple{Bool,String,Bool}} where T <: QComposite
    # First: coefficient
    coeff_sign, coeff_str = to_stringer(q.coeff_fun; braced=true, do_frac=do_frac, has_op=true, do_latex=do_latex)

    # Then: inner expression
    inner_sign, inner_str = QExpr2string(q.expr, where_acting; do_latex=do_latex, braced=braced, do_frac=do_frac, do_braket=do_braket)
    inner_str = inner_sign ? "-" * inner_str : inner_str

    # Wrap (exp, log, power, root, etc.)
    wrapped = wrap(q, inner_str; do_latex=do_latex)

    total_str = coeff_str * " " * wrapped
    return _ret(coeff_sign, total_str, return_if_braced; grouped=false)
end


function QComposites2string(terms::AbstractVector{<: QComposite}, where_acting::EnsembleBits; do_latex::Bool=false, braced::Bool=true, do_frac::Bool=true, separate_sign::Bool=false, do_braket::Bool=false)::String
    substrings::Vector{Tuple{Bool, String}} = [QComposite2string(t, where_acting; do_latex=do_latex, braced=braced, do_frac=do_frac, do_braket=do_braket) for t in terms]
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

import ..CFunctions: how_to_combine_Fs
function group_qAtomProducts(qs::Vector{QAtomProduct})::Vector{Union{QAtomProduct, Tuple{Union{CAtom, CSum}, Vector{QAtomProduct}}}}
    coeffs_funs::Vector{CFunction} = [q.coeff_fun for q in qs]
    coeff_groups, indices = how_to_combine_Fs(coeffs_funs)
    new_qs = []
    for (coeffs, indices) in zip(coeff_groups, indices)
        if !( coeffs isa Tuple )
            push!(new_qs, qs[indices[1]])
        else
            pre_F = coeffs[1]
            post_F = coeffs[2]
            post_q = QAtomProduct[]
            for (F, i) in zip(post_F, indices) 
                push!(post_q, modify_coeff(qs[i], F)) 
            end
            push!(new_qs, (pre_F, post_q)) 
        end
    end 
    return new_qs 
end 

function qAtomProduct_group2string(qs::QAtomProduct, where_acting::EnsembleBits; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, do_braket::Bool=false)::Tuple{Bool, String, String}
    curr_sign, operator_str = QComposite2string(qs, where_acting; do_latex=do_latex, braced=braced, do_frac=do_frac, do_braket=do_braket)
    return curr_sign, "", operator_str
end
function qAtomProduct_group2string(qs::Tuple{Union{CAtom, CSum}, Vector{QAtomProduct}}, where_acting::EnsembleBits; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, do_braket::Bool=false)::Tuple{Bool, String, String}
    # assume the qs can be simple grouped (see the functions: simple_combinable_Fs, group_Fs)
    F = qs[1]
    qs = qs[2]
    has_op = any([!isnumeric(s) for s in qs])
    f_sign, f_str = to_stringer(F; do_latex=do_latex, braced=braced, do_frac=do_frac, has_op=has_op)
    q_str = QComposites2string(qs, where_acting; do_latex=do_latex, braced=braced, do_frac=do_frac, separate_sign=false, do_braket=do_braket)  # don't worry about internal signs, this has already been taken care off by the sign handling of the grouping 
    return f_sign, f_str, q_str
end 

function allnegative(x::Tuple{Bool, String})::Bool 
    return x[1]
end
function allnegative(x::Vector{Tuple{Bool, String}})::Bool
    return all(allnegative, x)
end
function QExpr2string(q::QExpr, where_acting::EnsembleBits; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, return_grouping::Bool=false, do_braket::Bool=false)::Union{Tuple{Bool, String}, Tuple{Bool, String, Bool}}
    # outputs sign, string, {optional return_grouping:} single_group::Bool   => return grouping implies that the expression will be braced if it isn't already! , hence the outputted sign is handled differently 
    sort!(q)  # sorts by term
    if !braced # outside (not inside of a QComposite)
        if !return_grouping
            return QComposites2string(q.terms, where_acting; do_latex=do_latex, braced=braced, do_frac=do_frac, separate_sign=true, do_braket=do_braket)
        else
            first_sign, total_string = QComposites2string(q.terms, where_acting; do_latex=do_latex, braced=braced, do_frac=do_frac, separate_sign=true, do_braket=do_braket)
            return first_sign, total_string, false
        end
    else  # inside another QComposite
        # separate into QAtomProduct and other QComposite terms -> sorting puts qAtomProducts first 
        # then group qAtomProducts by their factors 
        first_non_qAtomProduct = findfirst(x -> !isa(x, QAtomProduct), q.terms)
        if first_non_qAtomProduct === nothing
            first_non_qAtomProduct = length(q.terms) + 1
        end
        qAtomProduct_terms::Vector{QAtomProduct} = q.terms[1:first_non_qAtomProduct-1]
        other_terms = q.terms[first_non_qAtomProduct:end]
        groups = group_qAtomProducts(qAtomProduct_terms)
        # create strings for each element 
        all_strings::Vector{Tuple{Bool, String}} = []
        for group in groups
            #if isonelike(group[1])  # => remove grouping with brace 
            curr_sign, first, second = qAtomProduct_group2string(group, where_acting; do_latex=do_latex, braced=braced, do_frac=do_frac, do_braket=do_braket)
            if isa(group, Tuple)
                push!(all_strings, (curr_sign, first*brace(second, do_latex=do_latex)))
            else
                push!(all_strings, (curr_sign, first*second))
            end
        end
        is_braced::Bool = false
        for term in other_terms
            curr_sign, curr_str, is_braced = QComposite2string(term, where_acting; do_latex=do_latex, braced=braced, do_frac=do_frac, return_if_braced=true, do_braket=do_braket)
            push!(all_strings, (curr_sign, curr_str))
        end
        
        # make QComposites2string better, by allowing a third output in case of single group. to traverse multiple layers of Composites within composites. 
        if return_grouping 
            # do we switch the sign? 
            if allnegative(all_strings) || (FLIP_IF_FIRST_TERM_NEGATIVE  && all_strings[1][1] )
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
            if single_group && length(groups) == 0
                # is other_term 
                single_group = true 
            end
            return first_sign, total_string, single_group
        else
            return first_sign, total_string
        end
    end
end

function diffQEq2string(eq::diffQEq; do_latex::Bool=true, braced::Bool=true, do_frac::Bool=true, do_braket::Bool=false)::String
    rhs_where_acting = which_ensemble_acting(eq.expr, do_abstract=true)
    curr_sign, curr_string = QExpr2string(eq.expr, rhs_where_acting; do_latex=do_latex, braced=braced, do_frac=do_frac, do_braket=do_braket)
    right_hand_side = curr_sign ? "-" * curr_string : curr_string

    lhs_braket = do_braket || eq.left_hand_side.braket
    lhs_where_acting = which_ensemble_acting(eq.left_hand_side, do_abstract=true)
    curr_sign, curr_string = QComposite2string(eq.left_hand_side, lhs_where_acting; do_latex=do_latex, braced=braced, do_frac=do_frac, do_braket=lhs_braket)
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
    string(eq::diffQEq) -> String

Returns a string representation of the QExpr, QAtomProduct or diffQEq object. The string is formatted in a human-readable way, but without LaTeX formatting.
"""
function string(eq::QExpr)::String
    # add default variables for do_Frac, braced and so on. take care of this by writing a single function called by every string and latex string function 
    where_acting = which_ensemble_acting(eq, do_abstract=true)
    curr_sign, curr_string = QExpr2string(eq, where_acting; do_latex=false, braced=DO_BRACED)
    total_string = curr_sign ? "-" * curr_string : curr_string
    return total_string
end
function string(eq::QCumulant)::String
    where_acting = which_ensemble_acting(eq, do_abstract=true)
    sign, total_string = QComposite2string(eq, where_acting; do_latex=false, braced=DO_BRACED)
    return sign ? "-" * total_string : total_string
end
function string(eq::QAtomProduct)::String
    # add default variables for do_Frac, braced and so on. take care of this by writing a single function called by every string and latex string function 
    where_acting = which_ensemble_acting(eq, do_abstract=true)
    sign, total_string = QComposite2string(eq, where_acting; do_latex=false, braced=DO_BRACED)
    total_string = sign ? "-" * total_string : total_string
    return total_string
end
function string(eq::diffQEq)::String
    return diffQEq2string(eq; do_latex=false, braced=DO_BRACED)
end

#### LaTeX-String ##########################################################################################################################
""" 
    latex_string(eq::QExpr) -> String
    latex_string(eq::QAtomProduct) -> String
    latex_string(eq::diffQEq) -> String

Returns a LaTeX string representation of the QExpr, QAtomProduct or diffQEq object. 
"""
function latex_string(eq::QExpr)::String
    where_acting = which_ensemble_acting(eq, do_abstract=true)
    curr_sign, curr_string = QExpr2string(eq, where_acting; do_latex=true, braced=DO_BRACED)
    total_string = curr_sign ? "-" * curr_string : curr_string
    return total_string
end
function latex_string(eq::QCumulant)::String
    where_acting = which_ensemble_acting(eq, do_abstract=true)
    sign, total_string = QComposite2string(eq, where_acting; do_latex=true, braced=DO_BRACED)
    return sign ? "-" * total_string : total_string
end
function latex_string(eq::QAtomProduct)::String
    # add default variables for do_Frac, braced and so on. take care of this by writing a single function called by every string and latex string function 
    where_acting = which_ensemble_acting(eq, do_abstract=true)
    sign, total_string = QComposite2string(eq, where_acting; do_latex=true, braced=DO_BRACED)
    total_string = sign ? "-" * total_string : total_string
    return total_string
end
function latex_string(eq::diffQEq)::String
    return diffQEq2string(eq; do_latex=true, braced=DO_BRACED)
end

#### Show off #########################################################################################################################
function show(io::IO, x::QExpr)
    print(io, string(x))
end
function show(io::IO, ::MIME"text/latex", x::QExpr)
    print(io, latexstring(latex_string(x)))
end

function show(io::IO, x::QCumulant)
    print(io, string(x))
end
function show(io::IO, ::MIME"text/latex", x::QCumulant)
    print(io, latexstring(latex_string(x)))
end

"""
    cumulant_string(q::QCumulant; do_latex=false, expanded=EXPAND_CUMULANTS) -> String 
    cumulant_string(q::QExpr; do_latex=false, expanded=EXPAND_CUMULANTS) -> String 

Prints the cumulant approximation (expanded or compact) and the operator it approximates. 
"""
function cumulant_string(q::QCumulant; do_latex::Bool=false, expanded::Bool=EXPAND_CUMULANTS)
    if expanded
        return _expanded_cumulant_string(q; do_latex=do_latex)
    else
        sign, body = _compact_cumulant_string(q; do_latex=do_latex)
        return sign ? "-" * body : body
    end
end

function cumulant_string(q::QExpr; do_latex::Bool=false, expanded::Bool=EXPAND_CUMULANTS)
    length(q.terms) == 1 || error("cumulant_string expects a QExpr with exactly one term; got $(length(q.terms)).")
    term = q.terms[1]
    term isa QCumulant || error("cumulant_string only supports QExpr whose single term is a QCumulant; got $(typeof(term)).")
    return cumulant_string(term; do_latex=do_latex, expanded=expanded)
end

function show(io::IO, x::QAtomProduct)
    print(io, string(x))
end
function show(io::IO, ::MIME"text/latex", x::QAtomProduct)
    print(io, latexstring(latex_string(x)))
end

function show(io::IO, q::diffQEq)
    print(io, string(q))
end
function show(io::IO, ::MIME"text/latex", q::diffQEq)
    print(io, latexstring(latex_string(q)))
end
