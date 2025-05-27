using Printf
using LaTeXStrings
import Base: string

export string, latex_string


function var_exponents2string(var_exponents::Vector{Int}, statespace::StateSpace; do_latex::Bool=false)::String
    monomial::String = ""
    for (i, exp) in enumerate(var_exponents)
        if exp != 0
            # Convert the variable to a string; you may adjust formatting if needed.
            if do_latex
                var = statespace.vars[i].var_latex
            else
                var = statespace.vars[i].var_str
            end
            monomial *= var
            if exp != 1
                if do_latex
                    monomial *= "^" * string(exp)
                else
                    monomial *= str2sup(string(exp))
                end
            end
        end
    end
    monomial = strip(monomial)
    return monomial
end

function operators2string(op_indices::Vector{Is}, statespace::StateSpace; do_latex::Bool=false, do_sigma::Bool=false)::String
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
                    op_str *= op_set.op2latex(curr_op_ind, key, do_sigma)
                else
                    op_str *= op_set.op2str(curr_op_ind, key)
                end
            end
        end
    end
    return op_str 
end

is_negative(x::Real) = x < 0
function is_negative(x::Complex)
    if !iszero(real(x))
        return real(x) < 0
    elseif !iszero(imag(x))
        return imag(x) < 0
    else
        return false
    end
end
function are_negative(x::Vector)
    for i in x
        if !is_negative(i)
            return false
        end
    end
    return true
end
function format_number(num::Number; precision::Int=4, n::Int=2, do_latex::Bool=false)
    # Special-case zero:
    if iszero(num)
        return "0"
    end

    # First, get a scientific notation string.
    sci_str = @sprintf("%.*e", precision, num)  # e.g. "1.2345e+03"
    parts = split(sci_str, "e")
    significand = parts[1]
    exponent_str = parts[2]
    exponent = parse(Int, exponent_str)

    # If exponent is within [-n, n], return fixed-point formatting.
    if -n <= exponent <= n
        return @sprintf("%.*f", precision, num)
    else
        if do_latex
            return significand * " \\times 10^{" * string(exponent) * "}"
        else
            return significand * "×10" * str2sup(string(exponent))
        end
    end
end
function format_number(num::Complex; precision::Int=4, n::Int=2, do_latex::Bool=false)
    re_val = real(num)
    im_val = imag(num)
    # If both parts are zero:
    if iszero(re_val) && iszero(im_val)
        return "0"
        # Only imaginary part nonzero:
    elseif iszero(re_val)
        im_str = format_number(abs(im_val); precision=precision, n=n, do_latex=do_latex)
        # Prepend a minus sign if needed.
        return (im_val < 0 ? "-" : "") * im_str * "i"
        # Only real part nonzero:
    elseif iszero(im_val)
        return format_number(re_val; precision=precision, n=n, do_latex=do_latex)
    else
        # Both parts nonzero: format each and wrap in parentheses.
        re_str = format_number(re_val; precision=precision, n=n, do_latex=do_latex)
        im_str = format_number(abs(im_val); precision=precision, n=n, do_latex=do_latex)
        # Use a minus sign if imaginary part is negative.
        sign_str = im_val < 0 ? " - " : " + "
        return "(" * re_str * sign_str * im_str * "i" * ")"
    end
end


function prefactor_to_string(coeff::Number, var_exponents::Vector{Int}, statespace::StateSpace, any_ops::Bool=false; do_latex::Bool=false)::Tuple{String,Bool}
    any_vars = any(var_exponents .!= 0)
    any_either = any_ops || any_vars
    varstr = var_exponents2string(var_exponents, statespace, do_latex=do_latex)
    coeff = crationalize(coeff)
    coeffstr = ""
    coeff_one = false
    if isone(abs(coeff)) && (isone(abs(real(coeff))) || isone(abs(imag(coeff))))
        coeff_one = true
        if isone(abs(imag(coeff)))
            if coeff.b / coeff.c > 0
                coeffstr *= "+i"
            else
                coeffstr *= "-i"
            end
        else
            if coeff.a / coeff.c > 0
                coeffstr *= "+"
            else
                coeffstr *= "-"
            end
            if !any_either
                coeffstr *= "1"
            end
        end
    elseif isa(coeff, ComplexRational)
        coeffstr *= complexrational2str(coeff, do_latex)
    else
        # use get_default(:FLOAT_DIGITS) to determine the digits of precision to use for the coefficient
        precision = get_default(:FLOAT_DIGITS)
        n = get_default(:EXP_DIGITS)
        coeffstr *= format_number(coeff, precision=precision, n=n, do_latex=do_latex)
    end
    if length(varstr) > 0
        if do_latex
            varstr *= raw"\,"
        else
            varstr *= "\u200A"
        end
    elseif length(coeffstr) > 0
        if do_latex
            coeffstr *= raw"\,"
        else
            coeffstr *= "\u200A"
        end
    end
    if !do_latex && any_vars && !coeff_one
        coeffstr *= "×"
    end
    final_str = strip(coeffstr * varstr)
    if !(startswith(final_str, "+") || startswith(final_str, "-"))
        final_str = "+" * final_str
    end
    do_either = any_vars || !coeff_one
    return final_str, do_either
end
function qobj_to_string(t::qTerm, statespace::StateSpace; do_latex::Bool=false, do_sigma::Bool=false, braket::Bool=false)::String
    # Build the monomial part
    var_exponents = t.var_exponents
    op_indices = t.op_indices
    op_str = operators2string(op_indices, statespace, do_latex=do_latex, do_sigma=do_sigma)
    any_ops = length(op_str) > 0
    prefactor_string, _ = prefactor_to_string(t.coeff, var_exponents, statespace, any_ops, do_latex=do_latex)
    if braket && any_ops
        if do_latex
            op_str = raw"\braket{" * op_str * raw"}"
        else
            op_str = "⟨" * op_str * "⟩"
        end
    end
    final_str = strip(prefactor_string * op_str)
    return final_str
end
function sum_symbol_str(s::qSum; do_latex::Bool=false)
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
function qobj_to_string(s::qSum, statespace::StateSpace; do_latex::Bool=false, do_sigma::Bool=false, braket::Bool=false)::String
    sum_str = sum_symbol_str(s; do_latex=do_latex)
    grouped_expressions = separate_into_common_subterms(s.expr.terms, statespace)
    allnegative = are_negative([g[1] for g in grouped_expressions])
    if allnegative
        new_grouped_expressions = [(-g[1], g[2], g[3]) for g in grouped_expressions]
        grouped_expressions = new_grouped_expressions
    end
    sep_strings = [brace(t, statespace, do_latex=do_latex, do_sigma=do_sigma, braket=braket) for t in grouped_expressions]
    term_string = strip(join(sep_strings))
    term_string = lstrip(term_string, '+')
    if length(grouped_expressions) > 1
        if do_latex
            term_string = raw"\left(" * term_string * raw"\right)"
        else
            term_string = "(" * term_string * ")"
        end
    end
    if !allnegative
        return "+" * sum_str * " " * term_string
    else
        return "-" * sum_str * " " * term_string
    end
end

"""
    best_common_factor(coeffs::Vector{Number}) -> Number

Given a vector of nonzero coefficients, try every candidate (each unique coefficient)
and return the one that, when used as a divisor for all coefficients, produces the most ones.
If there is a tie, the candidate with the lowest total absolute deviation from 1 among the
non‑1 factors is chosen.
"""
# no docstring here
function best_common_factor(coeffs::AbstractVector{<:Number})
    # Remove zeros (if any).
    coeffs_nonzero = filter(c -> !iszero(c), coeffs)
    if isempty(coeffs_nonzero)
        return one(promote_type(eltype(coeffs), Int))
    end
    candidates = unique(coeffs_nonzero)
    best_candidate = candidates[1]
    best_score = (-Inf, Inf)  # (count_ones, candidate_amplitude) 
    # We want to maximize count_ones, then minimize amplitude.
    for cand in candidates
        # Compute the factors obtained by dividing each coefficient by the candidate.
        factors = [c / cand for c in coeffs_nonzero]
        count_ones = count(f -> isone(f), factors)
        candidate_amplitude = real(complex(abs(cand)))
        score = (count_ones, -candidate_amplitude)
        if score > best_score
            best_score = score
            best_candidate = cand
        end
    end
    return best_candidate
end

"""
    find_and_apply_common_factor(qterms::Vector{qTerm}) -> (common, new_terms)

For a vector of qTerm's, first remove any with a zero coefficient, then find the best common
factor among the coefficients using `best_common_factor`. Returns a tuple:
- common: the chosen common factor.
- new_terms: the qTerm's updated with their coefficient divided by the common factor.
"""
# no docstring here
function find_and_apply_common_factor(qterms::Vector{qTerm}; remove_exponents::Bool=false)
    # implicitly assume each qExpr is a qTerm
    # Remove any terms with a zero coefficient.
    nonzero_terms = filter(t -> !iszero(t.coeff), qterms)
    # Collect the coefficients.
    coeffs = [t.coeff for t in nonzero_terms]
    common = best_common_factor(coeffs)
    # Divide each term's coefficient by the common factor.
    if remove_exponents
        zero_exponents = zeros(Int, length(nonzero_terms[1].var_exponents))
        new_terms = [qTerm(t.coeff / common, zero_exponents, t.op_indices) for t in nonzero_terms]
    else
        new_terms = [qTerm(t.coeff / common, t.var_exponents, t.op_indices) for t in nonzero_terms]
    end
    return common, new_terms
end
function separate_into_common_subterms(qexpr::Vector{qExpr}, statespace::StateSpace)::Vector{Tuple{Number,Vector{Int},Vector{qExpr}}}
    nonzero_qobj = filter(t -> !iszero(t), qexpr)
    sorted_q = sort(nonzero_qobj)
    sorted_expr = copy(sorted_q)
    if isempty(sorted_expr)
        return Tuple{Number,Vector{Int},Vector{qExpr}}[]
    elseif length(sorted_expr) == 1
        zero_exponents = zeros(Int, length(statespace.vars))
        common_factor = 1
        return [(common_factor, zero_exponents, sorted_expr)]
    end
    grouped_expressions::Vector{Tuple{Number,Vector{Int},Vector{qExpr}}} = Tuple{Number,Vector{Int},Vector{qExpr}}[]
    curr_qobj_qterm::Vector{qTerm} = []
    # first group the terms 
    curr_i = 1
    i = 1
    curr_var_exponents::Vector{Int} = var_exponents(sorted_expr[1])
    zero_exponents = zeros(Int, length(curr_var_exponents))
    curr_qobj::Vector{qExpr} = [sorted_expr[1]]
    while i < length(sorted_expr)
        i += 1
        if iszero(curr_var_exponents)
            push!(grouped_expressions, (1, zero_exponents, curr_qobj))
            curr_i = i
            curr_var_exponents = var_exponents(sorted_expr[i])
            curr_qobj = [sorted_expr[i]]
        elseif curr_var_exponents == var_exponents(sorted_expr[i])
            push!(curr_qobj, sorted_expr[i])
        else
            curr_qobj_qterm = [x::qTerm for x in curr_qobj]
            if length(curr_qobj_qterm) == 1
                push!(grouped_expressions, (1, zero_exponents, curr_qobj_qterm))
            else
                # find common factor
                common_factor, corrected_terms = find_and_apply_common_factor(curr_qobj_qterm, remove_exponents=true)
                push!(grouped_expressions, (common_factor, curr_var_exponents, corrected_terms))
            end
            curr_i = i
            curr_var_exponents = var_exponents(sorted_expr[i])
            curr_qobj = [sorted_expr[i]]
        end
    end
    if iszero(curr_var_exponents)
        common_factor = 1
        push!(grouped_expressions, (common_factor, zero_exponents, curr_qobj))
    else
        # turn into qTerms
        curr_qobj_qterm = [x::qTerm for x in curr_qobj]
        if length(curr_qobj_qterm) == 1
            push!(grouped_expressions, (1, zero_exponents, curr_qobj_qterm))
        else
            # find common factor
            common_factor, corrected_terms = find_and_apply_common_factor(curr_qobj_qterm, remove_exponents=true)
            push!(grouped_expressions, (common_factor, curr_var_exponents, corrected_terms))
        end
    end
    return grouped_expressions
end
# --- Functions to print an entire qExpr -----------------------------
function qeq_to_string_ungrouped(q::qExpr; do_latex::Bool=false, do_sigma::Bool=false, braket::Bool=true)::String
    if length(q) > 0
        term_string = strip(join([qobj_to_string(t, q.statespace, do_latex=do_latex, do_sigma=do_sigma, braket=braket) for t in q]))
    else
        term_string = "0"
    end
    term_string = lstrip(term_string, '+')
    return term_string
end
function brace(grouped_term::Tuple{Number,Vector{Int},Vector{qExpr}}, statespace::StateSpace; do_latex::Bool=false, do_sigma::Bool=false, braket::Bool=true)
    coeff::Number = grouped_term[1]
    var_exponents::Vector{Int} = grouped_term[2]
    qexprs::Vector{qExpr} = grouped_term[3]
    if length(qexprs) == 0
        return ""
    elseif length(qexprs) == 1
        if !isone(coeff) || !iszero(var_exponents)
            error("Invalid term: $grouped_term.")
        end
        term_string = qobj_to_string(qexprs[1], statespace, do_latex=do_latex, do_sigma=do_sigma, braket=braket)
        prefactorstring, do_either = prefactor_to_string(coeff, var_exponents, statespace, true, do_latex=do_latex)
        if do_either
            if do_latex
                return prefactorstring * "\\left(" * term_string * "\\right)"
            else
                return prefactorstring * "(" * term_string * ")"
            end
        else
            return term_string
        end
    else # need to brace 
        if length(qexprs) > 0
            term_string = strip(join([qobj_to_string(t, statespace, do_latex=do_latex, do_sigma=do_sigma, braket=braket) for t in qexprs]))
            term_string = lstrip(term_string, '+')
            if length(term_string) > 0
                prefactorstring, do_either = prefactor_to_string(coeff, var_exponents, statespace, true, do_latex=do_latex)
                if do_either
                    if do_latex
                        return prefactorstring * "\\left(" * term_string * "\\right)"
                    else
                        return prefactorstring * "(" * term_string * ")"
                    end
                else
                    return term_string
                end
            else
                return ""
            end

        else
            return ""
        end
    end
end
function qeq_to_string_grouped(q::qExpr; do_latex::Bool=false, do_sigma::Bool=false, braket::Bool=true)::String
    grouped_expressions = separate_into_common_subterms(q.terms, q.statespace)
    term_string = strip(join([brace(t, q.statespace, do_latex=do_latex, do_sigma=do_sigma, braket=braket) for t in grouped_expressions]))
    return lstrip(term_string, '+')
end

function qeq_to_string(q::qExpr; do_latex, do_sigma::Bool=false, braket::Bool=true, grouped::Bool=true)::String
    q_simplified = simplify(q)
    if grouped
        return qeq_to_string_grouped(q_simplified; do_latex=do_latex, do_sigma=do_sigma, braket=braket)
    else
        return qeq_to_string_ungrouped(q_simplified; do_latex=do_latex, do_sigma=do_sigma, braket=braket)
    end
end

function show(io::IO, x::qExpr)
    print(io, qeq_to_string(x, do_latex=false, do_sigma=false))
end
function show(io::IO, ::MIME"text/latex", x::qExpr)
    print(io, latexstring(qeq_to_string(x, do_latex=true, do_sigma=false)))
end


function diff_qEQ2string(q::diff_qEQ; do_latex::Bool=false, grouped::Bool=true)::String
    right_hand_side = qeq_to_string(q.right_hand_side, do_latex=do_latex, do_sigma=q.do_sigma, braket=q.braket, grouped=grouped)
    left_hand_side_op_str = qeq_to_string_ungrouped(qExpr([q.left_hand_side], q.statespace), do_latex=do_latex, do_sigma=q.do_sigma, braket=q.braket)
    left_hand_side_op_str = lstrip(left_hand_side_op_str, '+')
    if do_latex
        left_hand_side = raw"\frac{\text{d} \phantom{t}}{\text{d}t}\!" * left_hand_side_op_str
    else
        left_hand_side = "d(" * left_hand_side_op_str * ") / dt"
    end
    return left_hand_side * " = " * right_hand_side
end

function show(io::IO, q::diff_qEQ)
    print(io, diff_qEQ2string(q, do_latex=false))
end
function show(io::IO, ::MIME"text/latex", q::diff_qEQ)
    print(io, latexstring(diff_qEQ2string(q, do_latex=true)))
end

#### String ##########################################################################################################################
""" 
    string(eq::qExpr) -> String
    string(eq::qSum) -> String
    string(eq::diff_qEQ) -> String

Returns a string representation of the qExpr, qSum, or diff_qEQ object. The string is formatted in a human-readable way, but without LaTeX formatting.
"""
function string(eq::qExpr)::String
    return qeq_to_string(eq, do_latex=false, do_sigma=false)
end
function string(eq::qSum)::String
    return qobj_to_string(eq, eq.statespace, do_latex=false, do_sigma=false)
end
function string(eq::diff_qEQ)::String
    return diff_qEQ2string(eq, do_latex=false)
end
#### LaTeXString #####################################################################################################################
""" 
    latex_string(eq::qExpr) -> String
    latex_string(eq::qSum) -> String
    latex_string(eq::diff_qEQ) -> String

Returns a LaTeX string representation of the qExpr, qSum, or diff_qEQ object. 
"""
function latex_string(eq::qExpr)::String
    return qeq_to_string(eq, do_latex=true, do_sigma=false)
end
function latex_string(eq::qSum)::String
    return qobj_to_string(eq, eq.statespace, do_latex=true, do_sigma=false)
end
function latex_string(eq::diff_qEQ)::String
    return diff_qEQ2string(eq, do_latex=true)
end