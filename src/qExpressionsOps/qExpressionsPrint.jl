using Printf
using LaTeXStrings
import Base: string

export string, latex_string


function var_exponents2string(var_exponents::Vector{Int}, qspace::StateSpace; do_latex::Bool=false)::String
    monomial::String = ""
    for (i, exp) in enumerate(var_exponents)
        if exp != 0
            # Convert the variable to a string; you may adjust formatting if needed.
            if do_latex
                var = qspace.vars[i].var_latex
            else
                var = qspace.vars[i].var_str
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
function qexpr_to_string(t::qTerm, statespace::StateSpace; do_latex::Bool=false, do_sigma::Bool=false, braket::Bool=false)::String
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
function qexpr_to_string(s::qEQ, statespace::StateSpace; do_latex::Bool=false, do_sigma::Bool=false, braket::Bool=true)::String
    if length(s) > 0
        term_string = strip(join([qexpr_to_string(t, statespace, do_latex=do_latex, do_sigma=do_sigma, braket=braket) for t in s]))
    else
        term_string = "0"
    end
    term_string = lstrip(term_string, '+')
    return term_string
end

function show(io::IO, x::qEQ)
    print(io, qeq_to_string(x, do_latex=false, do_sigma=false))
end
function show(io::IO, ::MIME"text/latex", x::qEQ)
    print(io, latexstring(qeq_to_string(x, do_latex=true, do_sigma=false)))
end

function diff_qEQ2string(q::diff_qEQ; do_latex::Bool=false)::String
    right_hand_side = qeq_to_string(q.right_hand_side, do_latex=do_latex)
    left_hand_side_op_str = qeq_to_string_ungrouped(qEQ([q.left_hand_side], q.statespace), do_latex=do_latex)
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
