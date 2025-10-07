using LaTeXStrings
using ..CFunctions
export stringer, to_stringer, to_string

# --- Small, inlined helpers used everywhere ---
@inline connector(do_latex::Bool) = do_latex ? " " : ""

@inline function with_coeff(c::ComplexRational, body::AbstractString; do_latex::Bool=false)
    sig, c_str = sign_string(c, do_latex)
    if is_abs_one(c)
        # drop the visible "1" unless body is empty (pure ±1)
        return sig, isempty(body) ? c_str : body
    else
        return sig, c_str * connector(do_latex) * body
    end
end

# Use this when you've already accumulated signs from child terms
@inline function with_coeff_xor(acc_sig::Bool, c::ComplexRational, body::AbstractString; do_latex::Bool=false)
    sig, out = with_coeff(c, body; do_latex=do_latex)
    return xor(acc_sig, sig), out
end

function sign_string(c::ComplexRational, do_latex::Bool=false)::Tuple{Bool, String}
    neg = is_negative(c)
    body = string(neg ? -c : c, do_latex=do_latex)
    return neg, body
end
is_abs_one(c::ComplexRational)::Bool = (abs(c.a) == abs(c.c))
function is_abs_one(c::CFunction)
    if isnumeric(c)
        if isa(c, Union{CAtom, CAtomIndexed})
            return is_abs_one(c.coeff)
        elseif isa(c, CSum)
            if length(c) == 1
                return is_abs_one(c[1])
            else
                return false
            end
        else
            return false
        end
        return true
    else
        return false
    end
end

# superscript (plain text) or latex brace for exponent
_pow_sup_int(n::Int; do_latex::Bool=false) = do_latex ? "^{$n}" : str2sup(string(n))
_pow_sup_frac(p::Int, q::Int; do_latex::Bool=false) = do_latex ? "^\\{\\frac{$p}{$q}\\}" : "^(" * string(p) * "/" * string(q) * ")"

@inline function get_params(a::CFunction; do_latex::Bool=false)::Vector{String}
    return do_latex ? a.param_info.params_latex : a.param_info.params_str
end


# --- Generic string constructor ---
"""
    stringer(f::CFunction;
             do_latex::Bool = false,
             do_frac::Bool = true) -> (sign::Bool, body::String)

Internal helper that converts an `CFunction` into a signed string using provided variable names.
Used by `to_string`; returns a sign flag and a formatted string (in LaTeX or plain text).
"""
function stringer(f::CFunction; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=true)
    error("No stringer method for type $(typeof(f)) with variable names")
end

function _stringer_atom(coeff::ComplexRational, exps::SparseVector{Int,Int}, params::Vector{String};
                        do_latex::Bool, do_frac::Bool, is_numeric::Bool)
    if is_numeric
        return sign_string(coeff, do_latex)
    end

    if !do_frac
        param_str = join((int_exponent2str(b, x; do_latex=do_latex) for (b, x) in zip(params, exps)), "")
        return with_coeff(coeff, param_str; do_latex=do_latex)
    else
        pos_inds = findall(>(0), exps)
        neg_inds = findall(<(0), exps)

        pos_str = join((int_exponent2str(b, x; do_latex=do_latex) for (b, x) in zip(params[pos_inds], exps[pos_inds])), "")
        neg_str = join((int_exponent2str(b, abs(x); do_latex=do_latex) for (b, x) in zip(params[neg_inds], exps[neg_inds])), "")

        if isempty(neg_inds)
            return with_coeff(coeff, pos_str; do_latex=do_latex)
        else
            c_num = ComplexRational(coeff.a, coeff.b, 1)
            c_den = ComplexRational(coeff.c, 0, 1)

            _, c_pos = sign_string(c_num, do_latex)
            _, c_neg = sign_string(c_den, do_latex)

            if !is_abs_one(c_num) || isempty(pos_inds)
                c_pos *= connector(do_latex)
            else
                c_pos = ""
            end
            if !is_abs_one(c_den)
                c_neg *= connector(do_latex)
            else
                c_neg = ""
            end

            num_str   = c_pos * pos_str
            denom_str = c_neg * neg_str

            sig, _ = sign_string(coeff, do_latex)
            body = do_latex ? raw"\frac{" * num_str * "}{" * denom_str * "}" : num_str * "/(" * denom_str * ")"
            return sig, body
        end
    end
end

function stringer(a::CAtom; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=true)
    params = get_params(a, do_latex=do_latex)
    exps = a.var_exponents
    @assert length(params) == length(exps) "Number of symbols must match number of variables"
    return _stringer_atom(a.coeff, exps, params; do_latex=do_latex, do_frac=do_frac, is_numeric=isnumeric(a))
end

function stringer(a::CAtomIndexed; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=true)
    params = CFunctions._indexed_parameter_names(a, do_latex)
    exps = a.var_exponents
    @assert length(params) == length(exps) "Number of symbols must match number of variables"
    return _stringer_atom(a.coeff, exps, params; do_latex=do_latex, do_frac=do_frac, is_numeric=isnumeric(a))
end

function stringer(C::CAbstract; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    base = do_latex ? C.abstract_def.latex : C.abstract_def.name
    base = exponentdag2str(base, C.exponent, C.dag; do_latex=do_latex)
    return with_coeff(C.coeff, base; do_latex=do_latex)
end

function stringer(I::CIntegral; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    def = integral_definition(I)
    sig_int, body_int = stringer(def.expr; do_latex=do_latex, do_frac=do_frac, braced=true)
    base = do_latex ? raw"w_{\rho}" : "w_ρ"
    inner = (sig_int ? "-" : "") * body_int
    wrapper = base * "(" * inner * ")"
    return with_coeff(I.coeff, wrapper; do_latex=do_latex)
end

function stringer(C::CCustomType; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    base::String = if do_latex
        latex_str = C.ctype_def.latex
        startswith(latex_str, '\\') ? latex_str : (raw"\textrm{" * latex_str * "}")
    else
        C.ctype_def.plain
    end

    if C.ctype_def.has_abstract
        args = String[]
        for x in C.expr
            s, b = stringer(x; do_latex=do_latex, do_frac=do_frac, braced=braced)
            push!(args, (s ? "-" : "") * b)
        end
        if !isempty(args); base *= "(" * join(args, ",") * ")"; end
    else
        indexes, times = where_acting_to_index_strings(C; do_latex=do_latex)   #.ctype_def.fun
        base *= indexes2str(indexes; do_latex=do_latex)
        if !isempty(times); base *= "(" * join(times, ",") * ")"; end
    end

    return with_coeff(C.coeff, base; do_latex=do_latex)
end

function stringer(C::CCustomTypeIndexed; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    base::String = if do_latex
        latex_str = C.ctype_def.latex
        startswith(latex_str, '\\') ? latex_str : (raw"\textrm{" * latex_str * "}")
    else
        C.ctype_def.plain
    end

    if C.ctype_def.has_abstract
        args = String[]
        for x in C.expr
            s, b = stringer(x; do_latex=do_latex, do_frac=do_frac, braced=braced)
            push!(args, (s ? "-" : "") * b)
        end
        if !isempty(args)
            base *= "(" * join(args, ",") * ")"
        end
    else
        indexes, times = where_acting_to_index_strings(C; do_latex=do_latex)
        base *= indexes2str(indexes; do_latex=do_latex)
        if !isempty(times)
            base *= "(" * join(times, ",") * ")"
        end
    end

    suffix = CFunctions._indexes_suffix(C.indexes, do_latex)
    base *= suffix

    return with_coeff(C.coeff, base; do_latex=do_latex)
end

function stringer(s::CSum; do_latex::Bool=false, braced::Bool=false, do_frac::Bool=true)
    exprs = s.expr
    isempty(exprs) && return false, "0"

    if braced && (allnegative(s) || (FLIP_IF_FIRST_TERM_NEGATIVE && allnegative(s[1])))
        _, body = stringer(-s; do_latex=do_latex, do_frac=do_frac)
        return true, body
    end

    parts = String[]
    for (i, t) in enumerate(exprs)
        sig, body = stringer(t; do_latex=do_latex, do_frac=do_frac)
        push!(parts, (i == 1) ? (sig ? "-" * body : body) : (sig ? "-" * body : "+" * body))
    end
    return false, join(parts, "")
end

function stringer(r::CRational; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=true)
    n_sig, n_str = stringer(r.numer; do_latex=do_latex, braced=braced, do_frac=false)
    d_sig, d_str = stringer(r.denom; do_latex=do_latex, braced=braced, do_frac=false)
    n_sig = xor(n_sig, d_sig)  # warn if denom negative?
    if do_latex
        return n_sig, "\\frac{ $n_str }{ $d_str }"
    else
        n_wrapped = length(r.numer) > 1 ? "($n_str)" : n_str
        d_wrapped = length(r.denom) > 1 ? "($d_str)" : d_str
        return n_sig, "$n_wrapped/$d_wrapped"
    end
end

function stringer(r::CProd; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    sig_acc = false
    parts = String[]
    for x in r.expr
        s, b = to_stringer(x; do_latex=do_latex, do_frac=do_frac, braced=braced, has_op=true)
        if !isempty(b)
            push!(parts, b)
        end
        sig_acc = xor(sig_acc, s)
    end
    body = join(parts, "")
    return with_coeff_xor(sig_acc, r.coeff, body; do_latex=do_latex)
end

function stringer(e::CLog; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    sx, x_str = stringer(e.expr; do_latex=do_latex, do_frac=do_frac, braced=false)
    x_str_signed = sx ? "-" * x_str : x_str
    logger = do_latex ? "\\log " : "log"
    body = logger * brace(x_str_signed; do_latex=do_latex)
    return with_coeff(e.coeff, body; do_latex=do_latex)
end

function stringer(e::CExp; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    sx, x_str = stringer(e.expr; do_latex=do_latex, do_frac=do_frac, braced=false)
    x_str_signed = sx ? "-" * x_str : x_str

    use_fn = contains_vec_or_mat(e.expr)
    body = if use_fn
        (do_latex ? "\\exp " : "exp") * brace(x_str_signed; do_latex=do_latex)
    else
        do_latex ? " e^{" * x_str_signed * "}" : " exp(" * x_str_signed * ")"
    end

    return with_coeff(e.coeff, body; do_latex=do_latex)
end

function stringer(p::CPower; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    sx, bx = stringer(p.expr; do_latex=do_latex, do_frac=do_frac, braced=true)
    base = (sx ? "-" : "") * bx
    k = numerator(p.exponent)
    m = denominator(p.exponent)

    body = if m == 1 && k > 1
        brace(base; do_latex=do_latex) * _pow_sup_int(k; do_latex=do_latex)
    elseif k == 1 && m > 1 && do_latex
        (m == 2) ? raw"\sqrt{" * base * "}" : raw"\sqrt[" * string(m) * "]{" * base * "}"
    else
        pow = (m == 1) ? _pow_sup_int(k; do_latex=do_latex) : _pow_sup_frac(k, m; do_latex=do_latex)
        brace(base; do_latex=do_latex) * pow
    end

    return with_coeff(p.coeff, body; do_latex=do_latex)
end

function stringer(v::CVector; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    elems = Vector{String}(undef, length(v.expr))
    for i in eachindex(v.expr)
        si, bi = stringer(v.expr[i]; do_latex=do_latex, do_frac=do_frac, braced=false)
        elems[i] = (si ? "-" : "") * bi
    end

    body = if do_latex
        v.row ? (raw"\begin{bmatrix} " * join(elems, " & ") * raw" \end{bmatrix}") :
                (raw"\begin{bmatrix} " * join(elems, raw" \\ ") * raw" \end{bmatrix}")
    else
        v.row ? "[ " * join(elems, ", ") * " ]" :
                "[ " * join(elems, " ; ") * " ]"
    end

    return with_coeff(v.coeff, body; do_latex=do_latex)
end

function stringer(M::CMatrix; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    m, n = size(M.expr)
    rows = String[]
    if do_latex
        for i in 1:m
            cols = String[]
            for j in 1:n
                sij, bij = stringer(M.expr[i,j]; do_latex=do_latex, do_frac=do_frac, braced=false)
                push!(cols, (sij ? "-" : "") * bij)
            end
            push!(rows, join(cols, " & "))
        end
        body = raw"\begin{bmatrix} " * join(rows, raw" \\ ") * raw" \end{bmatrix}"
    else
        for i in 1:m
            cols = String[]
            for j in 1:n
                sij, bij = stringer(M.expr[i,j]; do_latex=do_latex, do_frac=do_frac, braced=false)
                push!(cols, (sij ? "-" : "") * bij)
            end
            push!(rows, "[ " * join(cols, ", ") * " ]")
        end
        body = "[" * join(rows, "; ") * "]"
    end

    return with_coeff(M.coeff, body; do_latex=do_latex)
end


function to_stringer(f::CFunction; do_latex::Bool=false, braced::Bool=false, do_frac::Bool=true, has_op::Bool=false)::Tuple{Bool, String}
    if braced && isa(f, CSum) && length(f) > 1
        sig, body = stringer(f; do_latex=do_latex, braced=braced)
        body = brace(body, do_latex=do_latex)
        done, pre_f, new_f = separate_CSum(f)
        if done
            sig_new, body_new = stringer(new_f; do_latex=do_latex, braced=braced, do_frac=do_frac)
            body_new = do_latex ? "\\left( $body_new \\right)" : "($body_new)"
            sig_outer, prefactor_outer = stringer(pre_f; do_latex=do_latex, do_frac=do_frac)
            sig_new = xor(sig_new, sig_outer)
            if !isonelike(pre_f)
                body_new = prefactor_outer * (do_latex ? " " : "") * body_new
            end
            if length(body_new) < length(body)
                return sig_new, body_new
            end
        end
        return sig, body
    else
        sig, body = stringer(f; do_latex=do_latex, do_frac=do_frac)
        if has_op && isnumeric(f)
            # Drop pure unit scalars in a product/sum context (incl. CRational == 1)
            if isonelike(f)
                return sig, ""
            end
        end
        return sig, body
    end
end



"""
    to_string(f::CFunction;
              do_latex::Bool = false,
              braced::Bool = false,
              optional_sign::Bool = true) -> String

Converts an `CFunction` into a human‐readable string using variable names in `params`.
- If `do_latex=true`, uses LaTeX syntax (e.g. `\frac{}` and superscripts).
- If `braced=true`, wraps sums in parentheses.
- If `optional_sign=false`, always prefixes a plus or minus sign.
"""
function to_string(f::CFunction; do_latex::Bool=false, braced::Bool=true, optional_sign::Bool=true, do_frac::Bool=true)::String
    #f = simplify(f)
    sig, body = to_stringer(f; do_latex=do_latex, braced=braced, do_frac=do_frac)
    if sig
        return "-" * body
    else
        if optional_sign
            return body
        end
        return "+" *body
    end
end

import Base: show

function show(io::IO, f::CFunction)
    print(io, latexstring(to_string(f, do_latex=false)))
end

function show(io::IO, ::MIME"text/latex", f::CFunction)
    print(io, latexstring(to_string(f, do_latex=true)))
end
