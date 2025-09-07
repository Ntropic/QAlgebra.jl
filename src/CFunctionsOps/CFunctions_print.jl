export stringer, to_stringer, to_string

function sign_string(c::ComplexRational, do_latex::Bool=false)::Tuple{Bool, String}
    if is_negative(c)
        return (true, string(c, do_latex=do_latex)[2:end])
    else
        return (false, string(c, do_latex=do_latex))
    end
end
is_abs_one(c::ComplexRational)::Bool = (abs(c.a) == abs(c.c))
function is_abs_one(c::CFunction)
    if isnumeric(c)
        if isa(c, CAtom)
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
    do_latex && return a.param_info.params_str 
    return a.param_info.params_latex
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

function stringer(a::CAtom; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=true)
    params = get_params(a, do_latex=do_latex)
    exps = a.var_exponents
    @assert length(params) == length(exps) "Number of symbols must match number of variables"
    connector = do_latex ? " " : ""
    if isnumeric(a)
        return sign_string(a.coeff) 
    else
        # build the variable part
        c = a.coeff
        if !do_frac
            varparts = String[]
            for (i, e) in enumerate(exps)
                if e == 0
                    continue
                elseif do_latex
                    push!(varparts, e == 1 ? params[i] : "$(params[i])^{$e}")
                else
                    push!(varparts, e == 1 ? params[i] : "$(params[i])"*str2sup(string(e)))
                end
            end
            param_str = isempty(varparts) ? "" : join(varparts, "")
            sign, c_str = sign_string(c) 
            if is_abs_one(c)
                return sign, param_str
            else
                return sign, c_str * connector * param_str
            end
        else
            # group terms with positive and negative exponents 
            pos_inds::Vector{Int} = Int[]
            neg_inds::Vector{Int} = Int[]
            for (i, e) in enumerate(exps)
                if e > 0
                    push!(pos_inds, i)
                elseif e < 0
                    push!(neg_inds, i)
                end
            end
            pos_parts = String[]
            neg_parts = String[]
            if do_latex
                for (i, e) in zip(pos_inds, exps[pos_inds])
                    push!(pos_parts, e == 1 ? params[i] : "$(params[i])^{$e}")
                end
                for (i, e) in zip(neg_inds, exps[neg_inds])
                    push!(neg_parts, e == -1 ? params[i] : "$(params[i])^{$(-e)}")
                end
            else
                for (i, e) in zip(pos_inds, exps[pos_inds])
                    push!(pos_parts, e == 1 ? params[i] : "$(params[i])"*str2sup(string(e)))
                end 
                for (i, e) in zip(neg_inds, exps[neg_inds])
                    push!(neg_parts, e == -1 ? params[i] : "$(params[i])"*str2sup(string(-e)))
                end
            end
            if length(neg_inds) == 0 
                sign, c_str = sign_string(c) 
                param_str = isempty(pos_parts) ? "" : join(pos_parts, "")
                if is_abs_one(c)
                    return sign, param_str
                else
                    if do_latex
                        return sign, c_str*" "*param_str
                    else
                        return sign, c_str*param_str
                    end
                end
            else 
                c_num = ComplexRational(c.a, c.b, 1)
                c_denom = ComplexRational(c.c, 0, 1)
                sign, c_pos = sign_string(c_num)
                c_pos *= connector 
                if is_abs_one(c_num) && !isempty(pos_parts)
                    c_pos = ""
                end
                _, c_neg = sign_string(c_denom) 
                c_neg *= connector 
    
                if is_abs_one(c_denom)
                    c_neg = ""
                end
                pos_str = isempty(pos_parts) ? "" : join(pos_parts, "")
                neg_str = isempty(neg_parts) ? "" : join(neg_parts, "")
                num_str = c_pos * pos_str
                denom_str = c_neg * neg_str
                if do_latex 
                    return sign, raw"\frac{"*num_str*" }{"*denom_str*"}"
                else
                    return sign, num_str*"/("*denom_str*")"
                end
            end
        end
    end
end

function stringer(s::CSum; do_latex::Bool=false, braced::Bool=false, do_frac::Bool=true) # braced must be false by default for this to work! 
    params = get_params(s, do_latex=do_latex)
    terms = s.terms
    if isempty(terms)
        return false, "0"
    end
    if braced && (allnegative(s)|| (FLIP_IF_FIRST_TERM_NEGATIVE && allnegative(s[1])))
        sig = true
        _, body = stringer(-s, params; do_latex=do_latex, do_frac=do_frac)
        return sig, body
    end

    parts = String[]
    for (i, t) in enumerate(terms)
        sig, body = stringer(t, params; do_latex=do_latex, do_frac=do_frac)
        if i == 1
            push!(parts, sig ? "-" * body : body)
        else
            push!(parts, sig ? "-" * body : "+" * body)
        end
    end

    out = join(parts, "")
    return false, out
end

function stringer(r::CRational, params::Vector{String}; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=true)
    n = r.numer
    d = r.denom

    n_sig, n_str = stringer(n, params; do_latex=do_latex, braced=braced, do_frac=false)
    d_sig, d_str = stringer(d, params; do_latex=do_latex, braced=braced, do_frac=false)
    if d_sig
        @warn "Denominator is negative in rational expression. Not meant to be"
        n_sig = xor(n_sig, d_sig)
    end
    if do_latex
        return n_sig, "\\frac{ $n_str }{ $d_str }"
    else
        n_wrapped = length(n) > 1 ? "($n_str)" : n_str 
        d_wrapped = length(d) > 1 ? "($d_str)" : d_str
        return n_sig, "$n_wrapped/$d_wrapped"
    end
end

function stringer(r::CProd; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    params = get_params(r, do_latex=do_latex)
    c = r.coeff
    sign, c_str = sign_string(c) 
    var_strings = [stringer(x, params; do_latex=do_latex, do_frac=do_frac, braced=braced) for x in r.terms]
    param_str = join(var_strings, "")
    connector = do_latex ? " " : ""
    if is_abs_one(c)
        return sign, param_str
    else
        return sign, c_str * connector * param_str
    end
end

function stringer(e::CLog; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    params = get_params(e, do_latex=do_latex)
    c = e.coeff
    sign, c_str = sign_string(c) 
    sign_x, x_str = stringer(e.x, params; do_latex=do_latex, do_frac=do_frac, braced=false)
    x_str_signed = sign_x ? "-"*x_str : ""*x_str
    connector = do_latex ? " " : ""
    logger = do_latex ? "\\log " : "log"
    if is_abs_one(c)
        return sign, logger*brace(x_str_signed, do_latex=do_latex)
    else
        return sign, c_str * connector * logger*brace(x_str_signed, do_latex=do_latex)
    end
end
# --- CExp (params) with conditional e^ vs \exp/exp ---
function stringer(e::CExp; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    params = get_params(e, do_latex=do_latex)
    c = e.coeff
    sign, c_str = sign_string(c, do_latex)

    sx, x_str = stringer(e.x, params; do_latex=do_latex, do_frac=do_frac, braced=false)
    x_str_signed = sx ? "-" * x_str : x_str

    use_fn = contains_vec_or_mat(e.x)

    connector = do_latex ? " " : ""

    body =
        if use_fn
            logger = do_latex ? "\\exp " : "exp"
            arg = brace(x_str_signed, do_latex=do_latex)
            logger * arg
        else
            if do_latex
                "e^{" * x_str_signed * "}"
            else
                "e^(" * x_str_signed * ")"
            end
        end

    if is_abs_one(e.coeff)
        return sign, body
    else
        return sign, c_str * connector * body
    end
end

# --- CPower (params) ---
function stringer(p::CPower; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    params = get_params(p, do_latex=do_latex)
    # external sign from coeff
    sig, c_str = sign_string(p.coeff, do_latex)
    sx, bx = stringer(p.x, params; do_latex=do_latex, do_frac=do_frac, braced=true)
    base = (sx ? "-" : "") * bx
    k = numerator(p.exponent)
    m = denominator(p.exponent)
    connector = do_latex ? " " : ""

    body::String = ""
    if m == 1 && k > 1
        # positive integer power => braced base with ^k
        base_b = brace(base, do_latex=do_latex)
        pow = _pow_sup_int(k; do_latex=do_latex)
        body = base_b * pow
    elseif k == 1 && m > 1 && do_latex
        # 1/n exponent -> n-th root (LaTeX)
        body = (m == 2) ? raw"\sqrt{" * base * "}" : raw"\sqrt[" * string(m) * "]{" * base * "}"
    else
        # general rational or non-LaTeX root
        base_b = brace(base, do_latex=do_latex)
        if m == 1
            pow = _pow_sup_int(k; do_latex=do_latex)
        else
            pow = _pow_sup_frac(k, m; do_latex=do_latex)
        end
        body = base_b * pow
    end

    if is_abs_one(p.coeff)
        return sig, body
    else
        return sig, c_str * connector * body
    end
end

# --- CVector (params) ---
function stringer(v::CVector; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    params = get_params(v, do_latex=do_latex)
    sig, c_str = sign_string(v.coeff, do_latex)
    elems = Vector{String}(undef, length(v.expr))
    for i in eachindex(v.expr)
        si, bi = stringer(v.expr[i], params; do_latex=do_latex, do_frac=do_frac, braced=false)
        elems[i] = (si ? "-" : "") * bi
    end

    body::String = ""
    if do_latex
        if v.row
            body = raw"\begin{bmatrix} " * join(elems, " & ") * raw" \end{bmatrix}"
        else
            body = raw"\begin{bmatrix} " * join(elems, raw" \\ ") * raw" \end{bmatrix}"
        end
    else
        body = v.row ? "[ " * join(elems, ", ") * " ]" : "[ " * join(elems, " ; ") * " ]"
    end

    connector = do_latex ? " " : " "
    if is_abs_one(v.coeff)
        return sig, body
    else
        return sig, c_str * connector * body
    end
end
# --- CMatrix (params) ---
function stringer(M::CMatrix; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    params = get_params(M, do_latex=do_latex)
    sig, c_str = sign_string(M.coeff, do_latex)
    m, n = size(M.expr)

    body::String = ""
    if do_latex
        rows = String[]
        for i in 1:m
            cols = String[]
            for j in 1:n
                sij, bij = stringer(M.expr[i,j], params; do_latex=do_latex, do_frac=do_frac, braced=false)
                push!(cols, (sij ? "-" : "") * bij)
            end
            push!(rows, join(cols, " & "))
        end
        body = raw"\begin{bmatrix} " * join(rows, raw" \\ ") * raw" \end{bmatrix}"
    else
        rows = String[]
        for i in 1:m
            cols = String[]
            for j in 1:n
                sij, bij = stringer(M.expr[i,j], params; do_latex=do_latex, do_frac=do_frac, braced=false)
                push!(cols, (sij ? "-" : "") * bij)
            end
            push!(rows, "[ " * join(cols, ", ") * " ]")
        end
        body = "[" * join(rows, "; ") * "]"
    end

    connector = do_latex ? " " : " "
    if is_abs_one(M.coeff)
        return sig, body
    else
        return sig, c_str * connector * body
    end
end

function to_stringer(f::CFunction; do_latex::Bool=false, braced::Bool=false, do_frac::Bool=true, has_op::Bool=false)::Tuple{Bool, String}
    params = get_params(f, do_latex=do_latex)
    if braced && isa(f, CSum) && length(f) > 1
        sig, body = stringer(f, params; do_latex=do_latex, braced=braced)
        # Apply braces if necessary
        body = brace(body, do_latex=do_latex)
        # braced attempt
        done, pre_f, new_f = separate_CSum(f)
        if done 
            sig_new, body_new = stringer(new_f, params; do_latex=do_latex, braced=braced, do_frac=do_frac)
            body_new = do_latex ? "\\left( $body_new \\right)" : "($body_new)"
            # add prefactor (i.e. base) if it isn'T trivial (isonelike)
            sig_outer, prefactor_outer = stringer(pre_f, params; do_latex=do_latex, do_frac=do_frac)
            sig_new = sig_new || sig_outer
            connector = do_latex ? " " : ""
            if !isonelike(pre_f)
                body_new = prefactor_outer*connector*body_new   # add an outer prefactor 
            end 
            if length(body_new) < length(body)
                body = body_new
                sig = sig_new
            end
        end
    else
        sig, body = stringer(f, params; do_latex=do_latex, do_frac=do_frac)
        if has_op && isnumeric(f) && !isa(f, CRational) 
            if isa(f, CAtom)
                f_pre =  f.coeff
            elseif isa(f, CSum)
                f_pre = f[1].coeff
            else 
                error("Unsupported type for CFunction")
            end
            if isonelike(f_pre)  
                return sig, "" 
            end
        end
    end
    return sig, body
end


"""
    to_string(f::CFunction, params::Vector{String};
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
    print(io, latexstring(to_string(f; do_latex=false)))
end

function show(io::IO, ::MIME"text/latex", f::CFunction)
    print(io, latexstring(to_string(f; do_latex=true)))
end