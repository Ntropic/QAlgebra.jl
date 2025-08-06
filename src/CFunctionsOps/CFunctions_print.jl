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


# --- Generic string constructor ---
"""
    stringer(f::CFunction, vars::Vector{String};
             do_latex::Bool = false,
             do_frac::Bool = true) -> (sign::Bool, body::String)

Internal helper that converts an `CFunction` into a signed string using provided variable names.
Used by `to_string`; returns a sign flag and a formatted string (in LaTeX or plain text).
"""
function stringer(f::CFunction)::Tuple{String,String}
    error("No fallback method for CFunction type: $(typeof(f))")
end

# --- CAtom ---
function stringer(a::CAtom; braced::Bool=false)::Tuple{Bool,String}  # true = minus, false = plus
    if isnumeric(a)
        return sign_string(a.coeff) 
    else
        c = a.coeff
        vec = string(a.var_exponents)
        sign, c_str = sign_string(c) 
        if is_abs_one(c)
            return sign, vec
        else
            return sign,  c_str*"*"*vec
        end
    end
end

# --- CSum ---
function stringer(s::CSum; braced::Bool=false) ::Tuple{Bool,String}
    # braced specifies whether the terms will be grouped, so that an external sign is needed
    terms = s.terms
    if braced 
        if allnegative(s) || (FLIP_IF_FIRST_TERM_NEGATIVE && allnegative(s[1]))
            # if all negative and braced, we can just negate the whole thing
            sig = true
            _, body = stringer(-s, braced=braced)
            return sig, body
        else
            sig = false
            body = join(stringer(t, braced=braced) for t in terms)
        end
    end

    # process each term into (sign, body)
    parts = String[]
    for (i, t) in enumerate(terms)
        sig, body = stringer(t, braced=braced)
        if i == 1
            # first term keeps its sign, but no space if positive
            push!(parts, sig == true ? "-$body" : body)
        else
            push!(parts, sig == true ? "-$body" : "+$body")
        end
    end
    return false, join(parts, "")
end


# --- CRational ---
function stringer(r::CRational; braced::Bool=false)::Tuple{Bool,String}
    n = r.numer
    d = r.denom

    # check first term in numerator
    n_sig, n_body = stringer(n, braced=true)
    d_sig, d_body = stringer(d, braced=braced)  # we ignore sign of denom

    if length(n) > 1
        n_body = "($n_body)"
    end
    if length(d) > 1
        d_body = "($d_body)"
    end

    return n_sig, "$n_body/$d_body"
end

function stringer(p::CProd; braced::Bool=false)::Tuple{Bool,String}
    signs, parts = [], []
    for ps in p.terms 
        sig, part = stringer(ps, braced=true)
        push!(signs, sig)
        push!(parts, part)
    end
    # for every sign = True (i.e. negative) flip the sign
    how_many_negatives = sum(signs)
    return how_many_negatives % 2 == 1 , join(parts, "*")
end

function stringer(e::CExp; braced::Bool=false)::Tuple{Bool,String}
    # always treat exp(x) as positive outside
    _, body = stringer(e.x, braced=braced)
    sign, c_str = sign_string(e.coeff) 
    if is_abs_one(e.coeff)
        c_str = ""
    end
    return sign, c_str*"exp(" * body * ")"
end

function stringer(l::CLog; braced::Bool=false)::Tuple{Bool,String}
    _, body = stringer(l.x, braced=braced)
    sign, c_str = sign_string(l.coeff)
    if is_abs_one(l.coeff)
        c_str = ""
    end
    return sign, c_str*"log(" * body * ")"
end

################# Show ###########################################################################
import Base: show

function show(io::IO, f::CFunction)
    sig, body = stringer(f)
    sig_str = sig ? "-" : ""
    print(io, sig_str * body)
end

################ Stringer with variable names ####################################################
# Generic fallback
function stringer(f::CFunction, vars::Vector{String}; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=true)
    error("No stringer method for type $(typeof(f)) with variable names")
end

function stringer(a::CAtom, vars::Vector{String}; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=true)
    exps = a.var_exponents
    @assert length(vars) == length(exps) "Number of symbols must match number of variables"
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
                    push!(varparts, e == 1 ? vars[i] : "$(vars[i])^{$e}")
                else
                    push!(varparts, e == 1 ? vars[i] : "$(vars[i])"*str2sup(string(e)))
                end
            end
            var_str = isempty(varparts) ? "" : join(varparts, "")
            sign, c_str = sign_string(c) 
            if is_abs_one(c)
                return sign, var_str
            else
                return sign, c_str * connector * var_str
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
                    push!(pos_parts, e == 1 ? vars[i] : "$(vars[i])^{$e}")
                end
                for (i, e) in zip(neg_inds, exps[neg_inds])
                    push!(neg_parts, e == -1 ? vars[i] : "$(vars[i])^{$(-e)}")
                end
            else
                for (i, e) in zip(pos_inds, exps[pos_inds])
                    push!(pos_parts, e == 1 ? vars[i] : "$(vars[i])"*str2sup(string(e)))
                end 
                for (i, e) in zip(neg_inds, exps[neg_inds])
                    push!(neg_parts, e == -1 ? vars[i] : "$(vars[i])"*str2sup(string(-e)))
                end
            end
            if length(neg_inds) == 0 
                sign, c_str = sign_string(c) 
                var_str = isempty(pos_parts) ? "" : join(pos_parts, "")
                if is_abs_one(c)
                    return sign, var_str
                else
                    if do_latex
                        return sign, c_str*" "*var_str
                    else
                        return sign, c_str*var_str
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

function stringer(s::CSum, vars::Vector{String}; do_latex::Bool=false, braced::Bool=false, do_frac::Bool=true) # braced must be false by default for this to work! 
    terms = s.terms
    if isempty(terms)
        return false, "0"
    end
    if braced && (allnegative(s)|| (FLIP_IF_FIRST_TERM_NEGATIVE && allnegative(s[1])))
        sig = true
        _, body = stringer(-s, vars; do_latex=do_latex, do_frac=do_frac)
        return sig, body
    end

    parts = String[]
    for (i, t) in enumerate(terms)
        sig, body = stringer(t, vars; do_latex=do_latex, do_frac=do_frac)
        if i == 1
            push!(parts, sig ? "-" * body : body)
        else
            push!(parts, sig ? "-" * body : "+" * body)
        end
    end

    out = join(parts, "")
    return false, out
end

function stringer(r::CRational, vars::Vector{String}; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=true)
    n = r.numer
    d = r.denom

    n_sig, n_str = stringer(n, vars; do_latex=do_latex, braced=braced, do_frac=false)
    d_sig, d_str = stringer(d, vars; do_latex=do_latex, braced=braced, do_frac=false)
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

function stringer(r::CProd, vars::Vector{String}; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    c = r.coeff
    sign, c_str = sign_string(c) 
    var_strings = [stringer(x, vars; do_latex=do_latex, do_frac=do_frac, braced=braced) for x in r.terms]
    var_str = join(var_strings, "")
    connector = do_latex ? " " : ""
    if is_abs_one(c)
        return sign, var_str
    else
        return sign, c_str * connector * var_str
    end
end

function stringer(e::CLog, vars::Vector{String}; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    c = e.coeff
    sign, c_str = sign_string(c) 
    sign_x, x_str = stringer(e.x, vars; do_latex=do_latex, do_frac=do_frac, braced=false)
    x_str_signed = sign_x ? "-"*x_str : ""*x_str
    connector = do_latex ? " " : ""
    logger = do_latex ? "\\log " : "log"
    if is_abs_one(c)
        return sign, logger*brace(x_str_signed, do_latex=do_latex)
    else
        return sign, c_str * connector * logger*brace(x_str_signed, do_latex=do_latex)
    end
end
function stringer(e::CExp, vars::Vector{String}; do_latex::Bool=false, do_frac::Bool=true, braced::Bool=false)
    c = e.coeff
    sign, c_str = sign_string(c) 
    sign_x, x_str = stringer(e.x, vars; do_latex=do_latex, do_frac=do_frac, braced=braced) 
    x_str_signed = sign_x ? "-"*x_str : ""*x_str
    connector = do_latex ? " " : ""
    logger = do_latex ? "\\exp " : "exp"
    if is_abs_one(c)
        return sign, logger*brace(x_str_signed, do_latex=do_latex)
    else
        return sign, c_str * connector * logger*brace(x_str_signed, do_latex=do_latex)
    end
end


function to_stringer(f::CFunction, vars::Vector{String}; do_latex::Bool=false, braced::Bool=false, do_frac::Bool=true, has_op::Bool=false)::Tuple{Bool, String}
    if braced && isa(f, CSum) && length(f) > 1
        sig, body = stringer(f, vars; do_latex=do_latex, braced=braced)
        # Apply braces if necessary
        body = brace(body, do_latex=do_latex)
        # braced attempt
        done, pre_f, new_f = separate_CSum(f)
        if done 
            sig_new, body_new = stringer(new_f, vars; do_latex=do_latex, braced=braced, do_frac=do_frac)
            body_new = do_latex ? "\\left( $body_new \\right)" : "($body_new)"
            # add prefactor (i.e. base) if it isn'T trivial (isonelike)
            sig_outer, prefactor_outer = stringer(pre_f, vars; do_latex=do_latex, do_frac=do_frac)
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
        sig, body = stringer(f, vars; do_latex=do_latex, do_frac=do_frac)
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
    to_string(f::CFunction, vars::Vector{String};
              do_latex::Bool = false,
              braced::Bool = false,
              optional_sign::Bool = true) -> String

Converts an `CFunction` into a human‐readable string using variable names in `vars`.
- If `do_latex=true`, uses LaTeX syntax (e.g. `\frac{}` and superscripts).
- If `braced=true`, wraps sums in parentheses.
- If `optional_sign=false`, always prefixes a plus or minus sign.
"""
function to_string(f::CFunction, vars::Vector{String}; do_latex::Bool=false, braced::Bool=true, optional_sign::Bool=true, do_frac::Bool=true)::String
    f = simplify(f)
    sig, body = to_stringer(f, vars; do_latex=do_latex, braced=braced, do_frac=do_frac)
    if sig
        return "-" * body
    else
        if optional_sign
            return body
        end
        return "+" *body
    end
end