module StringUtils

export str2sub, str2sup, indices2str, symbol2formatted, var_unsubstitution, t_suffix
export brace, braket, brace_separate, normalize_underscore_indices, format_normalized_indices, split_index
export int_exponent2str, exponentdag2str, normalize_label, reverse_var_substitution

const SUBSCRIPT_INDEXES = Dict('a' => "ₐ", 'h' => "ₕ", 'i' => "ᵢ", 'j' => "ⱼ", 'k' => "ₖ", 'l' => "ₗ", 'm' => "ₘ", 'n' => "ₙ", 
    'o' => "ₒ", 'p' => "ₚ", 'ρ' => "ᵨ", '1' => "₁", '2' => "₂", '3' => "₃", '4' => "₄", '5' => "₅", '6' => "₆", '7' => "₇", '8' => "₈", 
    '9' => "₉", '=' => "₌", '+' => "₊", '-' => "₋", '0' => "₀", 'x' => "ₓ", 'y' => "ᵧ", ',' => "ˏ", ' ' => " ", 
    '(' => "₍", ')' => "₎", 'e' => "ₑ", 'r' => "ᵣ", 's' => "ₛ", 't' => "ₜ", 'u' => "ᵤ", 'v' => "ᵥ")
const SUPERSCRIPT_INDEXES = Dict('a' => "ᵃ", 'b' => "ᵇ", 'c' => "ᶜ", 'd' => "ᵈ", 'e' => "ᵉ", 'f' => "ᶠ",
    'g' => "ᵍ", 'h' => "ʰ", 'i' => "ⁱ", 'j' => "ʲ", 'k' => "ᵏ", 'l' => "ˡ", 'm' => "ᵐ", 'n' => "ⁿ",
    'o' => "ᵒ", 'p' => "ᵖ", 'q' => "ᵠ", 'r' => "ʳ", 's' => "ˢ", 't' => "ᵗ", 'u' => "ᵘ", 'v' => "ᵛ",
    'w' => "ʷ", 'x' => "ˣ", 'y' => "ʸ", 'z' => "ᶻ", '2' => "²", '3' => "³", '4' => "⁴", '5' => "⁵", 
    '6' => "⁶", '7' => "⁷", '8' => "⁸", '9' => "⁹", '1' => "¹", '-' => "⁻", '=' => "⁼", "." => "·", 
    '(' => "⁽", ')' => "⁾", '+' => "⁺", '0' => "⁰", 'I' => "ᴵ", 'J' => "ᴶ", 'K' => "ᴷ", 'L' => "ᴸ", ',' => "ʼ", '/' => "𝄍", '∊' => "∊")
const VAR_SUBSTITUTION = Dict(
    "alpha"=>"α","beta"=>"β","gamma"=>"γ","delta"=>"δ","epsilon"=>"ε","zeta"=>"ζ","eta"=>"η","theta"=>"θ",
    "iota"=>"ι","kappa"=>"κ","lambda"=>"λ","mu"=>"μ","nu"=>"ν","xi"=>"ξ","omicron"=>"ο","pi"=>"π","rho"=>"ρ",
    "sigma"=>"σ","tau"=>"τ","upsilon"=>"υ","phi"=>"φ","chi"=>"χ","psi"=>"ψ","omega"=>"ω",
    "Alpha"=>"Α","Beta"=>"Β","Gamma"=>"Γ","Delta"=>"Δ","Epsilon"=>"Ε","Zeta"=>"Ζ","Eta"=>"Η","Theta"=>"Θ",
    "Iota"=>"Ι","Kappa"=>"Κ","Lambda"=>"Λ","Mu"=>"Μ","Nu"=>"Ν","Xi"=>"Ξ","Omicron"=>"Ο","Pi"=>"Π","Rho"=>"Ρ",
    "Sigma"=>"Σ","Tau"=>"Τ","Upsilon"=>"Υ","Phi"=>"Φ","Chi"=>"Χ","Psi"=>"Ψ","Omega"=>"Ω")
const VAR_SUBSTITUTION_LATEX = Dict(
    "alpha"=>raw"\alpha","beta"=>raw"\beta","gamma"=>raw"\gamma","delta"=>raw"\delta","epsilon"=>raw"\epsilon",
    "zeta"=>raw"\zeta","eta"=>raw"\eta","theta"=>raw"\theta","iota"=>raw"\iota","kappa"=>raw"\kappa",
    "lambda"=>raw"\lambda","mu"=>raw"\mu","nu"=>raw"\nu","xi"=>raw"\xi","omicron"=>raw"\omicron",
    "pi"=>raw"\pi","rho"=>raw"\rho","sigma"=>raw"\sigma","tau"=>raw"\tau","upsilon"=>raw"\upsilon",
    "phi"=>raw"\phi","chi"=>raw"\chi","psi"=>raw"\psi","omega"=>raw"\omega",
    "Alpha"=>raw"\Alpha","Beta"=>raw"\Beta","Gamma"=>raw"\Gamma","Delta"=>raw"\Delta","Epsilon"=>raw"\Epsilon",
    "Zeta"=>raw"\Zeta","Eta"=>raw"\Eta","Theta"=>raw"\Theta","Iota"=>raw"\Iota","Kappa"=>raw"\Kappa",
    "Lambda"=>raw"\Lambda","Mu"=>raw"\Mu","Nu"=>raw"\Nu","Xi"=>raw"\Xi","Omicron"=>raw"\Omicron",
    "Pi"=>raw"\Pi","Rho"=>raw"\Rho","Sigma"=>raw"\Sigma","Tau"=>raw"\Tau","Upsilon"=>raw"\Upsilon",
    "Phi"=>raw"\Phi","Chi"=>raw"\Chi","Psi"=>raw"\Psi","Omega"=>raw"\Omega",
    "α"=>raw"\alpha","β"=>raw"\beta","γ"=>raw"\gamma","δ"=>raw"\delta","ε"=>raw"\epsilon","ζ"=>raw"\zeta",
    "η"=>raw"\eta","θ"=>raw"\theta","ι"=>raw"\iota","κ"=>raw"\kappa","λ"=>raw"\lambda","μ"=>raw"\mu",
    "ν"=>raw"\nu","ξ"=>raw"\xi","ο"=>raw"\omicron","π"=>raw"\pi","ρ"=>raw"\rho","σ"=>raw"\sigma","τ"=>raw"\tau",
    "υ"=>raw"\upsilon","φ"=>raw"\phi","χ"=>raw"\chi","ψ"=>raw"\psi","ω"=>raw"\omega",
    "Α"=>raw"\Alpha","Β"=>raw"\Beta","Γ"=>raw"\Gamma","Δ"=>raw"\Delta","Ε"=>raw"\Epsilon","Ζ"=>raw"\Zeta",
    "Η"=>raw"\Eta","Θ"=>raw"\Theta","Ι"=>raw"\Iota","Κ"=>raw"\Kappa","Λ"=>raw"\Lambda","Μ"=>raw"\Mu",
    "Ν"=>raw"\Nu","Ξ"=>raw"\Xi","Ο"=>raw"\Omicron","Π"=>raw"\Pi","Ρ"=>raw"\Rho","Σ"=>raw"\Sigma","Τ"=>raw"\Tau",
    "Υ"=>raw"\Upsilon","Φ"=>raw"\Phi","Χ"=>raw"\Chi","Ψ"=>raw"\Psi","Ω"=>raw"\Omega")

const VAR_UNSUBSTITUTION = Dict(v => k for (k, v) in VAR_SUBSTITUTION)
const REVERSE_SUBSCRIPT = Dict(first(v) => string(k) for (k,v) in SUBSCRIPT_INDEXES)
const SUB_CHARS = Set(keys(REVERSE_SUBSCRIPT))  # Set{String} of single-codepoint subscripts

@inline function strip_underscores_and_parens(str::AbstractString)::String
    s = String(str)
    isempty(s) && return s
    s = replace(s, "\\" => "")        # drop LaTeX backslashes
    s = replace(s, r"\s+" => "")      # drop whitespace
    # keep everything up to (but not including) the first '_' or '('
    m = match(r"^[^_(]*", s)
    return m === nothing ? s : m.match
end
strip_underscores_and_parens(sym::Symbol) = strip_underscores_and_parens(String(sym))

function normalize_label(str::AbstractString)::String 
    t = strip_underscores_and_braces(str)
    isempty(t) && return t
    get(VAR_UNSUBSTITUTION, t, t)
end
normalize_label(sym::Symbol)::String = normalize_label(String(sym))

function var_unsubstitution(symbol::AbstractString)::String 
    return get(VAR_UNSUBSTITUTION, symbol, symbol) 
end 
var_unsubstitution(sym::Symbol) = Symbol(var_unsubstitution(String(sym))) 
function reverse_var_substitution(label::AbstractString)::String 
    text = String(label) 
    isempty(text) && return text 
    for (formatted, raw) in VAR_UNSUBSTITUTION 
        occursin(formatted, text) || continue 
        text = replace(text, formatted => raw) 
    end 
    return text 
end 
reverse_var_substitution(sym::Symbol) = reverse_var_substitution(String(sym))

"""
    str2sub(s::String) -> String

Converts the input string `s` into a string with Unicode subscript characters.
For characters not found in `SUBSCRIPT_INDEXES`, falls back to `_c` notation.
"""
function str2sub(s::String)::String
    new_str = ""
    for c in s
        if haskey(SUBSCRIPT_INDEXES, Char(c))
            new_str *= SUBSCRIPT_INDEXES[Char(c)]
        else
            #@warn "Character $c not found in SUBSCRIPT_INDEXES, printing as _$c instead. Avoid this by choosing one of the following characters: $SUBSCRIPT_INDEXES.keys()"
            new_str *= "_$c"
        end
    end
    return new_str
end

"""
    str2sup(s::String) -> String

Converts the input string `s` into a string with Unicode superscript characters.
For characters not found in `SUPERSCRIPT_INDEXES`, falls back to `^c` notation.
"""
function str2sup(s::String)::String
    new_str = ""
    for c in s
        if haskey(SUPERSCRIPT_INDEXES, Char(c))
            new_str *= SUPERSCRIPT_INDEXES[Char(c)]
        else
            #@warn "Character $c not found in SUPERSCRIPT_INDEXES, printing as ^$c instead. Avoid this by choosing one of the following characters: $SUPERSCRIPT_INDEXES.keys()"
            new_str *= "^$c"
        end
    end
    return new_str
end

"""
    symbol2formatted(symbol::String; indices::Vector{String}=String[], do_hat::Bool=false) -> Tuple

Returns a tuple of (`unicode_str`, `latex_str`) for the given `symbol`, using
variable substitution rules. Falls back to the raw `symbol` if no match is found.
Adds a hat on latex output if desired. Alternatively can also create indices. 
"""
function symbol2formatted(symbol::String, indices::Vector{String}=String[]; do_hat::Bool=false)
    # lookup substitutions (default: keep symbol itself)
    symbol_str   = get(VAR_SUBSTITUTION, symbol, symbol)
    symbol_latex = get(VAR_SUBSTITUTION_LATEX, symbol, symbol)

    # apply hat if requested
    if do_hat
        symbol_latex = raw"\hat{" * symbol_latex * "}"
    end

    # handle indices if provided
    if !isempty(indices)
        connector = all([length(i)==1 for i in indices]) ? "," : ""
        index_str_raw = join(indices, connector)
        symbol_str *= str2sub(index_str_raw)
        symbol_latex *= "_{" * index_str_raw * "}"
    end
    return symbol_str, symbol_latex
end

function t_suffix(t_ind::Int; do_latex::Bool=false)
    if t_ind == 0
        return "t"
    else 
        if do_latex
            return "t_$t_ind"
        else
            return "t"*str2sub(string(t_ind))
        end 
    end
end

function indices2str(indices::Vector{Symbol}; do_latex::Bool=false)::String 
    return indices2str(String.(indices), do_latex=do_latex)
end
function indices2str(indices::Vector{Int}; do_latex::Bool=false)::String 
    return indices2str(String.(indices), do_latex=do_latex)
end
function indices2str(indices::Int; do_latex::Bool=false)::String 
    return indices2str([String(indices)], do_latex=do_latex)
end
function indices2str(indices::Symbol; do_latex::Bool=false)::String 
    return indices2str([String.(indices)], do_latex=do_latex)
end
function indices2str(indices::Vector{String}; do_latex::Bool=false)::String 
    if !isempty(indices)
        connector = all([length(i)==1 for i in indices]) ? "," : ""
        index_str_raw = join(indices, connector)
        if do_latex 
            return "_{" * index_str_raw * "}"
        else
            return str2sub(index_str_raw)
        end
    else
        return ""
    end
end

""" 
    brace(x::String; do_latex::Bool=true)::String
Brace a string with parentheses. 
""" 
function brace(x::String; do_latex::Bool=true)::String
    if do_latex 
        return raw"\left(" * x * raw"\right)"
    else
        return "(" * x * ")"
    end
end
""" 
    braket(x::String; do_latex::Bool=true)::String
Brace a string with parentheses. 
""" 
function braket(x::String; do_latex::Bool=true)::String
    if do_latex 
        return raw"\braket{" * x * raw"}"
    else
        return "⟨" * x * "⟩"
    end
end

# Separates strings and the komme separated elements in their braces, so that 
#    "A(B,C)" -> ("A", ["B","C"])
#    "A" -> ("A", [])
function brace_separate(s::String; braces::Tuple{String, String} = ("(",")") )::Tuple{String, Vector{String}}
    s = strip(s)
    if length(braces[1]) != 1 || length(braces[2]) != 1
        error("Only supports single character braces!")
    end
    if occursin(braces[1], s)
        elements::Vector{String} = []
        pref::String = ""
        brace_ind = findall(braces[1], s)
        brace_ind2 = findall(braces[2], s)
        if length(brace_ind) != 1 || length(brace_ind2) != 1
            error("Only supports a single brace pair!")
        end

        brace_ind  = first(brace_ind[1])
        brace_ind2 = first(brace_ind2[1])

        if brace_ind2 != length(s)
            error("Does not support text after closing brace!")
        end
        pref = s[1:brace_ind-1]
        content = s[brace_ind+1:brace_ind2-1]
        # split by comma 
        elements = strip.(split(content, ","))
        return pref, elements 
    else 
        return s, String[]
    end
end

@inline function _split_trailing_args(base::String)
    if endswith(base, raw"\right)")
        m = match(r"^(.*?)(\\left\(.*\\right\))$", base)
        m === nothing || return (m.captures[1], m.captures[2])
    end
    m = match(r"^(.*?)(\([^()]*\))$", base)
    return m === nothing ? (base, "") : (m.captures[1], m.captures[2])
end
function int_exponent2str(base::String, exponent::Int, dag::Bool=false; do_latex::Bool=false)
    exponent == 0 && return ""
    exp_str = ""
    exponent != 1 && (exp_str *= do_latex ? string(exponent) : str2sup(string(exponent)))
    dag && (exp_str *= do_latex ? "*" : "'")
    if do_latex && !isempty(exp_str)
        exp_str = "^{$exp_str}"
    end
    core, suffix = _split_trailing_args(base)
    return core * exp_str * suffix
end
function exp2str(exponent_str::String; do_latex::Bool=false)::String 
    if do_latex && length(exponent_str) > 0
        exponent_str = "^{$exponent_str}"
    end
    return exponent_str
end
function exponentdag2str(base::String, exponent::Union{Int, Rational{Int}}, dag::Bool=false; do_latex::Bool=false)::String
    q = exponent isa Int ? exponent//1 : exponent
    num, den = numerator(q), denominator(q)
    if num == 0 
        return ""
    else
        exponent_str = ""
        if dag 
            exponent_str *= do_latex ? "*" : "'"
        end
        if num == 1 && den == 1   # 1 case
            return base*exp2str(exponent_str, do_latex=do_latex)
        end
        if do_latex 
            if num == 1   # n-root cases
                if den == 2
                    return raw"\sqrt{"*base*exp2str(exponent_str, do_latex=true)*"}"
                else 
                    return raw"\sqrt["*string(den)*"]{"*base*exp2str(exponent_str, do_latex=true)*"}"
                end
            elseif den == 1  # int case
                exponent_str *= " "*string(num)
                return base*exp2str(exponent_str, do_latex=true)
            else  # rational case
                exponent_str *= " "*raw"\tfrac{"*string(num)*"}{"*string(den)*"}"
                return base*exp2str(exponent_str, do_latex=true)
            end
        else 
            if den == 1 # int case 
                exponent_str *= str2sup(string(num))
            else
                exponent_str *= str2sup(string(num)*"/"*string(den))
            end
            return base * exponent_str
        end
    end
end                

include("StringUtilsOps/StringUntils_Indexes.jl")

end
