module StringUtils

export subscript_indexes, superscript_indexes, var_substitution, var_substitution_latex
export str2sub, str2sup, indexes2str, symbol2formatted, t_suffix, brace, braket, brace_separate, underscore_separate
export int_exponent2str, exponentdag2str
"""
    subscript_indexes::Dict{Char, String}

Contains the mapping from characters to their subscript representation for non-latex formatted outputs.
"""
const subscript_indexes = Dict('a' => "ₐ", 'h' => "ₕ", 'i' => "ᵢ", 'j' => "ⱼ", 'k' => "ₖ", 'l' => "ₗ", 'm' => "ₘ", 'n' => "ₙ", 
    'o' => "ₒ", 'p' => "ₚ", '1' => "₁", '2' => "₂", '3' => "₃", '4' => "₄", '5' => "₅", '6' => "₆", '7' => "₇", '8' => "₈", 
    '9' => "₉", '=' => "₌", '+' => "₊", '-' => "₋", '0' => "₀", 'x' => "ₓ", 'y' => "ᵧ", ',' => "ˏ", ' ' => " ", 
    '(' => "₍", ')' => "₎")

"""
    superscript_indexes::Dict{Char, String}

Contains the mapping from characters to their superscript representation for non-latex formatted outputs.
"""
const superscript_indexes = Dict('a' => "ᵃ", 'b' => "ᵇ", 'c' => "ᶜ", 'd' => "ᵈ", 'e' => "ᵉ", 'f' => "ᶠ",
    'g' => "ᵍ", 'h' => "ʰ", 'i' => "ⁱ", 'j' => "ʲ", 'k' => "ᵏ", 'l' => "ˡ", 'm' => "ᵐ", 'n' => "ⁿ",
    'o' => "ᵒ", 'p' => "ᵖ", 'q' => "ᵠ", 'r' => "ʳ", 's' => "ˢ", 't' => "ᵗ", 'u' => "ᵘ", 'v' => "ᵛ",
    'w' => "ʷ", 'x' => "ˣ", 'y' => "ʸ", 'z' => "ᶻ", '2' => "²", '3' => "³", '4' => "⁴", '5' => "⁵", 
    '6' => "⁶", '7' => "⁷", '8' => "⁸", '9' => "⁹", '1' => "¹", '-' => "⁻", '=' => "⁼", "." => "·", 
    '(' => "⁽", ')' => "⁾", '+' => "⁺", '0' => "⁰", 'I' => "ᴵ", 'J' => "ᴶ", 'K' => "ᴷ", 'L' => "ᴸ", ',' => "ʼ", '/' => "𝄍", '∊' => "∊")
const var_substitution = Dict("alpha" => "α", "beta" => "β", "gamma" => "γ", "delta" => "δ", "epsilon" => "ε", "zeta" => "ζ", "eta" => "η", "theta" => "θ", "iota" => "ι", "kappa" => "κ", "lambda" => "λ", "mu" => "μ", "nu" => "ν", "xi" => "ξ", "rho" => "ρ", "sigma" => "σ", "tau" => "τ", "phi" => "φ", "chi" => "χ", "psi" => "ψ", "omega" => "ω", "pi" => "π")
const var_substitution_latex = Dict("alpha" => raw"\alpha", "beta" => raw"\beta", "gamma" => raw"\gamma", "delta" => raw"\delta", "epsilon" => raw"\epsilon", "zeta" => raw"\zeta", "eta" => raw"\eta", "theta" => raw"\theta", "iota" => raw"\iota", "kappa" => raw"\kappa", "lambda" => raw"\lambda", "mu" => raw"\mu", "nu" => raw"\nu", "xi" => raw"\xi", "rho" => raw"\rho", "sigma" => raw"\sigma", "tau" => raw"\tau", "phi" => raw"\phi", "chi" => raw"\chi", "psi" => raw"\psi", "omega" => raw"\omega", "pi" => raw"\pi",
    "α" => raw"\alpha", "β" => raw"\beta", "γ" => raw"\gamma", "δ" => raw"\delta", "ε" => raw"\epsilon", "ζ" => raw"\zeta", "η" => raw"\eta", "θ" => raw"\theta", "ι" => raw"\iota", "κ" => raw"\kappa", "λ" => raw"\lambda", "μ" => raw"\mu", "ν" => raw"\nu", "ξ" => raw"\xi", "ρ" => raw"\rho", "σ" => raw"\sigma", "τ" => raw"\tau", "φ" => raw"\phi", "χ" => raw"\chi", "ψ" => raw"\psi", "ω" => raw"\omega", "π" => raw"\pi")

"""
    str2sub(s::String) -> String

Converts the input string `s` into a string with Unicode subscript characters.
For characters not found in `subscript_indexes`, falls back to `_c` notation.
"""
function str2sub(s::String)::String
    new_str = ""
    for c in s
        if haskey(subscript_indexes, Char(c))
            new_str *= subscript_indexes[Char(c)]
        else
            #@warn "Character $c not found in subscript_indexes, printing as _$c instead. Avoid this by choosing one of the following characters: $subscript_indexes.keys()"
            new_str *= "_$c"
        end
    end
    return new_str
end

"""
    str2sup(s::String) -> String

Converts the input string `s` into a string with Unicode superscript characters.
For characters not found in `superscript_indexes`, falls back to `^c` notation.
"""
function str2sup(s::String)::String
    new_str = ""
    for c in s
        if haskey(superscript_indexes, Char(c))
            new_str *= superscript_indexes[Char(c)]
        else
            #@warn "Character $c not found in superscript_indexes, printing as ^$c instead. Avoid this by choosing one of the following characters: $superscript_indexes.keys()"
            new_str *= "^$c"
        end
    end
    return new_str
end

"""
    symbol2formatted(symbol::String; indexes::Vector{String}=String[], do_hat::Bool=false) -> Tuple

Returns a tuple of (`unicode_str`, `latex_str`) for the given `symbol`, using
variable substitution rules. Falls back to the raw `symbol` if no match is found.
Adds a hat on latex output if desired. Alternatively can also create indexes. 
"""
function symbol2formatted(symbol::String, indexes::Vector{String}=String[]; do_hat::Bool=false)
    # lookup substitutions (default: keep symbol itself)
    symbol_str   = get(var_substitution, symbol, symbol)
    symbol_latex = get(var_substitution_latex, symbol, symbol)

    # apply hat if requested
    if do_hat
        symbol_latex = raw"\hat{" * symbol_latex * "}"
    end

    # handle indexes if provided
    if !isempty(indexes)
        connector = all([length(i)==1 for i in indexes]) ? "," : ""
        index_str_raw = join(indexes, connector)
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

function indexes2str(indexes::Vector{Symbol}; do_latex::Bool=false)::String 
    return indexes2str(String.(indexes), do_latex=do_latex)
end
function indexes2str(indexes::Vector{Int}; do_latex::Bool=false)::String 
    return indexes2str(String.(indexes), do_latex=do_latex)
end
function indexes2str(indexes::Int; do_latex::Bool=false)::String 
    return indexes2str([String(indexes)], do_latex=do_latex)
end
function indexes2str(indexes::Symbol; do_latex::Bool=false)::String 
    return indexes2str([String.(indexes)], do_latex=do_latex)
end
function indexes2str(indexes::Vector{String}; do_latex::Bool=false)::String 
    if !isempty(indexes)
        connector = all([length(i)==1 for i in indexes]) ? "," : ""
        index_str_raw = join(indexes, connector)
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

# Processes strings of the forms:
#   --> "pref_{i,j,k}" and returns ("pref", ["i","j","k"]) 
#   --> "pref_i" and returns ("pref", ["i"]) 
#   --> "pref" and returns ("pref", [])
function underscore_separate(s::String)
    s = string(strip(s)) 
    if occursin("_", s)
        pref::String = ""
        str_split = split(s, "_")
        if length(str_split) != 2
            error("Only supports a single underscore!")
        end
        pref = string(str_split[1])
        sub = string(str_split[2])
        # has {braces}? 
        if occursin("{", str_split[2])
            should_be_empty, indexes = brace_separate(sub, braces=("{","}")) 
            if length(should_be_empty) != 0
                error("If braces {} are used for multi indexing, they must follow the underscore immediately, found $(should_be_empty)!")
            end
            return pref, string.(indexes)
        else
            return pref, [sub]
        end
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
function int_exponents2str(bases::Vector{String}, exponents::Vector{Int}, dag::Bool=false; do_latex::Bool=false)::String
    return join([int_exponent2str(b, x, dag; do_latex=do_latex) for (b,x) in zip(bases, exponents)])
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

end