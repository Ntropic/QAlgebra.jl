module StringUtils

export subscript_indexes, superscript_indexes, var_substitution, var_substitution_latex
export str2sub, str2sup, term_pre_split, separate_terms, symbol2formatted, t_suffix, brace, match_indexed_pattern, brace_separate, underscore_separate

"""
    subscript_indexes::Dict{Char, String}

Contains the mapping from characters to their subscript representation for non-latex formatted outputs.
"""
const subscript_indexes = Dict('a' => "ₐ", 'h' => "ₕ", 'i' => "ᵢ", 'j' => "ⱼ", 'k' => "ₖ", 'l' => "ₗ", 'm' => "ₘ", 'n' => "ₙ", 
    'o' => "ₒ", 'p' => "ₚ", '1' => "₁", '2' => "₂", '3' => "₃", '4' => "₄", '5' => "₅", '6' => "₆", '7' => "₇", '8' => "₈", 
    '9' => "₉", '=' => "₌", '+' => "₊", '-' => "₋", '0' => "₀", 'x' => "ₓ", 'y' => "ᵧ", ',' => "ˏ", '∊' => "∊", ' ' => " ", 
    '(' => "₍", ')' => "₎")

"""
    superscript_indexes::Dict{Char, String}

Contains the mapping from characters to their superscript representation for non-latex formatted outputs.
"""
const superscript_indexes = Dict('a' => "ᵃ", 'b' => "ᵇ", 'c' => "ᶜ", 'd' => "ᵈ", 'e' => "ᵉ", 'f' => "ᶠ",
    'g' => "ᵍ", 'h' => "ʰ", 'i' => "ⁱ", 'j' => "ʲ", 'k' => "ᵏ", 'l' => "ˡ", 'm' => "ᵐ", 'n' => "ⁿ",
    'o' => "ᵒ", 'p' => "ᵖ", 'q' => "ᵠ", 'r' => "ʳ", 's' => "ˢ", 't' => "ᵗ", 'u' => "ᵘ", 'v' => "ᵛ",
    'w' => "ʷ", 'x' => "ˣ", 'y' => "ʸ", 'z' => "ᶻ", '2' => "²", '3' => "³", '4' => "⁴", '5' => "⁵", 
    '6' => "⁶", '7' => "⁷", '8' => "⁸", '9' => "⁹", '1' => "", '-' => "⁻", '=' => "⁼", "." => "·", 
    '(' => "⁽", ')' => "⁾", '+' => "⁺", '0' => "⁰", 'I' => "ᴵ", 'J' => "ᴶ", 'K' => "ᴷ", 'L' => "ᴸ", ',' => "ʼ")
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
    s = strip(s) 
    if occursin("_", s)
        pref::String = ""
        str_split = split(s, "_")
        if length(str_split) != 2
            error("Only supports a single underscore!")
        end
        pref = String(str_split[1])
        sub = String(str_split[2])
        # has {braces}? 
        if occursin("{", str_split[2])
            should_be_empty, indexes = brace_separate(sub, braces=("{","}")) 
            if length(should_be_empty) != 0
                error("If braces {} are used for multi indexing, they must follow the underscore immediately, found $(should_be_empty)!")
            end
            return pref, String.(indexes)
        else
            return pref, [sub]
        end
    else
        return s, String[]
    end
end
end
