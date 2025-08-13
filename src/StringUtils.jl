module StringUtils

export subscript_indexes, superscript_indexes, var_substitution, var_substitution_latex, str2sub, str2sup, term_pre_split, separate_terms, expstr_separate, symbol2formatted, brace

"""
    subscript_indexes::Dict{Char, String}

Contains the mapping from characters to their subscript representation for non-latex formatted outputs.
"""
const subscript_indexes = Dict('a' => "ₐ", 'h' => "ₕ", 'i' => "ᵢ", 'j' => "ⱼ", 'k' => "ₖ", 'l' => "ₗ", 'm' => "ₘ", 'n' => "ₙ", 
    'o' => "ₒ", 'p' => "ₚ", '1' => "₁", '2' => "₂", '3' => "₃", '4' => "₄", '5' => "₅", '6' => "₆", '7' => "₇", '8' => "₈", 
    '9' => "₉", '=' => "₌", '+' => "₊", '-' => "₋", '0' => "₀", 'x' => "ₓ", 'y' => "ᵧ", ',' => ",", '∊' => "∊", ' ' => " ", 
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
    '(' => "⁽", ')' => "⁾", '+' => "⁺", '0' => "⁰", 'I' => "ᴵ", 'J' => "ᴶ", 'K' => "ᴷ", 'L' => "ᴸ")
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
    symbol2formatted(symbol::String) -> Tuple{String, String}

Returns a tuple of (`unicode_str`, `latex_str`) for the given `symbol`, using
variable substitution rules. Falls back to the raw `symbol` if no match is found.
"""
function symbol2formatted(symbol::String; do_hat::Bool=false)::Tuple{String, String}
    if haskey(var_substitution, symbol)
        symbol_str = var_substitution[symbol]
    else
        symbol_str =  symbol 
    end
    if haskey(var_substitution_latex, symbol)
        symbol_latex = var_substitution_latex[symbol]
    else
        symbol_latex =  symbol 
    end
    if do_hat
        symbol_latex = raw"\hat{" * symbol_latex * "}"
    end
    return symbol_str, symbol_latex
end


function expstr_separate(expstr::String)::Tuple{String,Int}
    exp::Int = 1
    if occursin("^", expstr)
        expstr, b = split(expstr, "^")
        exp = parse(Int, b)
    end
    return expstr, exp
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

end