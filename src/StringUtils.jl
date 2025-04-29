module StringUtils

export subscript_indexes, superscript_indexes, var_substitution, var_substitution_latex, str2sub, str2sup, separate_terms, expstr_separate

const subscript_indexes = Dict('a' => "ₐ", 'h' => "ₕ", 'i' => "ᵢ", 'j' => "ⱼ", 'k' => "ₖ", 'l' => "ₗ", 'm' => "ₘ", 'n' => "ₙ", 'o' => "ₒ", 'p' => "ₚ", '1' => "₁", '2' => "₂", '3' => "₃", '4' => "₄", '5' => "₅", '6' => "₆", '7' => "₇", '8' => "₈", '9' => "₉", '=' => "₌", '+' => "₊", '-' => "₋", '0' => "₀", 'x' => "ₓ", 'y' => "ᵧ", ',' => ",", '∊' => "∊", ' ' => " ", '(' => "₍", ')' => "₎")
const superscript_indexes = Dict('a' => "ᵃ", 'b' => "ᵇ", 'c' => "ᶜ", 'd' => "ᵈ", 'e' => "ᵉ", 'f' => "ᶠ",
    'g' => "ᵍ", 'h' => "ʰ", 'i' => "ⁱ", 'j' => "ʲ", 'k' => "ᵏ", 'l' => "ˡ", 'm' => "ᵐ", 'n' => "ⁿ",
    'o' => "ᵒ", 'p' => "ᵖ", 'q' => "ᵠ", 'r' => "ʳ", 's' => "ˢ", 't' => "ᵗ", 'u' => "ᵘ", 'v' => "ᵛ",
    'w' => "ʷ", 'x' => "ˣ", 'y' => "ʸ", 'z' => "ᶻ", '2' => "²", '3' => "³", '4' => "⁴", '5' => "⁵", '6' => "⁶", '7' => "⁷", '8' => "⁸", '9' => "⁹", '1' => "", '-' => "⁻", '=' => "⁼")
const var_substitution = Dict("alpha" => "α", "beta" => "β", "gamma" => "γ", "delta" => "δ", "epsilon" => "ε", "zeta" => "ζ", "eta" => "η", "theta" => "θ", "iota" => "ι", "kappa" => "κ", "lambda" => "λ", "mu" => "μ", "nu" => "ν", "xi" => "ξ", "rho" => "ρ", "sigma" => "σ", "tau" => "τ", "phi" => "φ", "chi" => "χ", "psi" => "ψ", "omega" => "ω", "pi" => "π")
const var_substitution_latex = Dict("alpha" => raw"\alpha", "beta" => raw"\beta", "gamma" => raw"\gamma", "delta" => raw"\delta", "epsilon" => raw"\epsilon", "zeta" => raw"\zeta", "eta" => raw"\eta", "theta" => raw"\theta", "iota" => raw"\iota", "kappa" => raw"\kappa", "lambda" => raw"\lambda", "mu" => raw"\mu", "nu" => raw"\nu", "xi" => raw"\xi", "rho" => raw"\rho", "sigma" => raw"\sigma", "tau" => raw"\tau", "phi" => raw"\phi", "chi" => raw"\chi", "psi" => raw"\psi", "omega" => raw"\omega", "pi" => raw"\pi",
    "α" => raw"\alpha", "β" => raw"\beta", "γ" => raw"\gamma", "δ" => raw"\delta", "ε" => raw"\epsilon", "ζ" => raw"\zeta", "η" => raw"\eta", "θ" => raw"\theta", "ι" => raw"\iota", "κ" => raw"\kappa", "λ" => raw"\lambda", "μ" => raw"\mu", "ν" => raw"\nu", "ξ" => raw"\xi", "ρ" => raw"\rho", "σ" => raw"\sigma", "τ" => raw"\tau", "φ" => raw"\phi", "χ" => raw"\chi", "ψ" => raw"\psi", "ω" => raw"\omega", "π" => raw"\pi")

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

function separate_terms(input::String, total_tokens::Vector{String}, tokens_at_end::Vector{String}, tokens_at_begin::Vector{String})::Tuple{Vector{Vector{String}},Vector{Vector{String}},Vector{Vector{String}}}
    # Prepare output vectors:
    output_total = [String[] for _ in 1:length(total_tokens)]
    output_at_end = [String[] for _ in 1:length(tokens_at_end)]
    output_at_begin = [String[] for _ in 1:length(tokens_at_begin)]

    # Sort each list by descending length (to prioritize longer tokens when there is overlap).
    tokens_total_with_idx = collect(enumerate(total_tokens))
    sort!(tokens_total_with_idx, by=x -> -length(x[2]))

    tokens_at_end_with_idx = collect(enumerate(tokens_at_end))
    sort!(tokens_at_end_with_idx, by=x -> -length(x[2]))

    tokens_at_begin_with_idx = collect(enumerate(tokens_at_begin))
    sort!(tokens_at_begin_with_idx, by=x -> -length(x[2]))

    # Preprocess the input: remove underscores and replace spaces with asterisks.
    input = replace(replace(input, "_" => ""), " " => "*")
    # Split the input string by '*' signs.
    terms = split(input, "*")

    for term in terms
        if length(term) == 0
            continue
        end
        # First, check for a match with total_tokens.
        found_total = nothing
        for (orig_idx, token) in tokens_total_with_idx
            if startswith(term, token)
                rem = term[length(token)+1:end]
                # Accept if nothing follows...
                if isempty(rem)
                    found_total = (orig_idx, token, rem)
                    break
                    # ...or if the remainder begins with "^" and the rest is a valid number or with ' and then possibly ^ and a number 
                elseif startswith(rem, "^")
                    num_part = rem[2:end]
                    try
                        parse(Int, num_part)
                        found_total = (orig_idx, token, rem)
                        break
                    catch e
                        error("Not a valid number $rem in $token.")
                    end
                end
            end
        end
        if found_total !== nothing
            (idx, token, cleaned) = found_total
            push!(output_total[idx], cleaned)
            continue  # Process next term.
        end

        found_end = nothing   # For tokens that should appear at the end.
        found_begin = nothing # For tokens that should appear at the beginning.

        # Check if term ends with any token from tokens_at_end.
        for (orig_idx, token) in tokens_at_end_with_idx
            if endswith(term, token)
                found_end = (orig_idx, token)
                break
            end
        end

        # Check if term starts with any token from tokens_at_begin.
        for (orig_idx, token) in tokens_at_begin_with_idx
            if startswith(term, token)
                found_begin = (orig_idx, token)
                break
            end
        end

        # The term must match exactly one type.
        if (found_end === nothing) && (found_begin === nothing)
            error("No matching token found for term: \"$term\"")
        elseif (found_end !== nothing) && (found_begin !== nothing)
            error("Term matches both beginning and ending tokens: \"$term\"")
        end

        if found_end !== nothing
            (idx, token) = found_end
            # Remove the token from the end.
            cleaned = term[1:end-length(token)]
            push!(output_at_end[idx], cleaned)
        else
            (idx, token) = found_begin
            # Remove the token from the beginning.
            cleaned = term[length(token)+1:end]
            push!(output_at_begin[idx], cleaned)
        end
    end

    return output_total, output_at_end, output_at_begin
end
function expstr_separate(expstr::String)::Tuple{String,Int}
    exp::Int = 1
    if occursin("^", expstr)
        expstr, b = split(expstr, "^")
        exp = parse(Int, b)
    end
    return expstr, exp
end

end