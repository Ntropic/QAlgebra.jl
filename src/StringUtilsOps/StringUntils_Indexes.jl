# Processes strings of the forms:
#   --> "pref_{i,j,k}" and returns ("pref", ["i","j","k"]) 
#   --> "pref_i" and returns ("pref", ["i"]) 
#   --> "pref" and returns ("pref", [])
function underscore_index_split(s::AbstractString)
    # split at first "_" or first subscript char
    pos_us  = findfirst(==('_'), s)
    pos_sub = findfirst(c -> c in SUB_CHARS, s)
    startpos = pos_us === nothing ? pos_sub :
               pos_sub === nothing ? pos_us :
               (pos_us < pos_sub ? pos_us : pos_sub)
    if startpos === nothing
        return s, ""
    end
    ch = s[startpos]
    if ch == '_'
        head = startpos == firstindex(s) ? "" : s[1:prevind(s, startpos)]
        tail = startpos == lastindex(s)  ? "" : s[nextind(s, startpos):end]
        return head, tail
    else
        head = startpos == firstindex(s) ? "" : s[1:prevind(s, startpos)]
        tail = s[startpos:end]
        return head, tail
    end
end

_underscore_strip_braces(s::AbstractString) = begin
    t = strip(s)
    (startswith(t, "{") && endswith(t, "}")) ? t[2:end-1] : t
end

# --- unsubscribe using your REVERSE_SUBSCRIPT ---
function _underscore_unsubscript(s::AbstractString)
    out = IOBuffer()
    for c in s
        if haskey(REVERSE_SUBSCRIPT, c)
            write(out, REVERSE_SUBSCRIPT[c])
        else
            write(out, c)
        end
    end
    return String(take!(out))
end

# --- tokenize payload into components like ["i1", "i2"] ---
function _underscore_tokenize_indices(payload::AbstractString)
    p = _underscore_strip_braces(payload)
    p = replace(p, r"[_\s]+" => "")
    p = _underscore_unsubscript(p)
    parts = split(p, ',')
    comps = String[]
    for part in parts
        part = strip(part)
        isempty(part) && continue
        current = ""
        for c in part
            if isletter(c)
                if !isempty(current)
                    push!(comps, current)
                end
                current = string(c)
            elseif isdigit(c)
                current *= string(c)
            else
                if !isempty(current)
                    push!(comps, current)
                    current = ""
                end
            end
        end
        if !isempty(current)
            push!(comps, current)
        end
    end
    return comps
end

function normalize_underscore_indices(s::AbstractString)
    prefix, payload = underscore_index_split(s)
    if isempty(payload)
        return prefix, String[]
    end
    comps = _underscore_tokenize_indices(payload)
    return prefix, comps
end

"""
    format_normalized_indices(tokens::Vector{String}; do_latex::Bool=false)

Render index tokens produced by `normalize_underscore_indices` into a subscript
suffix suitable for plain text or LaTeX output.
"""
function format_normalized_indices(tokens::Vector{String}; do_latex::Bool=false)::String
    isempty(tokens) && return ""
    index_str_raw = join(tokens, ",")
    if do_latex
        return "_{" * index_str_raw * "}"
    else
        return str2sub(index_str_raw)
    end
end


# ==========> Matching components for SubSpaces <================================================

function split_index(idx::AbstractString)::Tuple{String, Int}
    pos = findfirst(isdigit, idx)
    if pos === nothing
        return idx, 0
    end
    letter = idx[1:prevind(idx, pos)]
    number = parse(Int, idx[pos:end])
    return letter, number
end

function find_allowed_index(idx::AbstractString, allowed::Vector{Vector{String}})::Tuple{Int, Int}
    for (outer_i, group) in enumerate(allowed)
        for (inner_i, name) in enumerate(group)
            if idx == name
                return outer_i, inner_i
            end
        end
    end
    error("idx not found in allowed indexes: $allowed")
end

function match_components(allowed::Vector{Vector{String}}, s::AbstractString)
    prefix, comps = normalize_underscore_indices(s)

    matches = Vector{Tuple{Int,Int, Int}}()
    for comp in comps
        base, num = split_index(comp)
        outer, inner = find_allowed_index(base, allowed)
        push!(matches, (outer, inner, num))
    end

    return prefix, matches
end
