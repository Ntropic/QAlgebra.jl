export Ladder

function cleanup_terms(terms::Vector{Tuple{T,S}})::Vector{Tuple{T,S}} where {T<:Number,S}
    # 1) sort once by index
    sort!(terms, by = x -> x[2])
    # 2) prealloc output to worst‑case length and scan in one pass
    n = length(terms)
    T0 = typeof(terms[1][1])
    S0 = typeof(terms[1][2])
    cleaned = Vector{Tuple{T0,S0}}(undef, n)
    cnt = 0
    i = 1
    @inbounds while i ≤ n
        sumc, idx = terms[i]           # destructure once
        j = i + 1
        # inner loop: accumulate identical idx
        @inbounds while j ≤ n && terms[j][2] == idx
            sumc += terms[j][1]
            j += 1
        end

        # push nonzero
        if sumc != zero(T0)
            cnt += 1
            cleaned[cnt] = (sumc, idx)
        end

        i = j
    end
    resize!(cleaned, cnt)                # trim unused slots
    return cleaned
end


@doc raw"""
    Ladder(; max_magnitude::Int=-1) -> OperatorSet

Creates the OperatorSet for a bosonic mode using creation and annihilation operators (Ladder operators: ``a^\dagger``, ``a``).
Provide `max_magnitude` to cap the highest occupation per index; the operator set records both that bound and the induced maximum magnitude. Pass `-1` (default) for unbounded.
"""
function Ladder(; max_magnitude::Int=-1)
    ops = [""]  # (Creation, Annihilation) -> removed annihilation
    base_ladder = [[0, 1]]
    max_vec = max_magnitude < 0 ? Int[-1, -1] : Int[max_magnitude, max_magnitude]
    function ladder_product_(a::Vector{Int}, b::Vector{Int})::Vector{Tuple{ComplexRational,Vector{Int}}}
        # using a^n a' = a' a^n + n a^{n-1}
        if a[2] > 0 && b[1] > 0
            first = ladder_product_([a[1] + 1, a[2]], [b[1] - 1, b[2]])
            second = ladder_product_([a[1], a[2] - 1], [b[1] - 1, b[2]])
            for s in second
                push!(first, (a[2] * s[1], s[2]))
            end
            return first
        else
            return [(ComplexRational(1,0,1), [a[1] + b[1], a[2] + b[2]])]
        end
    end
    function ladder_product(a::Vector{Int}, b::Vector{Int})::Vector{Tuple{ComplexRational,Vector{Int}}}
        return cleanup_terms(ladder_product_(a, b))
    end
    function ladder_dag(op::Vector{Int})::Vector{Tuple{ComplexRational,Vector{Int}}}
        return [(ComplexRational(1,0,1), [op[2], op[1]])]
    end
    function ladder2str(a::Vector{Int}, sym::String; formatted::Bool=true)::String
        curr_str = ""
        if formatted
            if a[1] > 0
                curr_str *= sym * "†"
                if a[1] > 1
                    curr_str *= str2sup(string(a[1]))
                end
            end
            if a[2] > 0
                curr_str *= sym
                if a[2] > 1
                    curr_str *= str2sup(string(a[2]))
                end
            end
        else
            if a[1] > 0
                curr_str *= sym * "'"
                if a[1] > 1
                    curr_str *= "^"* string(a[1])
                end
            end
            if a[2] > 0
                curr_str *= sym
                if a[2] > 1
                    curr_str *= "^"* str2supstring(a[2])
                end
            end
        end
        return curr_str
    end
    function ladder2latex(a::Vector{Int}, sym::String)::String
        curr_str = ""
        if a[1] > 0
            if a[1] > 1
                curr_str *= raw"\hat{" * sym * raw"}^{\dagger " * string(a[1]) * "}"
            else
                curr_str *= raw"\hat{" * sym * raw"}^\dagger"
            end
        end
        if a[2] > 0
            curr_str *= raw"\hat{" * sym * "}"
            if a[2] > 1
                curr_str *= "^" * string(a[2]) * " "
            end
        end
        return curr_str
    end
    function laddercommutes(a::Vector{Int}, b::Vector{Int})::Bool
        p, q = a
        r, s = b
        #       a neutral,        b neutral,       both only creation, both annihilation, both number operators 
        return (p==0 && q==0) || (r==0 && s==0) || (q==0 && s==0)  || (p==0 && r==0) ||  (p==q && r==s)
    end
    function ladder_operators_magnitude(op::Is)::Int
        return op[1] + op[2]
    end

    return OperatorSet("Ladder", "Boson", 2, Int[0, 0], base_ladder, ops, ladder_product, ladder_dag, ladder2str, ladder2latex;
                      commutes=laddercommutes, operator_magnitude=ladder_operators_magnitude,
                      min_ints=Int[0, 0], max_ints=max_vec, max_magnitude=max_magnitude) 
end
