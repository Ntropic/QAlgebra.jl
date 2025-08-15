export Ladder

@doc raw"""
    Ladder() -> OperatorSet

Creates the OperatorSet for a bosonic mode using creation and annihilation operators (Ladder operators: ``a^\dagger``, ``a``).
"""
function Ladder()
    ops = ["'"]  # (Creation, Annihilation) -> removed annihilation
    base_ladder = [[0, 1]]
    non_base_ops::Dict{String, Vector{Tuple{ComplexRational, Is}}} = Dict("n" => [(ComplexRational(1,0,1), [1,1])])
    function ladderstr2ind(str::String)::Vector{Tuple{ComplexRational,Vector{Int}}}
        str, exp = expstr_separate(str)
        res = findfirst(==(str), ops)
        if isnothing(res)
            error("Invalid Ladder string: $str, must be one of $ops, or empty (may contain ^integer). ")
        end
        if res == 1
            return [(ComplexRational(1,0,1), [exp, 0])]
        else
            return [(ComplexRational(1,0,1), [0, exp])]
        end
    end
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

    return OperatorSet("Ladder", "Boson", 2, Int[0, 0], base_ladder, non_base_ops, ops, ladder_product, ladder_dag, ladder2str, ladder2latex, laddercommutes) 
end
