export Ladder

"""
    Ladder() -> OperatorSet

Creates the OperatorSet for a bosonic mode using creation and annihilation operators (Ladder operators: $a^\dagger$, $a$).
"""
function Ladder()
    ops = ["'", ""]  # (Creation, Annihilation)
    base_ladder = [[0, 1]]
    function ladderstr2ind(str::String)::Vector{Tuple{Complex,Vector{Int}}}
        str, exp = expstr_separate(str)
        res = findfirst(==(str), ops)
        if res == nothing
            error("Invalid Ladder string: $str, must be one of $ops, or empty (may contain ^integer). ")
        end
        if res == 1
            return [(1.0 + 0im, [exp, 0])]
        else
            return [(1.0 + 0im, [0, exp])]
        end
    end
    function ladder_product_(a::Vector{Int}, b::Vector{Int})::Vector{Tuple{Complex,Vector{Int}}}
        # using a^n a' = a' a^n + n a^{n-1}
        if a[2] > 0 && b[1] > 0
            first = ladder_product_([a[1] + 1, a[2]], [b[1] - 1, b[2]])
            second = ladder_product_([a[1], a[2] - 1], [b[1] - 1, b[2]])
            for s in second
                push!(first, (a[2] * s[1], s[2]))
            end
            return first
        else
            return [(1.0 + 0im, [a[1] + b[1], a[2] + b[2]])]
        end
    end
    function ladder_product(a::Vector{Int}, b::Vector{Int})::Vector{Tuple{Complex,Vector{Int}}}
        return cleanup_terms(ladder_product_(a, b))
    end
    function ladder_dag(op::Vector{Int})::Vector{Tuple{Complex,Vector{Int}}}
        return [(1.0 + 0im, [op[2], op[1]])]
    end
    function ladderstr2ind(str::Vector{String})::Vector{Tuple{Complex,Vector{Int}}}
        if length(str) == 0
            return [(1.0 + 0im, [0, 0])]
        end
        new_exp::Vector{Tuple{Complex,Vector{Int}}} = ladderstr2ind(str[1])
        for s_str in str[2:end]
            s = ladderstr2ind(s_str)[1]
            new_new_exp::Vector{Tuple{Complex,Vector{Int}}} = []
            for exp in new_exp
                curr_terms = ladder_product(exp[2], s[2])
                for i in 1:length(curr_terms)
                    curr_terms[i] = (curr_terms[i][1] * s[1] * exp[1], curr_terms[i][2])
                end
                append!(new_new_exp, curr_terms)
            end
            new_exp = new_new_exp
        end
        return new_exp
    end
    function ladder2str(a::Vector{Int}, sym::String)::String
        curr_str = ""
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
        return curr_str
    end
    function ladder2str(a::Int, sym::String)::String
        if a == 1
            return ladder2str([1, 0], sym)
        else
            return ladder2str([0, 1], sym)
        end
    end
    function ladder2latex(a::Vector{Int}, sym::String, do_sigma::Bool=false)::String
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
    return OperatorSet("Ladder", false, 2, Int[0, 0], base_ladder, ops, ladder_product, ladder_dag, ladderstr2ind, ladder2str, ladder2latex)
end