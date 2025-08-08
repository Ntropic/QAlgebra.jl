# Test Helper Scripts 

using ComplexRationals
in_vsc() = any(startswith(key, "VSCODE_") for key in keys(ENV))
function line(msg::String=""; tail::String="=", heads::Vector{String}=[">", "<"], n::Int=-1)
    if n < 0
        if !in_vsc()
            n, _ = displaysize(stdout)
        else 
            n = 70
        end
        if n == 0
            n = 70
        end
    end
    m = length(msg)
    if m == 0
        println(tail^n)
        return
    end
    n2 = Int(floor((n - m - 4)/2)) 
    if n2 < 0
        n2 = 0
    end
    final_str = join([tail^n2, heads[1], " ", msg, " ", heads[2], tail^n2], "")
    if length(final_str) < n 
        final_str = final_str * tail 
    end
    println(final_str) 
    return
end

macro test_succeeds(expr, msg_func="")
    quote
        @test begin
            try
                $(esc(expr))
                true
            catch err
                line($(esc(msg_func)))
                line(string(typeof(err)))
                Base.showerror(stdout, err)
                line()
                println()
                rethrow(err)
                false
            end
        end
    end
end

### Combinatorics for AST Trees 
function CFunctions_gen_examples(::Type{T}, depth::Int, max_depth::Int, same_a::Bool=false) where T<:CFunction
    # -- one small collection of truly atomic CAtom examples:
    atomic = [
        CAtom(2,      [0,0]),
        CAtom(2,      [1,0]),
        CAtom(3//2,   [0,1]),
        CAtom(crationalize(1+2im), [1,1])
    ]
    # convenience:
    a = atomic[1]
    a2 = atomic[3]
    if same_a 
        a2 = atomic[2]
    end

    # at max depth, return exactly one “simple” example per type
    if depth >= max_depth
        if T === CAtom
            return [a]
        elseif T === CSum
            return [CSum([a, a2]), CSum([a, a2, a])]
        elseif T === CProd
            return [CProd([a, a2])]
        elseif T === CRational
            return [CRational(a, a2), CRational(a,a)]
        elseif T === CExp
            return [CExp(a), CExp(a2)]
        elseif T === CLog
            return [CLog(a), CLog(a2)]
        else
            error("Unhandled type $T")
        end
    end

    # otherwise depth < max_depth: build composites
    if T === CAtom
        return atomic
    elseif T === CSum
        ex = CFunction[]
        for U in CF_TYPES, V in CF_TYPES
            if U !== CSum && V !== CSum 
                e1 = CFunctions_gen_examples(U, depth+1, max_depth, same_a)[1]
                e2 = CFunctions_gen_examples(V, depth+1, max_depth, same_a)[1]
                push!(ex, CSum([e1,e2]))
            end
        end
        return ex
    elseif T === CProd
        ex = CFunction[]
        for U in CF_TYPES, V in CF_TYPES
            e1 = CFunctions_gen_examples(U, depth+1, max_depth, same_a)[1]
            e2 = CFunctions_gen_examples(V, depth+1, max_depth, same_a)[1]
            push!(ex, CProd([e1, e2]))
        end
        return ex
    elseif T === CRational
        ex = CFunction[]
        for U in CF_TYPES, V in CF_TYPES
            e1 = CFunctions_gen_examples(U, depth+1, max_depth, same_a)[1]
            e2 = CFunctions_gen_examples(V, depth+1, max_depth, same_a)[1]
            push!(ex, CRational(e1, e2))
        end
        return ex
    elseif T === CExp
        ex = CFunction[]
        for U in CF_TYPES
            e1 = CFunctions_gen_examples(U, depth+1, max_depth, same_a)[1]
            push!(ex, CExp(e1))
        end
        return ex
    elseif T === CLog
        ex = CFunction[]
        for U in CF_TYPES
            e1 = CFunctions_gen_examples(U, depth+1, max_depth, same_a)[1]
            push!(ex, CLog(e1))
        end
        return ex
    else
        error("Unhandled type $T")
    end
end
