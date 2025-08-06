using Test
using QAlgebra


# 1) List all your CFunction subtypes here:
const CF_TYPES = (CAtom, CExp, CLog, CRational, CProd, CSum)


function gen_examples(::Type{T}, depth::Int, max_depth::Int) where T<:CFunction
    # -- one small collection of truly atomic CAtom examples:
    atomic = [
        CAtom(1,      [0,0]),
        CAtom(2,      [1,0]),
        CAtom(3//2,   [0,1]),
        CAtom(crationalize(1+2im), [1,1])
    ]
    # convenience:
    a = atomic[1]

    # at max depth, return exactly one “simple” example per type
    if depth >= max_depth
        if T === CAtom
            return [a]
        elseif T === CSum
            return [CSum([a, a])]
        elseif T === CProd
            return [CProd(a, a)]
        elseif T === CRational
            return [CRational(a, a)]
        elseif T === CExp
            return [CExp(a)]
        elseif T === CLog
            return [CLog(a)]
        else
            error("Unhandled type $T")
        end
    end

    # otherwise depth < max_depth: build composites
    if T === CAtom
        return atomic
    elseif T === CSum
        ex = CFunction[]
        for U in EXAMPLE_TYPES, V in EXAMPLE_TYPES
            e1 = gen_examples(U, depth+1, max_depth)[1]
            e2 = gen_examples(V, depth+1, max_depth)[1]
            push!(ex, CSum([e1,e2]))
        end
        return ex
    elseif T === CProd
        ex = CFunction[]
        for U in EXAMPLE_TYPES, V in EXAMPLE_TYPES
            e1 = gen_examples(U, depth+1, max_depth)[1]
            e2 = gen_examples(V, depth+1, max_depth)[1]
            push!(ex, CProd(e1, e2))
        end
        return ex
    elseif T === CRational
        ex = CFunction[]
        for U in EXAMPLE_TYPES, V in EXAMPLE_TYPES
            e1 = gen_examples(U, depth+1, max_depth)[1]
            e2 = gen_examples(V, depth+1, max_depth)[1]
            push!(ex, CRational(e1, e2))
        end
        return ex
    elseif T === CExp
        ex = CFunction[]
        for U in EXAMPLE_TYPES
            e1 = gen_examples(U, depth+1, max_depth)[1]
            push!(ex, CExp(e1))
        end
        return ex
    elseif T === CLog
        ex = CFunction[]
        for U in EXAMPLE_TYPES
            e1 = gen_examples(U, depth+1, max_depth)[1]
            push!(ex, CLog(e1))
        end
        return ex
    else
        error("Unhandled type $T")
    end
end

@testset "CFunctions Tests" begin
    # First check that we are checking all subtypes 
    list_of_subtypes = subtypes(CFunction)
    for T in list_of_subtypes
        @assert T in list_of_subtypes "Type $T not included in test -> add to test!"
    end

    max_depth = 2
    # collect by type
    examples_by_type = Dict{DataType, Vector{CFunction}}()
    for T in EXAMPLE_TYPES
        examples_by_type[T] = gen_examples(T, 0, max_depth)
    end
    # flatten
    all_ex = reduce(vcat, values(examples_by_type))

    # a fixed small var‐name list for to_string
    const VARS = ["x", "y"]

    @testset "CFunctions smoke tests" begin

        # string generation
        for ex in all_ex
            @test try stringer(ex)    true catch _ false end
            @test try to_string(ex, VARS)    true catch _ false end
        end

        # binary arithmetic on every pair
        for ex1 in all_ex, ex2 in all_ex
            @test try _ = ex1 + ex2; true  catch _ false end
            @test try _ = ex1 * ex2; true  catch _ false end
            @test try _ = ex1 / ex2; true  catch _ false end
        end

        # simplify, sorting, recursive_sort!
        for ex1 in all_ex, ex2 in all_ex
            @test try _ = simplify(ex1);                   true catch _ false end
            @test try _ = simplify(ex1 + ex2);             true catch _ false end
            @test try recursive_sort!(ex1);                true catch _ false end
            @test try _ = sort!([ex1, ex2]);               true catch _ false end
        end

    end
end

