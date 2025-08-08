CF_TYPES = (CAtom, CExp, CLog, CRational, CProd, CSum)
VARS = ["x", "y"]

@testset "CFunction Tests" begin
    max_depth = 2
    # collect by type
    examples_by_type = Dict{DataType, Vector{CFunction}}()
    @testset "CFunction Examples Generation Tests" begin
        for T in CF_TYPES
            examples_by_type[T] = CFunctions_gen_examples(T, 0, max_depth)
        end
    end
    #line("Finished creating examples of all types." )

    # flatten
    all_ex = reduce(vcat, values(examples_by_type))

    # a fixed small var‐name list for to_string


    # string generation
    @testset "CFunction String Generation Tests" begin
        for ex in all_ex
            @test_succeeds to_string(ex, VARS)   "to_string($ex) failed"
            @test_succeeds to_string(ex, VARS, do_latex=true)   "to_string($ex, VARS, do_latex=true) failed"
        end
    end
    #line("Finished creating strings and LaTeXStrings of all examples.")

    # binary arithmetic on every pair
    @testset "CFunction Arithmetic Tests" begin
        for ex1 in all_ex, ex2 in all_ex
            @test_succeeds ex1 + ex2  "$ex1 + $ex2 failed"
            @test_succeeds ex1 - ex2  "$ex1 - $ex2 failed"
            @test_succeeds ex1 * ex2  "$ex1 * $ex2 failed"
            if !iszero(ex2) 
                @test_succeeds ex1 / ex2  "$ex1 / $ex2 failed"
            end
        end
    end

    # simplify, sorting, recursive_sort!
    @testset "CFunction Simplification Tests" begin
        for ex1 in all_ex, ex2 in all_ex
            @test_succeeds QAlgebra.CFunctions.simplify(ex1)    "simplify($ex1) failed"         
        end
    end

    # prepare a test point (all examples use 2 variables by construction)
    xv = [1.3, 0.7]
    xpows = build_xpows(xv, max_exponents(all_ex))

    #line("Testing evaluate(..., xv) and evaluate(..., xpows)")
    @testset "CFunction Evaluate Tests" begin
        for ex in all_ex
            @test_succeeds evaluate(ex, xv)      "evaluate($ex, xv) failed"
            @test_succeeds evaluate(ex, xpows)   "evaluate($ex, xpows) failed"
        end
    end
    #line("Testing expand modes")

    @testset "CFunction Expand Tests" begin
        for ex in all_ex
            # Taylor‐expand all CExp up to order 2
            @test_succeeds expand(ex, :Taylor, CExp, 2)     "expand(:Taylor, CExp) on $ex failed"
            # Taylor‐expand all CLog up to order 2
            @test_succeeds expand(ex, :Taylor, CLog, 2)     "expand(:Taylor, CLog) on $ex failed"
            # Distribute rationals over sums
            @test_succeeds expand(ex, :Rational, CRational) "expand(:Rational, CRational) on $ex failed"
            # Apply algebraic log rules
            @test_succeeds expand(ex, :Log, CLog)           "expand(:Log, CLog) on $ex failed"
        end
    end
end
