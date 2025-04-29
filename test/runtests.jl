using Test
using qAlgebra

@testset "qAlgebra Tests" begin

    # === SETUP ===
    qs = StateSpace("alpha", "beta(t)", "gamma_i", "delta_i", h=QubitPM(), i=(3, QubitPauli()), b=Ladder())

    var_dict, op_dict = base_operators(qs)

    xi, yi, zi = base_operators("i", qs)
    xj, yj, zj = base_operators("j", qs)
    xk, yk, zk = base_operators("k", qs)
    ph, mh, zh = base_operators("h", qs)
    I = base_operators("I", qs)
    b = base_operators("b", qs)
    alpha, beta, gamma_i, gamma_j, gamma_k, delta_i, delta_j, delta_k = base_operators("vars", qs)

    # === TESTS ===

    @testset "StateSpace Construction" begin
        @test qs isa StateSpace
    end

    @testset "Base Operators Extraction" begin
        @test var_dict isa Dict
        @test op_dict isa Dict
        @test haskey(var_dict, "alpha")
        @test haskey(op_dict, "b")
        @test xi isa qEQ
        @test b isa qEQ
    end

    @testset "Simple Expressions" begin
        A = 2 * alpha * im * xi
        B = alpha * (Dag(b) * xi * yi)
        @test A isa qEQ
        @test B isa qEQ

        expr1 = 2 * alpha * im * xi
        expr2 = 2 * alpha * im * xi
        @test expr1 == expr2
    end

    @testset "Sum and Nested Expressions" begin
        qsum_expr = Sum(["j", "k"], alpha * gamma_i * gamma_k * delta_k * xi * yi)
        @test qsum_expr isa qEQ
        @test qsum_expr.terms[1] isa qSum

        flat_expr = flatten(qsum_expr)
        @test flat_expr isa qEQ

        neq_expr = neq(qsum_expr)
        @test neq_expr isa qEQ
    end

    @testset "Differentiation Tests" begin
        diff_eq = d_dt(zi, alpha^2)
        @test diff_eq isa diff_qEQ
    end

    @testset "Commutator and Simplify Tests" begin
        a_term = alpha * beta^2 * xi * yi * Dag(b) * b
        b_term = alpha^2 * beta * zi * Dag(b)
        comm = simplify(Commutator(a_term, b_term))
        @test comm isa qEQ

        messy = alpha * xi + alpha^2 * zi + alpha * xi
        expected = 2 * alpha * xi + alpha^2 * zi
        @test simplify(messy) == expected
    end

    @testset "Pauli Algebra Rules" begin
        @test xi * yi == -yi * xi
        @test xi * yi == im * zi
        @test xi * yj == yj * xi
        @test xi * xi == 1
    end
    @testset "PM Basis Rules" begin
        @test mh * ph == 1 / 2 * (I - zh)
        @test ph * mh == 1 / 2 * (zh + I)
        @test ph * ph == 0
        @test ph' == mh
    end

    @testset "Ladder Operator Rules" begin
        @test b' * b - b * b' == -1
    end

end
