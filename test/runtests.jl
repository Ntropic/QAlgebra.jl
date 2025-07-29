using Test
using QAlgebra

@testset "QAlgebra Tests" begin

    # === SETUP ===
    qspace = StateSpace("alpha", "beta(t)", "gamma_i", "delta_i", h=QubitPM(), i=(3, QubitPauli()), b=Ladder())

    var_dict, op_dict, abstract_dict = base_operators(qspace)

    xi, yi, zi, _, _ = base_operators(qspace, "i", do_dict=false)
    xj, yj, zj, pj, mj = base_operators(qspace, "j", do_dict=false)
    xk, yk, zk, _, _ = base_operators(qspace, "k", do_dict=false)
    ph, mh, zh, _, _ = base_operators(qspace, "h", do_dict=false)
    b, n = base_operators(qspace, "b", do_dict=false)
    I = base_operators(qspace, "I")
    var_dict2 = base_operators(qspace, "vars")
    alpha = base_operators(qspace, "alpha")
    beta = base_operators(qspace, "beta")
    gamma_i, gamma_j, gamma_k = base_operators(qspace, "gamma", do_dict=false)
    delta_i, delta_j, delta_k = base_operators(qspace, "delta", do_dict=false)

    # === TESTS ===

    @testset "StateSpace Construction" begin
        @test qspace isa StateSpace
    end

    @testset "Base Operators Extraction" begin
        @test var_dict isa Dict
        @test var_dict2 isa Dict
        @test op_dict isa Dict
        @test abstract_dict isa Dict
        @test haskey(var_dict, "alpha")
        @test haskey(var_dict2, "alpha")
        @test haskey(op_dict, "b")
        @test haskey(op_dict, "x_i")
        @test haskey(abstract_dict, "A")
        @test xi isa QExpr
        @test pj isa QExpr
        @test b isa QExpr
        @test alpha isa QExpr
        @test I isa QExpr
    end

    @testset "Simple Expressions" begin
        As = 2 * alpha * im * xi
        Bs = alpha * (Dag(b) * xi * yi)
        @test As isa QExpr
        @test Bs isa QExpr

        expr1 = 2 * alpha * im * zi
        expr2 = 2 * alpha *  xi * yi 
        @test expr1 == expr2
    end

    @testset "Sum and Nested Expressions" begin
        qsum_expr = Sum(["j", "k"], alpha * gamma_i * gamma_k * delta_k * xi * yi)
        @test qsum_expr isa QExpr
        @test qsum_expr.terms[1] isa QSum

        flat_expr = flatten(qsum_expr)
        @test flat_expr isa QExpr

        neq_expr = neq(qsum_expr)
        @test neq_expr isa QExpr
    end

    @testset "Differentiation Tests" begin
        diff_eq = d_dt(zi, alpha^2)
        @test diff_eq isa diff_QEq
    end

    @testset "Commutator and Simplify Tests" begin
        a_term = alpha * beta^2 * xi * yi * Dag(b) * b
        b_term = alpha^2 * beta * zi * Dag(b)
        comm = simplify(Commutator(a_term, b_term))
        @test comm isa QExpr

        messy = alpha * xi + alpha^2 * zi + alpha * xi
        expected = simplify(alpha^2 * zi + 2 * alpha * xi )
        @test simplify(simplify(messy)) == expected
    end

    @testset "Pauli Algebra Rules" begin
        @test xi * yi == -yi * xi
        @test xi * yi == im * zi
        @test xi * yj == yj * xi
        @test xi * xi == I
    end
    @testset "PM Basis Rules" begin
        @test mh * ph == 1 / 2 * (I - zh)
        @test ph * mh == 1 / 2 * (zh + I)
        @test ph * ph == 0*I
        @test ph' == mh
    end

    @testset "Ladder Operator Rules" begin
        @test b' * b - b * b' == -1*I
    end

end
