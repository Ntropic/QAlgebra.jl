@testset "QAlgebra Tests" begin

    # === SETUP ===
    subspace_def = SubSpaceDefinitions(h=QubitPM(), i=Ensemble(3, 2, QubitPauli()), b=Ladder())
    op_def = OperatorDefinitions()
    param_def = ParameterDefinitions("alpha", "beta(t)", "gamma_i", "delta_i")
    qspace = QSpace(subspace_def, op_def, param_def)

    xi, yi, zi = base_operators(qspace, "i", by_ensemble=false)
    xj, yj, zj = base_operators(qspace, "j", by_ensemble=false)
    xk, yk, zk = base_operators(qspace, "k", by_ensemble=false)
    ph, mh, zh = base_operators(qspace, "h", by_ensemble=false)
    b = base_operators(qspace, "b", by_ensemble=false)
    I = base_operators(qspace, "I")
    alpha = base_operators(qspace, "alpha")
    beta = base_operators(qspace, "beta")
    gamma_i, gamma_j, gamma_k = base_operators(qspace, "gamma", by_ensemble=false)
    delta_i, delta_j, delta_k = base_operators(qspace, "delta", by_ensemble=false)

    # === TESTS ===

    @testset "QSpace Construction" begin
        @test qspace isa QSpace
        ensemble_cfg = qspace.subspaces[2].ensemble
        @test ensemble_cfg !== nothing
        @test ensemble_cfg.num_modes == -1
        @test :gamma in ensemble_cfg.parameter_groups
        @test length(qspace.ensembles) == 1
        @test qspace.ensembles[1] === ensemble_cfg
        @test ensemble_cfg.qspace_ref !== nothing
        @test ensemble_cfg.qspace_ref.value === qspace
    end

    @testset "Ensemble Naming" begin
        sub_def = SubSpaceDefinitions(i=Ensemble(3,3,QubitPauli("sigma")), j=QubitPM("beta"))
        q_tmp = QSpace(sub_def, OperatorDefinitions(), ParameterDefinitions())
        idx_i = findfirst(s -> s.key_symbol == :i, q_tmp.subspaces)
        idx_j = findfirst(s -> s.key_symbol == :j, q_tmp.subspaces)
        @test q_tmp.subspaces[idx_i].keys[1:3] == ["i0", "i1", "i2"]
        @test q_tmp.subspaces[idx_j].keys == ["j"]
    end

    @testset "Base Operators Extraction" begin
        alpha_expr = base_operators(qspace, "alpha")
        beta_expr = base_operators(qspace, "beta")
        ops = base_operators(qspace, ["b", "x_i"])
        @test alpha_expr isa QExpr
        @test beta_expr isa QExpr
        @test all(x -> x isa QExpr, ops)
        @test xi isa QExpr
        @test yj isa QExpr
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

        exp_bs = Expectation(Bs)
        @test all(term -> !(term isa QAtomProduct) || term.braket, exp_bs.terms)
    end

    @testset "Differentiation Tests" begin
        diff_eq = d_dt(zi, alpha^2)
        @test diff_eq isa diffQEq
        @test diff_eq.left_hand_side.braket
    end

    @testset "Pauli Algebra Rules" begin
        @test xi * yi == -yi * xi
        @test xi * yi == im * zi
        @test xi * yj == yj * xi
        @test xi * xi == I
    end
    @testset "Atom Product Ordering" begin
        base_atoms = base_operators(qspace, "i", by_ensemble=false)
        xi_expr = base_atoms[1]
        yi_expr = base_atoms[2]
        xi_term = xi_expr.terms[1].expr[1]
        yi_term = yi_expr.terms[1].expr[1]
        prod = QAtomProduct(qspace, QAtom[yi_term, xi_term])
        ordered = OrderedQAtomProduct(prod)
        @test ordered isa QAtomOrdered
        @test ordered.expr[1] == xi_term
        @test ordered.expr[2] == yi_term
        @test permutation(ordered) == [2, 1]

        expr_unordered = QExpr(qspace, QComposite[prod], Val(:nosimp))
        expr_ordered = OrderedQExpr(expr_unordered)
        @test expr_ordered isa QExpr
        ordered_term = expr_ordered.terms[1]
        @test ordered_term isa QAtomOrdered
        @test ordered_term.expr[1] == xi_term
        @test ordered_term.expr[2] == yi_term
        @test permutation(ordered_term) == [2, 1]

        rhs_simple = QExpr(qspace, QComposite[prod], Val(:nosimp))
        diff_unordered = diffQEq(qspace, prod, rhs_simple, Val(:raw))
        diff_ordered = OrderedDiffQEq(diff_unordered)
        @test diff_ordered isa diffQEqOrdered
        @test diff_ordered.left_hand_side isa QAtomOrdered
        @test diff_ordered.left_hand_side.expr[1] == xi_term
        @test permutation(diff_ordered.left_hand_side) == [2, 1]

        set_unordered = diffQEqSet([diff_unordered])
        set_ordered = OrderedDiffQEqSet(set_unordered)
        @test set_ordered isa diffQEqSetOrdered
        @test set_ordered.equations[1] == diff_ordered

        set_from_constructor = diffQEqSet([diff_unordered])
        @test set_from_constructor isa diffQEqSetOrdered
        @test set_from_constructor.equations[1] == diff_ordered
    end
    @testset "Concrete Index Attachments" begin
        param_info = qspace.param_info
        lengths = param_info.how_many_by_ensemble
        concrete = ConcreteIndexes(param_info)
        for (ensemble_idx, len) in enumerate(lengths)
            for inner in 1:len
                concrete.indexes[ensemble_idx][inner] = inner
            end
        end

        indexed_param = findfirst(!iszero, param_info.indexed_parameter_indexes)
        @test indexed_param !== nothing
        idx_val = indexed_param::Int
        tuples = parameter_index_tuples(param_info, idx_val)
        @test !isempty(tuples)
        @test all(t -> 1 ≤ t[1] ≤ length(lengths), tuples)
        for (ensemble, inner) in tuples
            @test concrete.indexes[ensemble][inner] == inner
        end

        exponents = zeros(Int, param_info.dims)
        exponents[idx_val] = 1
        coeff_atom = CAtom(param_info, exponents)
        indexed_atom = with_concrete_indexes(coeff_atom, concrete)
        @test indexed_atom isa CAtomIndexed
        @test indexed_atom.indexes.indexes == concrete.indexes

        base_atoms = base_operators(qspace, "i", by_ensemble=false)
        xi_term = base_atoms[1].terms[1].expr[1]
        yi_term = base_atoms[2].terms[1].expr[1]
        ordered = OrderedQAtomProduct(QAtomProduct(qspace, QAtom[yi_term, xi_term]))
        indexed_qatom = QAtomIndexed(ordered, concrete)
        @test indexed_qatom isa QAtomIndexed
        @test indexed_qatom.concrete_indexes.indexes == concrete.indexes
        @test indexed_qatom.ensemble_indexes == ordered.ensemble_indexes

        @test_throws ErrorException ConcreteIndexes(param_info, [[1]])
    end
    @testset "PM Basis Rules" begin
        @test mh * ph == 1 / 2 * (I - zh)
        @test ph * mh == 1 / 2 * (zh + I)
        @test ph * ph == 0*I
        @test ph' == mh
    end

    @testset "Ladder Operator Rules" begin
        comm_expr = b' * b - b * b'
        @test comm_expr isa QExpr
    end

    @testset "NeqConstraint Integration" begin
        expr_sum = xi + yi
        constraint = neq(:l, :m)

        sum_direct = ∑([:l, :m], expr_sum, constraint)
        sum_base = Base.sum([:l, :m], expr_sum, constraint)

        @test sum_direct isa QExpr
        @test sum_base isa QExpr
        base_cons = sum_base.terms[1].blocks[1].constraints
        @test sum_direct.terms[1].blocks[1].constraints == base_cons

        qsum = sum_direct.terms[1]
        block = qsum.blocks[1]
        @test length(block.indexes) == 2
        lhs_inner = block.ensemble_indexes[1]
        rhs_inner = block.ensemble_indexes[2]
        @test block.constraints[1][lhs_inner]
        @test block.constraints[2][rhs_inner]
        @test !block.constraints[1][rhs_inner]
        @test !block.constraints[2][lhs_inner]

        non_sum_constraint = neq(:l, :i)
        sum_with_fixed = ∑([:l], expr_sum, non_sum_constraint)
        single_block = sum_with_fixed.terms[1].blocks[1]
        non_sum_idx = SubSpaceIndex(:i, qspace.subspace_info)
        @test !single_block.constraints[1][non_sum_idx.inner]

        @test_throws ErrorException ∑([:l], expr_sum, constraint)
        @test_throws ErrorException ∑([:l, :m], expr_sum, neq(:l, :l))
        @test_throws ErrorException ∑([:l, :m], expr_sum, neq(:l, :h))
    end

end
