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
    @testset "CFunction Evaluate Tests" begin
        for ex in all_ex
            @test_succeeds evaluate(ex, xv)      "evaluate($ex, xv) failed"
            pv = ParameterValues(ex.param_info)
            for idx in 1:length(ex.param_info.params_name)
                set_param!(pv, idx, xv[(idx-1) % length(xv) + 1])
            end
            @test_succeeds evaluate(ex, pv)   "evaluate($ex, pv) failed"
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

@testset "Indexed evaluation" begin
    sub_def = SubSpaceDefinitions(i=Ensemble(2, 0, QubitPauli()))
    op_def = OperatorDefinitions()
    gamma_dist = QUniform(-1.0, 1.0, 32)
    param_def = ParameterDefinitions("alpha", "gamma_i" => gamma_dist)
    qspace = QSpace(sub_def, op_def, param_def)
    pinfo = qspace.param_info

    alpha_idx = findfirst(==(Symbol(:alpha)), pinfo.inner_labels_symbols_flat)
    gamma_idx = findfirst(!=0, pinfo.indexed_parameter_indexes)
    @test !isnothing(alpha_idx)
    @test !isnothing(gamma_idx)
    alpha_idx = alpha_idx::Int
    gamma_idx = gamma_idx::Int
    group_dists = pinfo.group_distributions
    @test group_dists[alpha_idx] === nothing
    @test group_dists[gamma_idx] isa QDistribution
    group_funcs = pinfo.group_functions
    @test group_funcs[alpha_idx] === nothing
    @test group_funcs[gamma_idx] === nothing

    exps_alpha = zeros(Int, pinfo.dims)
    exps_alpha[alpha_idx] = 1
    atom_alpha = CAtom(pinfo, exps_alpha)
    @test !has_indexed_parameters(atom_alpha)
    @test evaluate(atom_alpha, ones(pinfo.dims)) == 1.0

    exps_gamma = zeros(Int, pinfo.dims)
    exps_gamma[gamma_idx] = 1
    atom_gamma = CAtom(pinfo, exps_gamma)
    @test has_indexed_parameters(atom_gamma)
    @test_throws ErrorException evaluate(atom_gamma, ones(pinfo.dims))

    concrete = ConcreteIndexes(pinfo)
    concrete.indexes[1] = [2, 1]
    atom_gamma_indexed = CAtomIndexed(pinfo, exps_gamma, concrete)

    pv = ParameterValues(pinfo)
    set_param!(pv, :alpha, 1.0)
    set_time!(pv, 0.0)

    @test evaluate(atom_alpha, pv) == 1.0

    gamma_group = findfirst(==(Symbol("gamma")), pinfo.outer_labels_symbols)
    gamma_params = pinfo.params_by_group[gamma_group]
    gamma_values = [5.0, 6.0, 7.0]
    for (val_idx, param_idx) in enumerate(gamma_params)
        set_param!(pv, param_idx, gamma_values[val_idx])
    end
    @test evaluate(atom_gamma_indexed, pv) == 6.0

    indexed_sum = Indexed(CSum(pinfo, [atom_alpha, atom_gamma]), concrete)
    @test indexed_sum isa CSum
    @test all(term -> term isa CAtomIndexed, indexed_sum.expr)
    @test evaluate(indexed_sum, pv) isa Float64

    prod_expr = CProd(pinfo, ComplexRational(1, 0, 1), [atom_alpha, atom_gamma], Val(:nosimp))
    indexed_prod = Indexed(prod_expr, concrete)
    @test indexed_prod isa CProd
    @test all(term -> term isa CAtomIndexed, indexed_prod.expr)
    @test evaluate(indexed_prod, pv) isa Float64

    vec_expr = CVector(pinfo, [atom_alpha, atom_gamma])
    indexed_vec = Indexed(vec_expr, concrete)
    @test indexed_vec isa CVector
    @test all(entry -> entry isa CAtomIndexed, indexed_vec.expr)

    indexed_again = Indexed(atom_gamma_indexed, concrete)
    @test indexed_again isa CAtomIndexed
    @test indexed_again.indexes === atom_gamma_indexed.indexes

    @test_throws ArgumentError Indexed(atom_alpha, [[1, 2]])
end
