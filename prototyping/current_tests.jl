using Pkg
Pkg.instantiate()
using QAlgebra

using BenchmarkTools

subspace_def = SubSpaceDefinitions(h=QubitPM("beta"), i=(3, 3, QubitPauli("sigma")), b=Ladder())
op_def = OperatorDefinitions("A(i,t)", "B(U,H,i)")
var_def = ParameterDefinitions("alpha", "beta(t)", "gamma_i", "delta_{i,j}(t)")
qspace = StateSpace(subspace_def, op_def, var_def, max_t_ind=0)


alpha, beta, gamma, delta = base_operators(qspace, ["alpha", "beta", "gamma", "delta"], do_fun=true)
t = base_operators(qspace, :t)
ph, mh, zh = base_operators(qspace, "h")
sigma = base_operators(qspace, "i", by_ensemble=true, do_fun=true)
xi,yi,zi = base_operators(qspace, "i")
xj, yj, zj = base_operators(qspace, "j")
xk, yk, zk = base_operators(qspace, "k")
xl, yl, zl = base_operators(qspace, "l")
xm, ym, zm = base_operators(qspace, "m")
xn, yn, zn = base_operators(qspace, "n")
b = base_operators(qspace, "b")
I = base_operators(qspace, "I")
A = base_operators(qspace, "A", do_fun=true)
B = base_operators(qspace, "B", do_fun=true)

# See if this works
@define(qspace, theta, 1/(2*1im)*delta[:i,:j,:t0])
@define_basics(qspace)
println(sin(theta)*xi*beta)

expr = sin(theta) * xi * beta
term = expr.terms[1]

#∑(:l, xl) * ∑(:l, xl)
#∑(:m, xm) * ∑(:n, xn)
#exp(0*xi)

#dx_dt = d_dt(xi, log(Sum(:k, gamma[:k]  * yk))) # wrong gamma not deleted from before 
#expr = log(Sum("j", alpha * A()'^2 * yj + Sum("k", beta * gamma[:j] * gamma[:k] * A()^2 * xj * xk)))
#substitute(expr, A() --> xi)
