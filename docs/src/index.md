# QAlgebra.jl

```@meta
CurrentModule = QAlgebra
```

**QAlgebra.jl** provides symbolic tools for constructing and manipulating quantum operator expressions on structured composite systems of qubits, spins, and bosonic modes.

The core abstraction is the `QSpace`, which defines:
- The symbolic variables (e.g. coupling constants),
- The operator bases (e.g. Pauli, ladder, raising/lowering),
- And the indexed structure of composite systems.

The examples below walk through defining a tensor-product space, the
associated parameters and abstract operators and some basic manipulations of quantum expressions built from them.

---

## Example: Constructing a State Space

The state space is defined by describing each subsystem—either as a single mode or
an ensemble of replicated modes—declaring the symbolic operators, and listing the
classical parameter groups. Operator sets can be custom objects or picked from the
built-in Pauli, ladder, and PM collections. Ensembles can also be marked as
continuum approximations when their parameter ranges are sampled densely.

```@example qalgebra
using QAlgebra
subspace_def = SubSpaceDefinitions( h=QubitPM("beta"), 
                                    i=Ensemble(3, 3, QubitPauli("sigma"), as_continuum=true), 
                                    b=Ladder(max_magnitude=4))
op_def = OperatorDefinitions("A(i,t)", "B(U,H,i)")
var_def = ParameterDefinitions( "alpha", 
                                "beta(t)" => t->t^2, 
                                "delta_i" => QUniform(0,1,10), 
                                "gamma_{i,j}(t,delta_i,delta_j)" => (t,gi,gj)->t*(gi-gj)) 
qspace = QSpace(subspace_def, op_def, var_def, max_t_ind=2)
```

```@setup qalgebra
alpha, beta, gamma, delta = base_operators(qspace, ["alpha", "beta", "gamma", "delta"], do_fun=true)
t0, t1 = base_operators(qspace, :t)
ph, mh, zh = base_operators(qspace, "h")
sigma = base_operators(qspace, "i", by_ensemble=true, do_fun=true)  # general constructor for all ensemble indexes
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
```

```@example qalgebra
alpha, beta[:t0], delta[:j], gamma[:i,:j,:t]
```
```@example qalgebra
A(), A(1)
```

```@example qalgebra
xi, yi, zi, ph, mh, zh
```

### Ensemble Configuration

When defining ensemble subsystems use [`Ensemble`](@ref) to collect the relevant metadata. The positional
arguments specify how many operator indexes and summation indexes are reserved alongside the `OperatorSet`.
Additional keyword arguments can be supplied, for instance to record the number of physically realised modes or to
pre-register sample points that span the ensemble parameter space:

```julia
Ensemble(2, 1, QubitPauli("sigma"); num_modes=4)
```

During [`QSpace`](@ref) construction each ensemble automatically records which parameter groups (e.g. `:gamma`) act on
it via the `param_groups` field, while `samples` and `distribution` remain available for later use when sampling
ensemble parameter spaces.


## Building, Printing and Modifying Expressions

You can build symbolic expressions using variables and operators:

```@example qalgebra
expr = alpha * xi + Dag(b) * yi
expr
```

You can get the `string` and `latex_string` expressions via 
```@example qalgebra
string(expr)
```

```@example qalgebra
latex_string(expr)
```

We provide a number of useful operations, such as sums and products:

```@example qalgebra
expr = ∑(:l, delta[:l] * yl) * ∑(:l, delta[:l] * yl * A())
expr
```
The automatic decollision step rewrites the second summation so each factor
acts on a unique index. You can still enforce explicit distinctness using
`neq`:

```@example qalgebra
neq(expr)
```

Integrals follow the same interface:

```@example qalgebra
qint = ∫("l", alpha * xl)
qint
```

You can define differential equations of operator expectation values
```@example qalgebra
qs = ∑(:l, alpha * xl)
diff_eq = d_dt(zi, alpha * qs + qs)
diff_eq
```
`d_dt` automatically normalises time indexes: if the equation only references a
single time slot `t_k`, it is internally substituted to `t_0` so subsequent
manipulations operate on a canonical representation.


### Operator Functions
We provide the following operator functions:
```@example qalgebra
exp(expr)
```

```@example qalgebra
log(expr + zi)
```

```@example qalgebra
power(expr, 2)
```

```@example qalgebra
root(expr + zi, 2)
```

---

## Author

- Michael Schilling
