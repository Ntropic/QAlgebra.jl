# QAlgebra.jl

```@meta
CurrentModule = QAlgebra
```

**QAlgebra.jl** organises symbolic quantum modelling around three layers:

- **Operator spaces (`QSpaces`)** describe how quantum modes are arranged. A `QSpace`
  stitches together `OperatorSet`s as single subspaces or replicated ensembles and
  records whether an ensemble is treated as a discrete collection or a continuum.
- **Classical functions (`CFunctions`)** manage scalar parameters, time dependence,
  and ensemble-driven coefficient functions. They remain independent of the operator
  algebra so you can reuse them across different state-space layouts.
- **Quantum expressions (`QExpressions`)** build operator-valued formulas inside a
  given `QSpace`, combine them algebraically, and expose utilities for reordering,
  summing, and integrating with respect to ensemble indexes.

The examples below walk through defining a tensor-product space, instantiating the
associated parameters, and manipulating quantum expressions that depend on them.

---

## Example: Constructing a State Space

The state space is defined by describing each subsystem—either as a single mode or
an ensemble of replicated modes—declaring the symbolic operators, and listing the
classical parameter groups. Operator sets can be custom objects or picked from the
built-in Pauli, ladder, and PM collections. Ensembles can also be marked as
continuum approximations when their parameter ranges are sampled densely.

```@example qalgebra
using QAlgebra
subspace_def = SubSpaceDefinitions(
    h = QubitPM("beta"),
    i = Ensemble(3, 3, QubitPauli("sigma")),
    b = Ladder()
)
op_def = OperatorDefinitions("A(!i)", "B(U,H,i)")
param_def = ParameterDefinitions(
    "alpha",
    "beta(t)",
    "gamma(t)"
)
qspace = QSpace(subspace_def, op_def, param_def, max_t_ind=2)
qspace
```

```@setup qalgebra
ph, mh, zh = base_operators(qspace, "h")
xi, yi, zi = base_operators(qspace, "i")
xj, yj, zj = base_operators(qspace, "j")
xk, yk, zk = base_operators(qspace, "k")
xl, yl, zl = base_operators(qspace, "l")
b = base_operators(qspace, "b")
A = base_operators(qspace, "A", do_fun=true)
alpha, beta, gamma_time = base_operators(qspace, ["alpha", "beta", "gamma"], do_fun=true)
delta = QExprLookup([[:l], [:m]], [alpha, alpha])
t0, t1 = base_operators(qspace, :t)
```

```@example qalgebra
alpha, beta, gamma_time
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
it via the `parameter_groups` field, while `samples` and `distribution` remain available for later use when sampling
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
