# qAlgebra.jl

**qAlgebra.jl** provides symbolic tools for constructing and manipulating quantum operator expressions  
on structured composite systems of qubits, spins, and bosonic modes.

The core abstraction is the [`StateSpace`](@ref), which defines:
- the symbolic variables (e.g. coupling constants),
- the operator bases (e.g. Pauli, ladder, raising/lowering),
- and the indexed structure of composite systems.

---

## Example: Constructing a State Space

```@example qalgebra
using qAlgebra
qspace = StateSpace("alpha", "beta(t)", "gamma_i", "delta_i",
                    operators=["A(!i)", "B(U,H,i)"],
                    h=QubitPM(), i=(3, QubitPauli()), b=Ladder())
```

You can access variables and operators from this space in several ways.
By creating dictionaries of all variables, operators and abstract operators
```@example qalgebra ; output = false
var_dict, op_dict, abstract_dict = base_operators(qspace)
```
Individual variables 
```@example qalgebra ; output = false
alpha = base_operators(qspace, "alpha")
beta  = base_operators(qspace, "beta")
```

Individual subsystem bases (`do_dict=true` would return the entries as a dictionary):
```@example qalgebra ; output = false
ph, mh, zh = base_operators(qspace, "h", do_dict=false)
xi, yi, zi, pi, mi = base_operators(qspace, "i", do_dict=false)
b, n = base_operators(qspace, "b", do_dict=false)
```

Or specific abstract operators:
```@example qalgebra ; output = false
A = base_operators(qspace, "A")
```

---

## Building and Simplifying Expressions

You can build symbolic expressions using variables and operators:

```@example qalgebra
expr = 2 * (alpha + beta) * im * xi + alpha * Dag(b) * xi * yi
```

If you want to render it as math:

```@jldoctest
eq = latexify(simplify(alpha * (beta + beta^2)))
"$$" * String(eq) * "$$"
```

---

## Indexed and Nested Sums

```@jldoctest
qsum = Sum("j", alpha * yi * yj + Sum("k", beta * alpha^2 * xi * xj * xk))
flatten(qsum)
```

Force index distinction (i.e. $j \ne k$):

```@jldoctest
neq_sum(qsum)
```

---

## Differential Equations

```@jldoctest
diff_eq = d_dt(zi, alpha * expr + qsum)
```

---

## Operator Functions

```@jldoctest
d_dt(A, xi * Dag(A) * A + alpha * b * Dag(b))
```

```@jldoctest
simplify(qCommutator(Sum("i", alpha * ph * xi * yi) + zj, zh))
```

```@jldoctest
exp(Sum("i", alpha * ph * xi * yi) + zj) + zh
```

```@jldoctest
log(Sum("i", alpha * ph * xi * yi) + zj)
```

```@jldoctest
power(Sum("i", alpha * ph * xi * yi) + zj, 2) + zh
```

```@jldoctest
root(Sum("i", alpha * ph * xi * yi) + zj, 2) + zh
```

```@jldoctest
Sum("i", alpha * ph * xi * yi * A * xi + xi * yi) * Sum("j", zi)
```

---

## Author

- Michael Schilling
