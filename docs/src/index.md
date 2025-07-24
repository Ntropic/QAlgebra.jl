# QAlgebra.jl

**QAlgebra.jl** provides symbolic tools for constructing and manipulating quantum operator expressions  
on structured composite systems of qubits, spins, and bosonic modes.

The core abstraction is the [`StateSpace`](@ref), which defines:
- the symbolic variables (e.g. coupling constants),
- the operator bases (e.g. Pauli, ladder, raising/lowering),
- and the indexed structure of composite systems.

---

## Example: Constructing a State Space

```@example qalgebra
using QAlgebra
qspace = StateSpace("alpha", "beta(t)", "gamma_i", "delta_i",
                    operators=["A(!i)", "B(U,H,i)"],
                    h=QubitPM(), i=(3, QubitPauli()), b=Ladder())
```

You can access variables and operators from this space in several ways.
By creating dictionaries of all variables, operators and abstract operators
```@example qalgebra 
var_dict, op_dict, abstract_dict = base_operators(qspace)
nothing #hide
```
individual variables 
```@example qalgebra 
alpha = base_operators(qspace, "alpha")
beta  = base_operators(qspace, "beta")
nothing #hide
```
individual subsystem bases (`do_dict=true` would return the entries as a dictionary):
```@example qalgebra 
ph, mh, zh = base_operators(qspace, "h", do_dict=false)
xi, yi, zi, pi, mi = base_operators(qspace, "i", do_dict=false)
xj, yj, zj, pj, mj = base_operators(qspace, "j", do_dict=false)
xk, yk, zk, pk, mk = base_operators(qspace, "k", do_dict=false)
b, n = base_operators(qspace, "b", do_dict=false)
nothing #hide
```
or specific abstract operators (with or without indexes or as a function accepting an index):
```@example qalgebra
A = base_operators(qspace, "A")
A1 = base_operators(qspace, "A_1")
As = base_operators(qspace, "A", do_fun=true)
A, A1, As(2)
```

---

## Building, Printing and Modifying Expressions

You can build symbolic expressions using variables and operators:

```@example qalgebra
expr = 2 * (alpha + beta) * 1im * xi + alpha * Dag(b) * xi * yi
```

You can get the `string` and `latex_string` expressions via 
```@example qalgebra
string(expr)
latex_string(expr)
nothing #hide
```

We provide a number of useful operations, such as sums

```@example qalgebra
qsum = Sum("j", alpha * yi * yj + Sum("k", beta * alpha^2 * xi * xj * xk))
```
Nested sums can be flattened
```@example qalgebra
flatten(qsum)
```
and you can enforce index distinction (i.e. $j \ne k$):
```@example qalgebra
neq(qsum)
```

You can define differential equations of operator expectation values
```@example qalgebra
diff_eq = d_dt(zi, alpha * expr + qsum)
```


### Operator Functions
We provide the following operator functions, with as of now limited support
```@example qalgebra
QCommutator(Sum("i", alpha * ph * xi * yi) + zj, zh)
```

```@example qalgebra
exp(Sum("i", alpha * ph * xi * yi) + zj) + zh
```

```@example qalgebra
log(Sum("i", alpha * ph * xi * yi) + zj)
```

```@example qalgebra
power(Sum("i", alpha * ph * xi * yi) + zj, 2) + zh
```

```@example qalgebra
root(Sum("i", alpha * ph * xi * yi) + zj, 2) + zh
```

### Simplify, substitute and reorder Expressions
You can simplify expressions with 
```@example qalgebra
simplify(xi*yi+yi*xi)
```
---

## Author

- Michael Schilling
