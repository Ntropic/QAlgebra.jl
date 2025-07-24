# QAlgebra.jl

[![CI](https://github.com/Ntropic/QAlgebra.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Ntropic/QAlgebra.jl/actions/workflows/CI.yml)
[![Julia](https://img.shields.io/badge/Julia-1.10%2B-blue.svg)](https://julialang.org/)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://Ntropic.github.io/QAlgebra.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)


**QAlgebra.jl** provides symbolic tools for constructing and manipulating quantum operator expressions on structured composite systems of qubits, spins, and bosonic modes.

---

## Features

- **Flexible State Spaces**  
  Define composite systems with symbolic variables, indexed subsystems, and a variety of operator bases (Pauli, ladder, raising/lowering, …).

- **Symbolic Expressions**  
  Build, combine, and manipulate expressions of operators and parameters.

- **Index‑summations & Constraints**  
  Generate sums over indices, nest and flatten them, and enforce index‑distinction (e.g. \(i \neq j\)).

- **Operator Algebra**  
  Compute commutators, exponentials, logarithms, powers, and roots of operator expressions.

- **Differential Equations**  
  Define and manage operator differential equations for expectation values.

- **Simplification & Substitution**  
  Simplify operator products, reorder terms, and substitute symbols or sub‑expressions.

---

## Installation

```julia
import Pkg
Pkg.add("QAlgebra")
```


## Quick Start

```julia
using QAlgebra
```

### Constructing a State Space

Define a composite system with symbolic coupling constants and multiple operator bases:

```julia
qspace = StateSpace(
    "alpha", "beta(t)", "gamma_i", "delta_i";
    operators = ["A(!i)", "B(U,H,i)"],
    h = QubitPM(),           # single-qubit plus/minus basis
    i = (3, QubitPauli()),   # three Pauli‑qubit subsystems
    b = Ladder()             # single bosonic ladder mode
)
```

Access the building blocks:

```julia
# All variables, operators, and abstract operators
var_dict, op_dict, abstract_dict = base_operators(qspace)

# Individual scalar variables
alpha = base_operators(qspace, "alpha")
beta  = base_operators(qspace, "beta")

# Operators for subsystem “h” (returns a tuple [σ⁺, σ⁻, σᶻ])
ph, mh, zh = base_operators(qspace, "h", do_dict=false)

# Pauli operators on the i,j,k‑indexed qubits
xi, yi, zi, pi, mi = base_operators(qspace, "i", do_dict=false)
xj, yj, zj, pj, mj = base_operators(qspace, "j", do_dict=false)
xk, yk, zk, pk, mk = base_operators(qspace, "k", do_dict=false)

# Bosonic mode: ladder and number operators
b, n = base_operators(qspace, "b", do_dict=false)
```

You can also pull out parameterized or indexed abstract operators:

```julia
A   = base_operators(qspace, "A")       # abstract A(!i)
A1  = base_operators(qspace, "A_1")     # concrete index 1
As  = base_operators(qspace, "A", do_fun=true)  # function A(i)
A, A1, As(2)
```

### Building and Printing Expressions
Compose symbolic expressions just like you would in standard Julia:
```julia
expr = 2 * (alpha + beta) * im * xi + alpha * Dag(b) * xi * yi

# View as plain string or LaTeX
string(expr)
latex_string(expr)
```

### Summations & Index Constraints
```julia
# Nested sums over indices j and k
qsum = Sum("j",
    alpha * yi * yj +
    Sum("k", beta * alpha^2 * xi * xj * xk)
)

# Flatten nested sums
flatten(qsum)

# Enforce i ≠ j where appropriate
neq(qsum)
```

### Differential Equations
Define the time‑derivative of an operator expectation value:

```julia
diff_eq = d_dt(zi, alpha * expr + qsum)
```

### Operator Functions
Apply common functions to operator expressions:

```julia
# Commutator [·,·]
QCommutator(Sum("i", alpha * ph * xi * yi) + zj, zh)

# Exponential
exp(Sum("i", alpha * ph * xi * yi) + zj) + zh

# Logarithm
log(Sum("i", alpha * ph * xi * yi) + zj)

# Power and root
power(Sum("i", alpha * ph * xi * yi) + zj, 2) + zh

root(Sum("i", alpha * ph * xi * yi) + zj, 2) + zh
```

### Simplify, Substitute, and Reorder
We can simplify expressions with
```julia
# Simplify anti‑commuting Pauli products
simplify(xi*yi + yi*xi)
``` 

## Author 
- [Michael Schilling](https://github.com/Ntropic)