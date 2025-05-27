# qAlgebra.jl

[![CI](https://github.com/Ntropic/qAlgebra.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Ntropic/qAlgebra.jl/actions/workflows/CI.yml)
[![Julia](https://img.shields.io/badge/Julia-1.10%2B-blue.svg)](https://julialang.org/)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://Ntropic.github.io/qAlgebra.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)


**qAlgebra.jl** provides a symbolic framework for operator algebra on composite quantum systems.

The operators are defined via a **StateSpace**, on which both the *variables* and *operators* are defined.  
On these spaces, one can construct and manipulate *quantum expressions* and construct *differential equations* for expectatiuon values of (composite) operators.

Key features include:
- Structured representation of operators across subsystems
- Algebraic manipulation and simplification of quantum expressions
- Support for sums over bath indexes with distinct or overlapping indices
- Pretty-printing and LaTeX rendering of expressions

## Installation
You can install the package directly from GitHub:
```julia
import Pkg
Pkg.add(url="https://github.com/Ntropic/qAlgebra.jl")
```

## Quick Usage Example
```julia
using qAlgebra

# Generate the State Space consisting of a spin h using the raising and lowering basis (QubitPM), a spin bath in the Pauli Bases (QubitPauli) and a bosonic mode (Ladder).
qs = StateSpace("alpha", "beta(t)", "gamma_i", "delta_i", h=QubitPM(), i=(3, QubitPauli()), b=Ladder())

# We can return a dictionary of the operators in the system
var_dict, op_dict = base_operators(qs)

# Or individually for each subsystem (or subsystem index, i.e. i,j, and k) and variable
xi, yi, zi = base_operators("i", qs)
xj, yj, zj = base_operators("j", qs)
xk, yk, zk = base_operators("k", qs)
ph, mh, zh = base_operators("h", qs)
b = base_operators("b", qs)
I = base_operators("I", qs)     # Identity operator
alpha, beta, gamma_i, gamma_j, gamma_k, delta_i, delta_j, delta_k = base_operators("vars", qs)

# A simple expression can then be constructed from the operators 
expr = 2 * alpha * im * xi + alpha * Dag(b) * xi * yi

# It can be simplified via
simplified_expr = simplify(expr)

# Sums can be constructed via 
qsum = Sum("j", alpha*yi*yj+Sum("k", beta*alpha^2*xi*xj*xk))

# and the nested sum flattened via 
flat_qsum = flatten(qsum)

# The sum still covers all combinations of indexes j,k
# We can transform it into a neq sum, in which the indexes j and k are distinct. the following function then expands into all possible cases
neq_sum = neq_sum(qsum) # this also flattens the sum

# A differential equation of expectation values can be constructed via
diff_eq = d_dt(zi, alpha*expr+sum) # this simplifies, flattens and neq's the quantum equation
```

## Author 
- [Michael Schilling](https://github.com/Ntropic)