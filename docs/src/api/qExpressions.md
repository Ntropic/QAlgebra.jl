# qExpressions

The `qExpressions` module provides the symbolic quantum expressions built on top of the operator and state spaces defined in `qSpace`.  
Expressions can represent single operator terms, sums of terms, indexed sums over subsystems, and full differential equations.

qAlgebra organizes expressions into:

- Core struct types
- Functions for constructing and simplifying expressions
- Pretty-printing functions for strings and LaTeX output
- Algebraic operations like commutators and symbolic sums

---

## 1. Expression Types

```@docs
qExpr
qTerm
qEQ
qSum
diff_qEQ
```

- [`qExpr`](@ref) is the abstract base type for all symbolic expressions.
- [`qTerm`](@ref) represents a single operator term with a coefficient.
- [`qEQ`](@ref) represents a sum of expressions.
- [`qSum`](@ref) represents a symbolic summation over indices.
- [`diff_qEQ`](@ref) represents a differential equation \( d/dt \langle \text{Op} \rangle \).

---

## 2. Basic Construction and Simplification

```@docs
term
simplify
base_operators
Sum
∑
flatten
neq
d_dt
is_numeric
```

- [`term`](@ref) constructs a single `qTerm` manually.
- [`simplify`](@ref) combines like terms and removes near-zero terms.
- [`base_operators`](@ref) extracts symbolic variables and operators from a [`StateSpace`](@ref).
- [`Sum`](@ref) and [`∑`](@ref) construct symbolic sums.
- [`flatten`](@ref) expands nested sums.
- [`neq`](@ref) enforces distinctness between summation indices.
- [`d_dt`](@ref) builds time derivatives (`diff_qEQ`).
- [`is_numeric`](@ref) checks if an expression simplifies to a pure number or neutral operator.

---

## 3. Printing and Output Formatting

```@docs
string
latex_string
```

- [`string`](@ref) returns a Unicode string representation of an expression (for console output).
- [`latex_string`](@ref) returns a LaTeX-compatible string representation for notebook/HTML rendering.

Pretty-printing happens automatically when using `display` or `println` in environments like Jupyter or Pluto.

---

## 4. Algebraic Operations
We overload the most common algebraic operations, such as `+`, `*`, `-`, and `/`, to work with `qTerm`, `qEQ`, `qSum` and `diff_qEQ` types. We forthermore add the following functions:
```@docs
Dag
Commutator
```

- [`Dag`](@ref) computes the Hermitian conjugate (dagger) of an expression.
- [`Commutator`](@ref) computes the commutator \([A, B] = AB - BA\).

You can also construct commutators using vector notation:

```julia
A = alpha * xi
B = beta * yi
comm = A + [A, B]  # appends the commutator [A, B] to expression A
```