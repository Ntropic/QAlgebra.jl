# qExpressions

The `qExpressions` module provides the symbolic quantum expressions built on top of the operator and state spaces defined in `statespace`.  
Expressions can represent single operator terms, sums of terms, indexed sums over subsystems, and full differential equations.

qAlgebra organizes expressions into:

- Core struct types
- Functions for constructing and simplifying expressions
- Pretty-printing functions for strings and LaTeX output
- Algebraic operations like commutators and symbolic sums

---

## Expression Types
We define the abstract types:
```@docs
qObj
qAtom
qComposite
```
Where `qObj` describes the main class of all quantum objects, with the subtypes `qAtom` and `qComposite` representing the simplest and more complex quantum objects, respectively.

For the atomic operators we have 
```@docs
qTerm
qAbstract
```

- [`qTerm`](@ref) specifies a concrete operator, with components on one or multiple of the subsystems defined in a state space.
- [`qAbstract`](@ref) is an abstract type for all quantum operators, which can later be substituted for a concrete operator, but is now only known by its properties, such as its name, types (Hermitian, unitary, and the subspaces on which it acts etc.).

For composite operators we have
```@docs
qAtomProduct
qCompositeProduct
qMultiComposite
qSum
```

Finally we have quantum expressions that are constructed as linear combinations of composite operators and differential equations in time \( d/dt \langle \text{Op} \rangle \) which define the time derivative (of the expectation value) of a composite operator as a quantum expression
```@docs
qExpr
diff_qEQ
```

---

## Basic Construction and Simplification

```@docs
term
qAlgebra.qExpressions.simplify
base_operators
Sum
∑
flatten
neq
d_dt
qAlgebra.qExpressions.is_numeric
```

---

##  Printing and Output Formatting

```@docs
string
latex_string
```

- [`string`](@ref) returns a Unicode string representation of an expression (for console output).
- [`latex_string`](@ref) returns a LaTeX-compatible string representation for notebook/HTML rendering.

Pretty-printing happens automatically when using `display` or `println` in environments like Jupyter or Pluto.

---

## 4. Algebraic Operations
We overload the most common algebraic operations, such as `+`, `*`, `-`, and `/`, to work with `qTerm`, `qAbstract`, `qExpr`, `qSum` and `diff_qEQ` types. We forthermore add the following functions:
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