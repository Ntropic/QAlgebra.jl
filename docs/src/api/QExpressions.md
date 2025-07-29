# QExpressions

The `QExpressions` module provides the symbolic quantum expressions built on top of the operator and state spaces defined in `statespace`.  
Expressions can represent single operator terms, sums of terms, indexed sums over subsystems, and full differential equations.

QAlgebra organizes expressions into:

- Core struct types
- Functions for constructing and simplifying expressions
- Pretty-printing functions for strings and LaTeX output
- Algebraic operations like commutators and symbolic sums

---

## Expression Types
We define the abstract types:
```@docs
QObj
QAtom
QComposite
```
Where `QObj` describes the main class of all quantum objects, with the subtypes `QAtom` and `QComposite` representing the simplest and more complex quantum objects, respectively.

For the atomic operators we have 
```@docs
QTerm
QAbstract
```

- [`QTerm`](@ref) specifies a concrete operator, with components on one or multiple of the subsystems defined in a state space.
- [`QAbstract`](@ref) is an abstract type for all quantum operators, which can later be substituted for a concrete operator, but is now only known by its properties, such as its name, types (Hermitian, unitary, and the subspaces on which it acts etc.).

For composite operators we have
```@docs
QAtomProduct
QCompositeProduct
QMultiComposite
QSum
```

Finally we have quantum expressions that are constructed as linear combinations of composite operators and differential equations in time \( d/dt \langle \text{Op} \rangle \) which define the time derivative (of the expectation value) of a composite operator as a quantum expression
```@docs
QExpr
diff_QEq
```

---

## QObj Construction 

Quantum expressions can be constructed via 
```@docs
term
base_operators
Sum
∑
d_dt
``` 
## QObj Modification
General modifications of quantum objectscan be performed with
```@docs
flatten
neq
QAlgebra.QExpressions.simplify
QAlgebra.QExpressions.substitute
QAlgebra.QExpressions.reorder!
```

## QObj Checks and Properties
The following scripts perform basic checks on the quantum expressions
```@docs
QAlgebra.QExpressions.is_numeric
QAlgebra.QExpressions.contains_abstract
QAlgebra.QExpressions.which_continuum_acting
QAlgebra.QExpressions.are_indexes_defined
```

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
We overload the most common algebraic operations, such as `+`, `*`, `-`, and `/`, to work with `QTerm`, `QAbstract`, `Expr`, `QSum` and `diff_QEq` types. We forthermore add the following functions:
```@docs
Dag
Commutator
```

You can also construct commutators using vector notation:

```julia
A = alpha * xi
B = beta * yi
comm = A + [A, B]  # appends the commutator [A, B] to expression A
```