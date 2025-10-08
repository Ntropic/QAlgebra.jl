# QExpressions

```@meta
CurrentModule = QAlgebra.QExpressions
```

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
QCompositeN
QMultiComposite
```

For the atomic operators we have 
```@docs
QTerm
QAbstract
```

For concrete composites we have
```@docs
QAtomProduct
QSum
QCompositeProduct
```
and specialised symbolic wrappers
```@docs
QExp
QLog
QPower
QRoot
QCommutator
```
```@docs
power
root
```


Finally we have quantum expressions that are constructed as linear combinations of composite operators and differential equations in time \( d/dt \langle \text{Op} \rangle \) which define the time derivative (of the expectation value) of a composite operator as a quantum expression
```@docs
QExpr
diffQEq
```

Construction normalises time indexes automatically: equations carrying a single
explicit index `t_k` are rewritten to `t_0`, ensuring downstream substitutions
and reorderings operate on a consistent baseline.

---

## QObj Construction 

Quantum expressions can be constructed via 
```@docs
base_operators
∑
d_dt
```
## QObj Modification
General modifications of quantum objectscan be performed with
```@docs
neq
flatsums
complexsums
reorder
reorder_time
reorder_full
substitute
```

### Ordering Helpers

```@docs
OrderedQAtomProduct
OrderbyOperator
```

### Substitution Helpers

```@docs
Substitution
Substitution_t
Substitution_index
```

### Analysis Helpers

```@docs
which_summations_acting
is_local
contains_time
contains_non_simple
contains_non_simple_QObj
contains_which_t_indexes
```

### Aggregators and Advanced Utilities

```@docs
AbstractQSum
QAtomOrdered
QAtomIndexed
Expectation
Sum2Int
decompose_sorted_blocks
which_abstracts
tree_iter_composite
modify_coeff_funs_tree_composite
∫
decollision_QSum
are_all_neq
NeqConstraint
```

## QObj Checks and Properties
The following scripts perform basic checks on the quantum expressions
```@docs
contains_abstract
which_ensemble_acting
are_indexes_defined
```

##  Printing and Output Formatting
```@docs
string
latex_string
```
```@docs
cumulant_string
```

Pretty-printing happens automatically when using `display` or `println` in environments like Jupyter or Pluto.

---

## 4. Algebraic Operations
We overload the most common algebraic operations, such as `+`, `*`, `-`, and `/`, to work with our introduced types. We furthermore add the following functions:
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
