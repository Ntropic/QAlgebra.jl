# Operator Spaces

The `statespace` submodule defines the operator and variable structure of quantum systems. These spaces serve as the foundation for building symbolic quantum expressions in QAlgebra. It centers around the code struct `StateSpace`, which holds the entire system structure: subspaces, variables, time dependence, and indexing.

## Core Types

The following additional types represent the building blocks of a composite quantum system:
- `OperatorSet` defines an algebra (e.g., Pauli, Ladder, PM) including how operators multiply, conjugate, and display.
- `SubSpace` represents one component of the full Hilbert space — each associated with an `OperatorSet`.
- `Parameter` defines symbolic variables (e.g., $\alpha$, $\beta(t)$, $\gamma_i$) that appear in expressions.

```@docs
OperatorSet
SubSpace
Parameter
OperatorType
StateSpace
```

---

## Predefined Operator Sets

The package includes three standard operator sets for immediate use:

```@docs
Ladder
QubitPM
QubitPauli
```

They each provide appropriate `op_product`, `op_dag`, and rendering functions, and can be passed as values to the `StateSpace` constructor.