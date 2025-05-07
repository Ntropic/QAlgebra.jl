# Operator Spaces

The `qSpace` submodule defines the operator and variable structure of quantum systems. These spaces serve as the foundation for building symbolic quantum expressions in qAlgebra. It centers around the code struct `StateSpace`, which holds the entire system structure: subspaces, variables, time dependence, and indexing.

## Core Types

The following additional types represent the building blocks of a composite quantum system:
- `OperatorSet` defines an algebra (e.g., Pauli, Ladder, PM) including how operators multiply, conjugate, and display.
- `SubSpace` represents one component of the full Hilbert space — each associated with an `OperatorSet`.
- `Parameter` defines symbolic variables (e.g., $\alpha$, $\beta(t)$, $\gamma_i$) that appear in expressions.

```@docs
OperatorSet
SubSpace
Parameter
StateSpace
```

---

## Predefined Operator Sets

The package includes three standard operator sets for immediate use:

- `QubitPauli` defines a qubit in the standard Pauli basis ($\sigma_x$, $\sigma_y$, $\sigma_z$)
- `QubitPM` defines a qubit with aising/lowering basis (`σ⁺`, `σ⁻`, `σᶻ`)
- `Ladder` defines a bosonic mode via ladder operators (`b`, `b†`)

These are loaded from:

```julia
include("OperatorSets/Qubit_Pauli.jl")
include("OperatorSets/Qubit_PM.jl")
include("OperatorSets/Ladder.jl")
```

They each provide appropriate `op_product`, `op_dag`, and rendering functions, and can be passed as values to the `StateSpace` constructor.