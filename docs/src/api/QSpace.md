# Operator Spaces

```@meta
CurrentModule = QAlgebra.QSpaces
```

The `QSpaces` submodule defines the operator and variable structure of quantum systems. These spaces serve as the foundation for building symbolic quantum expressions in QAlgebra. It centers around the [`QSpace`](@ref) type, which holds the entire system structure: subspaces, variables, time dependence, and indexing.

## Core Types

The following types represent the building blocks of a composite quantum system:
- `OperatorSet` defines an algebra (e.g., Pauli, Ladder, PM) including how operators multiply, conjugate, and display.
- `SubSpace` represents one component of the full Hilbert space — each associated with an `OperatorSet`.
- `Ensemble` stores metadata for replicated subsystems with shared sampling behaviour.
- `Parameter` defines symbolic variables (e.g., $\alpha$, $\beta(t)$, $\gamma_i$) that appear in expressions.

```@docs
OperatorSet
SubSpace
Ensemble
Parameter
OperatorType
QSpace
SubSpaceDefinitions
ParameterDefinitions
```

## Parameter Groups

```@docs
ParameterGroupKind
ParameterGroup
resolve_param!
```

---

## Predefined Operator Sets

The package includes three standard operator sets for immediate use:

```@docs
Ladder
QubitPM
QubitPauli
```

They each provide appropriate `op_product`, `op_dag`, and rendering functions, and can be passed as values to the `QSpace` constructor.

---

## Sampler Helpers

```@docs
Main.QAlgebra.EnsembleSamples.DiscreteSamples
Main.QAlgebra.EnsembleSamples.ContinuousSamples
```

```@docs
Main.QAlgebra.Sampler.build_discrete_samples
Main.QAlgebra.Sampler.build_continuous_samples
```
