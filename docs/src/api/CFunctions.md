# CFunctions

```@meta
CurrentModule = QAlgebra.CFunctions
```

The `CFunctions` submodule handles symbolic **c-number** objects: scalars,
polynomials, and generic expressions that act as coefficients in quantum
formulas. It provides the infrastructure needed to define custom parameter
families, manipulate them algebraically, and integrate them with
`QExpressions`.

## Core Types

```@docs
CFunction
CAtom
CSum
CProd
CRational
CExp
CLog
CPower
CVector
CMatrix
CAbstract
CAbstractDefinition
CTypeDefinition
ParameterInfo
ParameterValues
CCustomType
```

---

## Constructors and Registration

Build custom coefficient types and register abstract symbols or numeric
wrappers that can be reused across expressions.

```@docs
define_cabstract
define_ctype
list_cabstracts
list_ctypes
define_cintegral
list_cintegrals
```

---

## Properties and Analysis

These helpers inspect coefficient expressions, return structural metadata,
and drive simplification passes shared with `QExpressions`.

```@docs
coeff
var_exponents
contains_non_simple_CFunction
evaluate
unique_first_terms
```

---

## Transformation Utilities

Standard algebraic routines for expanding and substituting inside coefficient
expressions mirror the operations available on the quantum side.

```@docs
expand
stringer
to_stringer
to_string
```

---

## Parameter Stores

```@docs
recompute_functions!
```

---

## Indexed Parameters

```@docs
CAtomIndexed
CIntegral
CIntegralDefinition
Indexed
param_index_tuples
which_params_acting
where_acting
separate_by_cond
ConcreteIndexes
```

---

## Convenience Macros

```@docs
Main.QAlgebra.@define
```
