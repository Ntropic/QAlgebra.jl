# FFunctions

The `FFunctions` module provides the structs and operations needed to express functions of classical variables. 

```@meta
CurrentModule = QAlgebra.FFunctions
```

## Types / Structs

```@docs
FFunction
FAtom
FSum
FRational
```



## Interface & Utilities

```@docs
isnumeric
iszero
max_exponents
build_xpows
evaluate
to_string
stringer
sort_key
how_to_combine_Fs
```

## Arithmetic & Operations
The module overloads the following arithmetic operations for objects of the structs defined above: `+`, `-`, `*`, `/`, `^`, `==`, `ìnv`, `sort!`, `sort`, `copy`, `gcd`, `adjoint`, `conj`
You can also use `FFunction.simplify` to simplify your expressions, but this function is not explicitely exported.
