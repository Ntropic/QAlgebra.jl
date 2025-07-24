# Defaults and Precision Settings

QAlgebra allows configuration of the number formatting used for coefficients, including float and exponential precision.  
These defaults are persistent across sessions using `Preferences.jl`.

## API

```@docs
get_default
set_flip_if_first_term_negative
set_do_braced
```