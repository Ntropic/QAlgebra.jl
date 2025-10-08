# Other Functions

```@meta
CurrentModule = QAlgebra
```

The following are currently exposed functions and constants 

```@docs
symbol2formatted
str2sub
str2sup
QAlgebra.StringUtils.superscript_indexes
QAlgebra.StringUtils.subscript_indexes
QAlgebra.StringUtils.brace
QAlgebra.StringUtils.braket
```

---

## Sampling Utilities

```@docs
QDistribution
QNormal
QUniform
QEnsembleFunction
pdf
QAlgebra.Sampler.pdf2cdf
QAlgebra.Sampler.cdf2inverse
QInterpolator
nodes
basis_values
basis_values!
eval_interpolation
QIntegrator
integrate_node_funs
eval_integration
```

---

## Indexing Utilities

```@docs
BinomialCache
EnsembleRankWorkspace
MultiEnsembleWorkspace
index_number_for_ensemble
combined_index_for_ensembles!
multi_ensemble_iterator
```

---

## Default Preferences

```@docs
get_default
set_flip_if_first_term_negative
set_do_braced
set_expand_cumulants
```

---
