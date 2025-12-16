# TODO – Subspace Refactor Follow Ups

## Update Ensemble Call Sites
- test/runtests_qalgebra.jl (multiple `Ensemble(…)` invocations)
- test/runtests_cfunctions.jl:81
- docs/src/index.md (example code snippets)
- prototyping/test_parametervalues.ipynb
- prototyping/test_qexpressions.ipynb
- prototyping/test_plotting.ipynb

## Replace Legacy SubSpace Fields
- QSpaceOps/QSpace_abstract.jl – stop using `ss_inner_ind`/`ss_outer_ind`.
- QExpressionsOps/QExpressions_algebra.jl – update loops over `ss.ss_inner_ind`.
- QExpressionsOps/QExpressions_print.jl – adjust printing routines removing `ss_inner_ind`.
- Any remaining references to `num_operator_indices` / `num_sum_indices` (search project-wide).

## Refactor SubSpaceInfo Consumers
- SubSpaceInfo now only exposes `subspaces`, ensemble mappings, and label metadata; update downstream code to stop expecting old attributes (`inner_labels_symbols`, `how_many_*`, etc.).
- QSpaceOps/QSpace_parameters.jl – relies on `inner_labels_symbols` and `how_many_*` fields.
- QExpressionsOps/QExpressions_reorder.jl & related helpers – expect counts from `SubSpaceInfo`.
- QExpressionsOps/QExpressions_welldefined.jl & ConstrainedIndexes.jl – depend on `how_many_*` arrays.

## Printing / Construction
- Integrate new label handling in CFunctions/QExpressions printing once particle-based indices are in place.
- Ensure new SubSpace/Ensemble displays propagate to documentation and tutorials.
- QSpace.jl: revisit the temporary neutral-element reconstruction once new data flows are in place.
