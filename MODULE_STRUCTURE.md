# QAlgebra Module Structure

`QSpace` is the keystone object in QAlgebra. It binds **subspaces** (and their `OperatorSet`s), **abstract operators**, and **parameter collections** into a single runtime state that everything else reads from.

**Constructor inputs**

- `SubSpaceDefinitions` → supplies the subspace list, associated ensembles, and neutral operators.
- `OperatorDefinitions` → turns into `OperatorType` metadata scoped to the previously defined subspaces.
- `ParameterDefinitions` → produces `ParameterGroup`s, `Parameter`s, `ParameterInfo`, and initial payloads that `ParameterValues` stores.

Within that hub the parameter layer touches ensemble subspaces in two complementary ways:
- **Abstract linkage** uses `SubSpaceInfo`, `SubSpaceIndex`, and `ParameterInfo` to map symbolic outer/inner indices to expanded positions.
- **Sample linkage** relies on `ParameterValues` plus sampler objects to populate concrete values for those abstract slots.

Keeping the distinction between abstract indices (symbolic coordinates used during expression construction) and sample indices (actual draws or grid nodes) is crucial when moving between definition code and evaluation routines.

## Core Concept Distinctions
- **Parameter groups vs parameters** (`ParameterGroup` in `src/ParameterGroups.jl:32` vs `Parameter` in `src/QSpaceOps/QSpace_parameters.jl:14`): groups capture the canonical signature, dependencies, and payload for a family of values, while individual parameters realise that group across concrete time and ensemble index combinations.
- **Abstract ensemble indexes vs sample indexes** (`SubSpaceIndex` in `src/Helper.jl:1` vs runtime sample containers such as `ConcreteIndexes` and `ParameterValues`): abstract indexes describe symbolic positions (outer, inner, expanded) used during algebraic construction; sample indexes enumerate actual draws or grid points when evaluating ensemble-backed parameters.
- **Subspaces, ensemble subspaces, ensemble samples** (`SubSpace` in `src/QSpaceOps/QSpace_subspaces.jl:120`, `Ensemble` in `src/QSpaceOps/QSpace_subspaces.jl:15`, and `DiscreteSamples`/`ContinuousSamples` in `src/QSpaceOps/Sampler.jl:22`): a subspace is a named operator context, an ensemble subspace is a subspace tied to sampling metadata, and ensemble samples materialise the random or grid values backing distribution-defined parameter groups.
- **Operator sets vs operator types** (`OperatorSet` in `src/QSpace.jl:37` vs `OperatorType` in `src/QSpaceOps/QSpace_abstract.jl:5`): operator sets encode the algebra of concrete indexed operators, whereas operator types declare abstract families (with hermitian/unitary flags and subspace support) used when building expressions.
- **Parameter values vs abstract index parameters** (`ParameterValues` in `src/CFunctionsOps/ParameterValues.jl:18` vs `AbstractIndexParameters` in `src/CFunctionsOps/ParameterValues_abstract.jl:6`): the former stores realised data per group, time slot, and sample combination; the latter provides an index-aware snapshot that can be updated independently for symbolic manipulations.
- **Outer/inner/expanded indexes** (managed by `SubSpaceInfo` in `src/QSpaceOps/QSpace_subspaces.jl:269`): outer indexes select a subspace, inner indexes resolve into ensemble multiplicities, and expanded indexes linearise the pair for fast lookups across the algebra and parameter systems.
- **Distribution payloads vs ensemble functions** (handled via `ParameterGroupDistribution` and `ParameterGroupEnsembleFunction` in `src/ParameterGroups.jl:10`): a distribution payload provides marginal random variables per group, while an ensemble function maps sampled arguments (and optional time) to derived parameter values.

## Layered Overview
`src/QAlgebra.jl` knits together helper modules, sampling utilities, parameter infrastructure, subspace management, operator metadata, and the expression layer. The structs collaborate in layered fashion: subspace definitions establish the operator canvas, parameter definitions populate values over that canvas, operator definitions describe abstract families, and the `QSpace` aggregate binds everything for use by `CFunctions` and `QExpressions`.

### Parameter Layer
- `ParameterDefinitions` (see `src/QSpaceOps/QSpace_parameters.jl:57`) parses user declarations into `ParameterGroupDefinition` records, infers each group's kind, validates optional payloads, and records dependencies. Payloads can already contain literals, time functions, distributions, or ensemble functions, so initial values are present as soon as the definitions are parsed.
- `ParameterDefinitions2Parameters` (ibid. `:379`) expands group definitions into concrete `Parameter` instances, initialises `ParameterGroup` metadata, and assembles `ParameterInfo` (`src/CFunctions.jl:186`) with indexing maps for both algebraic and sampling contexts.
- During that expansion each group records `group_outer_indices` and `group_index_outer_subspaces`, aligning parameter arguments with the ensemble subspaces captured in `SubSpaceInfo`. This step ties abstract ensemble indices to the correct subspace blocks.
- `ParameterValues` (`src/CFunctionsOps/ParameterValues.jl:18`) allocates storage per group, tracks definition status, and encodes dependency-aware update waves. `WhereWhichParamGroup` (`src/ParameterGroups.jl:46`) partitions groups by kind so downstream code can quickly fetch “all distributions”, “all time functions”, and similar slices.
- `AbstractIndexParameters` (`src/CFunctionsOps/ParameterValues_abstract.jl:6`) wraps `ParameterInfo` with copy-on-write storage suitable for manipulations that reason about abstract ensemble indexes rather than realised samples.
- `ParameterDicts` and `ParameterIndexes` (`src/CFunctions.jl:140` and `:91`) provide lookup helpers from names to parameter indices and vice versa, enabling both runtime evaluation and expression-level introspection.

### Subspace & Ensemble Layer
- `SubSpaceDefinitions` (`src/QSpaceOps/QSpace_subspaces.jl:162`) collects named subspaces or ensembles, deriving label sets and neutral elements while guarding against symbol collisions.
- Each definition expands to a `SubSpace` (`src/QSpaceOps/QSpace_subspaces.jl:120`) that stores operator keys, index mappings, ensemble affiliation, and the attached `OperatorSet`. Ensemble subspaces carry an `Ensemble` (`src/QSpaceOps/QSpace_subspaces.jl:15`) with sampling strategy, payload tracking, and runtime sampler references (`AbstractEnsembleSample`).
- `SubSpaceInfo` (`src/QSpaceOps/QSpace_subspaces.jl:269`) consolidates the indexing topology, providing outer/inner/expanded conversions reused throughout parameter, operator, and expression code.
- `SubSpaceIndex` and `ConcreteIndexes` (`src/Helper.jl:1` and `:15`) offer lightweight wrappers for referencing subspace coordinates, bridging symbolic and concrete enumeration of ensemble positions.
- When an ensemble carries distribution-backed parameter groups, `assign_ensemble_samples!` (see `src/QSpace.jl:93`) pulls the abstract mapping from `ParameterInfo` and instantiates `DiscreteSamples`/`ContinuousSamples`, producing the sample-index realisations tied back into `ParameterValues`.

### Operator Layer
- `OperatorSet` (`src/QSpace.jl:37`) records the concrete algebra (products, daggers, formatting, magnitude bounds) for each family of indexed operators included in a subspace.
- `OperatorDefinitions` and `OperatorType` (`src/QSpaceOps/QSpace_abstract.jl:24` and `:5`) declare abstract operator families with optional hermitian/unitary/time properties and subspace support filters. `OperatorDefinitions2OperatorType` ties those declarations to the subspace layout.
- `OperatorTypeInfo` (`src/QSpaceOps/QSpace_abstract.jl:48`) computes commutation relationships (via `gen_commutes_function`) and stores operator metadata for downstream use in ordering, validation, and expression rewriting.
- `AbstractOperatorDicts` (`src/QSpace.jl:18`) index operator types by name, and are populated when `QSpace` is constructed.

### Sampling Runtime
- The `Sampler` module (`src/QSpaceOps/Sampler.jl`) defines `DiscreteSamples` and `ContinuousSamples` containers, utilities for CDF inversion and grid generation, and `build_discrete_samples` / `build_continuous_samples` factories.
- `assign_ensemble_samples!` (`src/QSpace.jl:93`) inspects each ensemble, verifies whether all distribution-backed groups have payloads, and materialises samplers that feed the parameter layer.

### QSpace Aggregator
- `QSpace` (`src/QSpace.jl:267`) is the hub that stores subspaces, `SubSpaceInfo`, operator types, `Parameter` vectors, `ParameterInfo`, `ParameterValues`, lookup dictionaries, and neutral `CFunction` instances. The constructor orchestrates the pipeline: subspaces first, then operator types, parameters, ensemble sampling, and finally per-operator/per-subspace dictionaries.
- Each ensemble captured in `QSpace.ensembles` receives a back-reference (`WeakRef`) so sampling state can refresh when parameter payloads change.

## Key Struct Connectivity
| Anchor struct | Connected structs | Purpose |
| --- | --- | --- |
| `ParameterDefinitions` | `ParameterGroupDefinition`, `ParameterGroup`, `Parameter` | Parses user declarations and expands groups into concrete parameters. |
| `ParameterGroup` | `ParameterValues`, `Ensemble` | Supplies runtime storage layout and registers ensemble affiliations for sampling. |
| `ParameterInfo` | `ParameterValues`, `ParameterIndexes`, `CFunctions` abstractions | Central index map for evaluation, substitution, and integral construction. |
| `SubSpaceDefinitions` | `SubSpace`, `SubSpaceInfo`, `OperatorSet`, `Ensemble` | Builds the operator canvas and attaches sampling metadata. |
| `OperatorDefinitions` | `OperatorType`, `OperatorTypeInfo`, `AbstractOperatorDicts` | Declares abstract operator families and their commutation properties. |
| `QSpace` | `ParameterValues`, `SubSpaceDicts`, `AbstractOperatorDicts`, `ReducedCumulantList` | Aggregates every layer for consumption by `QExpressions` and evaluation routines. |

## Construction Flow inside QSpace
1. `SubSpaceDefinitions` are validated and expanded into `SubSpace` vectors and `SubSpaceInfo`, recording ensemble-specific indexing (`src/QSpace.jl:280`).
2. `OperatorDefinitions2OperatorType` maps declared operator families onto the newly built subspaces, producing `OperatorType` and `OperatorTypeInfo`.
3. `ParameterDefinitions2Parameters` generates `Parameter`, `ParameterGroup`, and `ParameterInfo` objects, allocating `ParameterValues` storage and name dictionaries (`src/QSpaceOps/QSpace_parameters.jl:379`).
4. `assign_ensemble_samples!` materialises samplers for ensembles whose distribution-backed groups now have payloads (`src/QSpace.jl:93`).
5. The constructor finalises lookup dictionaries (`build_subspace_dicts`, `build_operator_dicts`), prepares neutral `CFunction` atoms, and links ensembles and parameter storage back to the owning `QSpace`.

## Supporting Index Utilities
- `SparsePermutationTools` (`src/Helper.jl:56`) stores sparse permutation matrices used when translating between different index orderings across subspaces and parameter groups.
- `ParameterIndexes` (`src/CFunctions.jl:91`) records which parameters introduce specific ensemble or time indexes, aiding inspectors like `where_acting` and `which_params_acting`.
- `ConcreteIndexes` and `SubSpaceIndex` bridge abstract expressions and concrete sampling loops, ensuring ensemble-aware algorithms receive both symbolic and numeric coordinates.
