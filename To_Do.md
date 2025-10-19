### Aspects to change:
#### BS Removal 

#### High Priority
 - Add a check, if all parameter group payloads have been resolved before creating `diff_QEq` sets.
 - Implement `CIntegral` construction of `QIntegral` 
 - Add a reference to `CIntegral`s in `ParameterValues` with a `cintegrals_by_group::Vector{BitVector}` vector field specifying which parameter groups act on eahc of the `CIntegral`s. This allows us to automatically update affected `CIntegral`s on calls of `update_t!` or `resolve_param!` if one or more of the parameter groups, that they depend upon was changed. 
    - In order to reduce the number of such occurances, and assuming that the most common update is `update_t!`, we attempt to separate parameters of time in the automatic construction of the `CIntegral` set for `diff_QEq`'s.

#### Medium Priority
 - Add a Integral Product simplification option. the CIntegral, should maybe be transformed, so that the integrand is given by both an `EnsembleIndex` and the index of the `ParameterGroup`. 

#### Low Priority
 - Change integrands/summands in `QSum` and `CIntegral` expressions, from `SubSpaceIndex`es to `EnsembleIndex`es. This should remove the need to determine the ensemble indexes from `SubSpaceIndex`es. 
 - Add new `base_operator` alternative, qspace `getindex!` method to get the operators that replaces the old `base_operators`. 
 - for `diff_QEq` add a `is_hermitian_lhs::Bool` field, this allows us, to perform a consistency check --> for hermitian operators the expectation values must remain real!
 - Add new `IFunction` module, for Index functions, supporting only the following notation a*abstractindex+b where a and b are integers. 
 - Make string construction optional for Operator Sets, with a default variant.  

### For CIntegral implementation

I am considering a refactor of `CIntegralDefinition`. Please don't change anything, for mow we discuss the ideas. 
- at construction the integral interpolator should not be built! It needs to be built, deliberately later on. This is to ensure we only calculate integral weights once all parameter groups payloads have been deployed. 
- I would like to not enly specify the `SubSpaceIndex`es, but more importantly the distribution parameters of the `Ensemble`, belonging to the outer index of `SubSpaceIndex`es. In the end we anyways integrate over those parameter-abstract subspace index pairs. 
- By separating each integrand like this, we can also check if there are factorizations of the integral into other integrals, leading to dimensionality reductions of the integrals and hence significant performance gains. 
- Ideally such a check is performed at CIntegral definition. We already have a function that attempts to separate CFunctions into sums and products of CFunctions, so that some parameters are only in one or the other term. I guess this could be expanded into a more general framework (this function is placed in test_qexpessions.ipynb and has as of yet not been integrated into the main package).   

#### Added: 
we will not use ensure integrator or anything like this. I will take care of that later. but for clarity, I will add a reference to the cintegraldefinitions to PArameterValues, and trigger automatic recomputing of the integrals when payloads or times are changed. 
but the ensuring that the values are calculated is not necessary right now. we can take cre of it later. 

I think we should for now focus on writing the refactor of how we create `CIntegral`, so that it adds the parameters of the ensemble subspaces. what do we need to pass to the consructor for this? do we pass the subspaces vector to the constructor for this? we should be able to import subspaces, or does it itself depend on QSpace_parameters.jl?  

 #### Prompt
 Review the `QAlgebra.jl` julia package . We have added a `ParameterValues` struct in order to facilitate the evaluation of `CFunction` expressions. For `CIntegral`s, I wish to evaluate the integrals on the sample nodes, at which the acting Ensembles are defined. We have a `QInterpolator` and `QIntegrator` module that will help with that. 
 - It is important to note how the integration is performed. The integrands are specified as `SubSpaceIndex`es. Yet this is only a proxy, for the integration we integrate over the parameter space associated with the `Ensemble`s specified by those indexes. The integration volume is the tensor product space of the parameter spaces of each of the integrand `SubSpaceIndex`es, i.e. for indexes `i` and `j` belonging to the same `Ensemble` with parameters `a` and `b`, we integrate over the product space of $d a_i \times d b_i \times d a_j \times d b_j$. We integrate from the `minimum` to the `maximum` values of the `QDistribution` given as `payload` for the `ParameterGroup` of each of the parameters. The integration is performed in parallel for all basis functions of the sample node based interpolator, but requies that for each parameter value combination requested by the internal `quadgk` call, we map the parameters to a `ParameterValues{AbstractIndexMode}` (which we can do via `resolve_distribution_values!` using the `SubSpaceIndex`es of the integrand) and evaluate the `CFunction` for these values via`abstract_evaluate`. 
 
 Please develop a succinct plan to construct the `QInterpolator` for a `CIntegral` expression, to determine the weights for each sample node combination, so that the integral can be transformed into a weighted sum of the nodes. 