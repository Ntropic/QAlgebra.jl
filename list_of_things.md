## List of things we need to do later on, to fix module again: 
wait. why are you changing stuff in `ParameterValues`? I asked for a change to the import, because it was broken. that should be a single line edit! and then I asked for a change to the `default_group_signature` function's way of printing indexes. it uses a way to create the indexes. I want the index generation function to have two modes, one for a standalone index, and then one for using it as an index. because an index that itself has an index i₁ (in latex: i_{1}) when used as an index becomes: ᵢ₁ (in latex: _{i_{1}}). just make that with a minimal change. i undid the last changes. then use that way of indexing on both the parameter groups and their arguments. when defining a parameter group via the string: "gamma_{i,j}(t, delta_i, delta_j", I want it to detect that  so that we get things like: γᵢ₁ⱼ₁(t,δᵢ₁,δⱼ₁) 


For the abstract mode of ParameterValues we need to determine the array sizes for the distributions. when we provide/update distribution values together with an AbstractIndex specifying the subspace/ensemble and slot in that ensemble, we have to store the provided distribution values into that index. for this we need to store similarly to the time index, the max_index_by_distribution_group, and resize the distribution arrays (and correspondigly the array sizes of the groups that depend on it. 

---

- I think in ParameterValues, we should update the number of time indexes, so that the array has enough space for the highest set time index. if that makes sense
- For abstractParameterValues, we again have to resize whenever a higher slot value is used when setting a value. 

---

- Update `QSpace_get_types.jl` for the new parameter style. it will no longer find params, but only parameter groups. 

- We need to get rid of `SubSpaceIndex` and `EnsembleIndex` and replace them with `AbstractIndex`

--- 

- I want us to split base_operators into multiple functions, by default we want to use the 

### Other 
- Concrete Index needs a rewrite, the expected lengths, need to be determined at construction, not based on any subspace parameters, but instead by checking for the maximum slot by ensemble value from some equations. 

### Removed functions: 
Fix call sited for:
- `map_by_tindex`
- `map_by_subspace`

To Do: 
symbol2formatted add 2nd version with do_latex::Bool, that only outputs one of the two. and 
