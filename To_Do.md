# To Do List 03.09.2025

## High Priority 
- [x] Add time_indexes to QAtoms, then update: 
    - [x] base_operators
    - [x] QAtomProduct multiplication  & simplification
    - [x] commutes
    - [x] sort / isless 
    - [x] print => to print `t_{time_index}`
- [ ] Distinguish summation indexes from non-summation indexes
    - [x] `which_summation_acting` to track only summation indexes
    - [x] Refactored QSum
    - [x] Removed indexing
        - [ ] remove `indexes` from `QExpressions_simplify_composites
        - [ ] rewrite `QSum_modify` 
        - [ ] rewrite `QSum_repartition`
    - [x] add `parent` to QComposites and QExpr 
        - [x] make parents get tracked automatically when optimizing. 
    - [ ] repartition among sum indexes. 
    - [ ] Shift/repartition within summation indexes 
    - [ ] Implement a repartition variants only for sum indexes!

- [ ] Implement *Unitary Conjugation*
    - [x] Create *Conjugate (Krylov) Cycle* via Guassian elimination 
    - [x] Add new CFunction types `Cpower`, `CVector` and `CMatrix` 
        - [x] Add modular Function definitions for substitutions, i.e. `sin(x)`, or symbols `
    - [ ] Change `is_local` to `local_in_ensembles`!
    - [ ] Diagonalize CFunction Matrix to construct Operator transformations
    - [ ] Implement Statespace transformation
    - [ ] Remove (sum index of current operator) if present in H0 => `make_index_known` function, for example for `H0 = Sum(:i, gamma[:i]*zi)` 

- [ ] Add *EquationSet*
    - [ ] Construct *EquationSet* from *Loss Function*
        - [ ] QExpr to Function 
    
    - [ ] Implement Cumulants Properly for Terms outside of Operator Set 
        - [x] Add *Operator Moment Orders* functions (`max_moment_of_terms`)

    - [ ] Index EquationSet and Cumulants 

- [ ] Add *QInt* with integral solver also for time. Should only accept `t_index=0`
    - [ ] Numerical Solver
    - [ ] Sampler  

## Medium Priority 
- [x] is_t_var
- [ ] **Simplify custom CTypes** -> expansion and pattern detection
- [ ] Additional **substitution** checks. 
    - [ ] add time_index check to substitutions, apply time_index transformations to operators at different time points!
- [ ] Extend `expand` to QExpressions, to expand cFunctions through QExpressions, or alternatively expand QExpressions directly 
## Low Priority 
 
- Additional **substitution** checks. 
    - [ ] are commutation properties of an optional *commutation_fun* fulfilled, for multiple parallel substitutions 
    - [ ] add time_index check to substitutions, apply time 