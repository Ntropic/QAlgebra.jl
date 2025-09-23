# To Do List 03.09.2025

## High Priority 
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
    - [ ] Finder of QInt types (which CFunctions are present) -> for this the CFunctions need to be reorder_full 
        -> and detection be permutation invariant between index names.

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
