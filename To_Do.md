# To Do List 03.09.2025

## High Priority 
- [x] Add time_indexes to QAtoms, then update: 
    - [x] base_operators
    - [x] QAtomProduct multiplication  & simplification
    - [x] commutes
    - [x] sort / isless 
    - [x] print => to print `t_{time_index}`

- [ ] Implement *Unitary Conjugation*
    - [x] Create *Conjugate (Krylov) Cycle* via Guassian elimination 
    - [x] Add new CFunction types `Cpower`, `CVector` and `CMatrix` 
        - [ ] Add modular Function definitions for substitutions, i.e. `sin(x)`, or symbols `
            - Add core parameter infos to `DataBase`
            - [x] `substitute` for c_function 
            - [x] grouping in helper
            - [x] update term_equal_indexes
            - [x] Functions to add abstracts and functions to ParamInfo 
                - [x] `contains_which_abstract`
            - [x] Evaluate
            - [x] Iterators
            - [x] print with index -> Distinguish CCustomTypes by whether they have CAbstracts or just a function hidden by a symbol -> the trivial CType 
            - [x] expand CCustomTypes 
            - [x] add cabstract to base_operators, 
            - [x] update CAlgebra
            - [x] add custom ctypes constructor function for every newly constructed one. i.e. if we define :sin, sin(x) should be a valid operation now. 
            - [x] Fix CFunction printing --> unnecessary coeff prints 
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
- [ ] Simplify custom CTypes -> expansion and pattern detection
- [ ] Additional **substitution** checks. 
    - [ ] add time_index check to substitutions, apply time_index transformations to operators at different time points!
- [ ] Extend `expand` to QExpressions, to expand cFunctions through QExpressions, or alternatively expand QExpressions directly 
## Low Priority 
 
- Additional **substitution** checks. 
    - [ ] are commutation properties of an optional *commutation_fun* fulfilled, for multiple parallel substitutions 
    - [ ] add time_index check to substitutions, apply time 