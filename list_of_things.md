# Start
We are in the process of refactoring the QAlgebra package. In the past we were tracking QTerm's and CAtoms as fixed length Vectors of Vectors of Integers describing the operators and a fixed length Vector of a short list of possible parameters, respectively. 

The refactor now introduced QParticle and CParticle, and constructs the Atoms from variable numbers of these particles, that specify the positions and types. We still have a long way to go, many functions are not yet ported to the new way things work, please do not try to fix everything, I will always instruct you, which part we are changing, so that I can verify the changes properly at every stage. 


# Middle
I want to make sure, that isless and isequal work for both CAtoms and QAtoms, I think they already work for CParticles and QParticles. 



# Final
First I want you to not change anything, I want you to explain to me what you understood first. I will confirm when I think you are ready to proceed!

