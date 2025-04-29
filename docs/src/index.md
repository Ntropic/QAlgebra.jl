# qAlgebra.jl

**qAlgebra.jl** provides symbolic tools for building and manipulating quantum operator expressions  
on structured composite systems of qubits, spins, or bosonic modes.

The core idea is to first define a [`StateSpace`](@ref), which encodes the variables, subsystems, and operator bases.  
From this space, symbolic expressions can be constructed and manipulated. Our expression types are

- [`qTerm`](@ref): individual operator terms
- [`qEQ`](@ref): additive symbolic expressions
- [`qSum`](@ref): symbolic sums over indices
- [`diff_qEQ`](@ref): symbolic time evolution equations

## Key Features

- Composite operator algebra over structured subsystems
- Indexed sums with support for distinct (`neq`) indices
- Expression simplification, and expansion
- LaTeX output

## Example

```julia
using qAlgebra

# Generate the State Space consisting of a spin h using the raising and lowering basis (QubitPM), a spin bath in the Pauli Bases (QubitPauli) and a bosonic mode (Ladder).
qs = StateSpace("alpha", "beta(t)", "gamma_i", "delta_i", h=QubitPM(), i=(3, QubitPauli()), b=Ladder())

# We can return a dictionary of the operators in the system
var_dict, op_dict = base_operators(qs)

# Or individually for each subsystem (or subsystem index, i.e. i,j, and k) and variable
xi, yi, zi = base_operators("i", qs)
xj, yj, zj = base_operators("j", qs)
xk, yk, zk = base_operators("k", qs)
ph, mh, zh = base_operators("h", qs)
b = base_operators("b", qs)
I = base_operators("I", qs)     # Identity operator
alpha, beta, gamma_i, gamma_j, gamma_k, delta_i, delta_j, delta_k = base_operators("vars", qs)

# A simple expression can then be constructed from the operators 
expr = 2 * alpha * im * xi + alpha * Dag(b) * xi * yi

# It can be simplified via
simplified_expr = simplify(expr)

# Sums can be constructed via 
qsum = Sum("j", alpha*yi*yj+Sum("k", beta*alpha^2*xi*xj*xk))

# and the nested sum flattened via 
flat_qsum = flatten(qsum)

# The sum still covers all combinations of indexes j,k
# We can transform it into a neq sum, in which the indexes j and k are distinct. the following function then expands into all possible cases
neq_sum = neq_sum(qsum) # this also flattens the sum

# A differential equation of expectation values can be constructed via
diff_eq = d_dt(zi, alpha*expr+sum) # this simplifies, flattens and neq's the qEQ
```

To exemplify the LaTeX rendering, we can display an expression via
```julia
display(diff_eq)
```
Which will return the LaTeX expression 
$$
\frac{\text{d} \phantom{t}}{\text{d}t}\!\,\braket{\hat{\sigma}_{z}^{(i)}} = \sum_{j}^{\neq} \left(\alpha\,\braket{\hat{\sigma}_{y}^{(i)}\hat{\sigma}_{y}^{(j)}}+2\alpha^2\beta(t)\,\left(\,\braket{\hat{\sigma}_{x}^{(i)}}+\,\braket{\hat{\sigma}_{x}^{(j)}}\right)\right)+\sum_{(j,k)}^{\neq} \alpha^2\beta(t)\,\braket{\hat{\sigma}_{x}^{(i)}\hat{\sigma}_{x}^{(j)}\hat{\sigma}_{x}^{(k)}}+\sum_{k}^{\neq} \alpha^2\beta(t)\,\braket{\hat{\sigma}_{x}^{(k)}}+\alpha\,+i\alpha^2\,\left(2\,\braket{\hat{\sigma}_{x}^{(i)}}+\,\braket{\hat{\sigma}_{z}^{(i)}\hat{b}^\dagger}\right)
$$

