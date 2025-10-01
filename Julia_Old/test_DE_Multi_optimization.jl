# Generate sets of differential equations for single qubits and multiple qubits in a cavity up to n'th operator order using the hierarchical equations of motion approach
# Finally use cumulants to make the problem solvable
include("operator_terms.jl")
include("combinatorics.jl")
include("cumulants.jl")
include("diff_Eq.jl")
include("preprocessing.jl")
include("preprocessing2.jl")
include("indexing.jl")
include("initial_states.jl")
include("pulses.jl")
include("sampled_integrals.jl")
include("diff_Eq_solver.jl")
include("plotting.jl")
include("saving.jl")
println("Included all files")
using Random
rng = MersenneTwister(1234)

#using ForwardDiff
#using Optimization, Ipopt, OptimizationOptimJL 
# Pulse parameters ---------------------------------------------------------------------------
max_A::Float64 = 0.0
T::Float64 = 10
n_t::Int = 250
how_many::Int = 1
pulse_param::Pulse_Param_Struct = Random_Pulse_Param_Sin(max_A, T, rng=rng, how_many=how_many)
sin_pulses = sin_pulse_generator(pulse_param)
val_dict::Dict = Dict("\\beta" => sin_pulses, "\\Gamma" => 2 * pi * 0.0, "\\gamma" => 2 * pi * 0.0, "\\kappa" => 2 * pi * 0.0) # Gamma is Purcell decay rate, gamma is the spin dephasing rate, kappa is the cavity decay rate

# --- Samples and Weights --------------------------------------------------------------------
var_strs = ["g", "\\Delta"]
sample_order = [10, 10] #[15,15]
g0 = 2 * pi * 150 / 10^6  # 150Hz
delta_g = g0 * 0.20
delta0 = 0.0
delta_delta = 2 * pi * 0.1 # 100kHz 
sigmas = 3.0
N = 10^6
max_order = 3

spin_state = "single_non_excitation"

circle = false
gaussian = (x, mu, sigma) -> exp(-(x - mu)^2 / sigma^2 / 2) * 1 / sqrt(2 * pi * sigma^2)
prob_fun_vec = [x -> gaussian(x, g0, delta_g), x -> gaussian(x, delta0, delta_delta)]
range_vec = [[g0 - sigmas * delta_g, g0 + sigmas * delta_g], [-delta_delta * sigmas, delta_delta * sigmas]]
samples = prepare_samples_and_weights(prob_fun_vec, var_strs, ["g"], N, range_vec, sample_order; prepare_analysis=1, max_comb_order=-1, on_borders=false, circle=circle, reltol=1e-15, abstol=1e-15, extra_points=0)

vals2coeffs = samples.vals2coeffs
coeffs2vals = samples.coeffs2vals
fprintln("Fit Matrix Deviation: ", maximum(vals2coeffs * coeffs2vals - Diagonal(ones(size(vals2coeffs, 1)))), " at n=", length(samples.locations), " samples.")

spin_param = Dict("N" => 10^6, "A" => 1.0, "phi" => even_odd_vector(samples, :even, :complex, rng=MersenneTwister(1234)))
cavity_param = Dict("alpha" => 0.0, "Theta" => 0.0)
# --- Generate all equations and terms --------------------------------------------------------

system = prepare_indexed_eqs_from_samples(max_order, samples, do_lower_order_cumulants=true, less_spins=1, printing=true)

# --- Initial Conditions ---------------------------------------------------------------------
param_generator, constant_value_indexes = value_generator(pulse_param, val_dict, system) # Pulse Values
initial_conditions = operators2initialstates(system, samples, spin_state, "coherent", spin_param, cavity_param)
println("Done preparing initial conditions")

# --- Solving the System ---------------------------------------------------------------------
#all_eqs_indexed_reduced = remove_zero_terms(all_eqs_indexed, param_generator, constant_value_indexes; how_many_eps=5)
diff_problem = DifferentialProblem(param_generator, system)  # all_eqs_indexed_reduced
solution = evolve_system(initial_conditions, diff_problem, T, solvetype=:nonstiff, n_t=n_t, reltol=1e-15, abstol=1e-15, save_cummulants=0, save_lower_cummulants=1, threaded=true, print_dt=T / 5)
println("Done")

# Save the Solution
# create a filename 
filename = "free_evolution_$(spin_state)_T$(T)_order_$(max_order)_$(less_spins)_sin_pulse_$(how_many)_pulses"
# save the solution
save_solution(solution, filename)