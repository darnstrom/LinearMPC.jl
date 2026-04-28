using LinearMPC

# Load an example and run its first predefined scenario
example = LinearMPC.mpc_example("invpend")
sim = LinearMPC.Simulation(example, 1)

## Plot the trajectory 
using Plots
plot(sim)
