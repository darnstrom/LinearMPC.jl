using LinearMPC
# Load inverted pendulum example and certify the iteration complexity
example = LinearMPC.mpc_example("invpend")
mpc,range = example.mpc, example.range
result = LinearMPC.certify(mpc;range)

using Plots
plot(result;parameters=[:x2,:u1p])
