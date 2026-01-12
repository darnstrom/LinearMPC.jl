using LinearMPC
# Load inverted pendulum example and certify the iteration complexity 
mpc,range = LinearMPC.mpc_examples("invpend")
result = LinearMPC.certify(mpc;range)

using Plots
plot(result;parameters=[:x2,:u1p])
