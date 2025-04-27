using LinearMPC
# Load inverted pendulum example and certify the iteration complexity 
mpc,range = LinearMPC.mpc_examples("invpend")
result = LinearMPC.certify(mpc;range)
LinearMPC.plot(result,:x2,:u1p)
