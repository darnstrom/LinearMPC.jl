using LinearMPC
# Load inverted pendulum example and generate code 
mpc,range = LinearMPC.mpc_examples("invpend")
LinearMPC.codegen(mpc;dir="codegen_impc")
# Compute explicit 
empc = LinearMPC.ExplicitMPC(mpc;range)
LinearMPC.codegen(empc;dir="codegen_empc")
