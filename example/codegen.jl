using LinearMPC
# Load inverted pendulum example and generate code 
example = LinearMPC.mpc_example("invpend")
mpc,range = example.mpc, example.range
LinearMPC.codegen(mpc;dir="codegen_impc")

# Compute explicit 
empc = LinearMPC.ExplicitMPC(mpc;range)
LinearMPC.codegen(empc;dir="codegen_empc")
