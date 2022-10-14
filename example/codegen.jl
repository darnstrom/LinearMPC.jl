using LinearMPC
mpQP,P_theta,mpc = LinearMPC.mpc_examples("invpend",double_sided=true)
LinearMPC.codegen(mpc)
