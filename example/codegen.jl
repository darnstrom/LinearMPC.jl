using LinearMPC
mpQP,P_theta,mpc = LinearMPC.mpc_examples("invpend",double_sided=true)
LinearMPC.codegen(mpc)

## BnB
mpQP,P_theta,mpc = LinearMPC.mpc_examples("invpend_contact",5,5,double_sided=true)
LinearMPC.codegen(mpc)
