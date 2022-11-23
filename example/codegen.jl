using LinearMPC
## Inverted pendulum  
mpQP,P_theta,mpc = LinearMPC.mpc_examples("invpend")
mpc.settings.QP_double_sided = true
LinearMPC.codegen(mpc)

## Inverted pendulum with contact forces 
mpQP,P_theta,mpc = LinearMPC.mpc_examples("invpend_contact",5,5)
LinearMPC.codegen(mpc)
