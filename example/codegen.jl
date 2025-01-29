using LinearMPC
## Inverted pendulum  
mpQP,TH,mpc = LinearMPC.mpc_examples("invpend",50,10)
#LinearMPC.codegen(mpc)
empc = LinearMPC.ExplicitMPC(mpc;TH)
#ParametricDAQP.codegen(empc.solution)
#LinearMPC.codegen(empc)
fig = LinearMPC.plot_regions(empc)
#fig = LinearMPC.plot_feedback(empc;u_id = 1)
#pgfsave("invpend_regs5.tex", fig)

## Inverted pendulum with contact forces 
mpQP,P_theta,mpc = LinearMPC.mpc_examples("invpend_contact",5,5)
LinearMPC.codegen(mpc)
