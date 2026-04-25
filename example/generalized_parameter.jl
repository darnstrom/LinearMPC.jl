using LinearMPC,Plots
# Basic system
F,G = [1 0.1; 0 1], [0.005;0.1;;] # double integrator with Ts=0.1 

mpc = LinearMPC.MPC(F,G,C=[1.0 0])
mpc.settings.parameter_preview = true
set_objective!(mpc;Eu=[1.0;;])
add_constraint!(mpc,Au = [1;;], Ap =[1.0;;], ub=[1.0],ks=1:10)
add_constraint!(mpc,Au = [1;;], Ap =[-1.0;;], lb=[-1.0])

p = zeros(1,200) 
p[:,1:100] .= 1.0 

sim = Simulation(mpc;r=[1.0],p=p,N=200)

plot(sim)

mpqp = LinearMPC.mpc2mpqp(mpc)
