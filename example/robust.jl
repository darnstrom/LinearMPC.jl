using LinearMPC,Plots

# Basic system
F,G = [1 0.1; 0 1], [0.005;0.1;;] # double integrator with Ts=0.1 
true_dynamics = (x,u,d) -> F*x + G*u+0.01*(rand(2).-0.5);

# Nominal MPC
mpc_nominal = LinearMPC.MPC(F,G;Ts=0.1,Np=25,C=[1 0;])
set_bounds!(mpc_nominal;umin=[-0.2],umax=[0.2],ymin=[-0.5],ymax=[0.5]) 

sim_nominal = Simulation(true_dynamics,mpc_nominal;r=[0.5])
hline([0.5],label="Constraint bound", linestyle=:dash)
plot!(sim_nominal.ys[1,:],xlabel="Time step", ylabel="Position [m]",label="Nominal MPC")

# MPC with constraint tightening  
mpc_robust = LinearMPC.MPC(F,G;Ts=0.1,Np=25,C=[1 0;])
set_prestabilizing_feedback!(mpc_robust)
set_disturbance!(mpc_robust,[-0.005;-0.005],[0.005;0.005])
set_bounds!(mpc_robust;umin=[-0.2],umax=[0.2],ymin=[-0.5],ymax=[0.5]) 

sim_robust = Simulation(true_dynamics,mpc_robust;r=[0.5])
plot!(sim_robust.ys[1,:], label="Robust MPC")

# Rerun with the worst possible disturbance
worst_case_dynamics = (x,u,d) -> F*x + G*u+0.005*ones(2);
sim_robust_wc = Simulation(worst_case_dynamics,mpc_robust;r=[0.5])
hline([0.5],label="Constraint bound", linestyle=:dash)
plot!(sim_nominal_wc.ys[1,:],xlabel="Time step", ylabel="Position [m]",label="Nominal MPC")
plot!(sim_robust_wc.ys[1,:], label="Robust MPC")
