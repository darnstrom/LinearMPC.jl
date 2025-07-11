using LinearMPC

# Get inverted pendulum MPC 
mpc,_= LinearMPC.mpc_examples("invpend")
# Run simulation 
x0 = [0.0,0.0,0.0,0.0]; 
rs = [zeros(2,20) repeat([10;0],1,780) repeat([9;0],1,10)];
dynamics = (x,u,d) -> mpc.model.F*x + mpc.model.G*u
sim = LinearMPC.Simulation(dynamics, mpc; x0, N = 1000, r=rs)

## Plot the trajectory 
using Plots
plot(sim)
