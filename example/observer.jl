using LinearMPC

mpc,th_range = LinearMPC.mpc_examples("invpend",100)
move_block!(mpc,[1,1,5,10,10])
LinearMPC.set_state_observer!(mpc;Q=[1e-3,1,1e-3,1],R=[1,0.1])

N=1000;
x0 = zeros(4)
r = [zeros(2,20) repeat([10;0],1,N)]
dynamics = (x,u,d) -> mpc.model.F*x + mpc.model.G*u + [0 0; 0.05 0; 0 0;0 0.005]*randn(2)
get_measurement = (x,d) -> mpc.state_observer.C*x + [0.1;0.01].*randn(2)

sim = Simulation(dynamics,mpc;r,x0,get_measurement)

using Plots
plot(sim)
