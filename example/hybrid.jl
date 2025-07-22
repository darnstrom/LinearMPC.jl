using LinearMPC
mpc,_ = LinearMPC.mpc_examples("satellite",20)
mpc.settings.reference_preview=false
x0, N = zeros(3), 20; 
rs = [zeros(1,5) 0.5*ones(1,N-5);
      zeros(2,N)];
dynamics = (x,u,d) -> mpc.model.F*x + mpc.model.G*u
sim = LinearMPC.Simulation(dynamics, mpc; x0,N, r=rs)
using Plots
plot(sim)
