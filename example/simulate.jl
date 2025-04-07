## Inverted Pendulum
using LinearMPC

opts = LinearMPC.MPCSettings()
opts.QP_double_sided=true
mpc,range = LinearMPC.mpc_examples("invpend",settings=opts)
mpc.weights.Q*=10; 
empc = LinearMPC.ExplicitMPC(mpc;range,build_tree=true)

x0 = [0.0,0.0,0.0,0.0]; 
rs = [zeros(2,20) repeat([10;0],1,730) repeat([9;0],1,80)];
dynamics = (x,u) -> mpc.F*x + mpc.G*u
@time xs,us,rs = LinearMPC.simulate(dynamics, empc,x0,1000;r=rs)

## Plotting
using Plots
plot(xs[1,:])
plot!(rs[1,:])
