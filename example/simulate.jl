## Inverted Pendulum
using Revise
using LinearMPC

opts = LinearMPC.MPCSettings()
opts.QP_double_sided=true
mpQP,TH,mpc = LinearMPC.mpc_examples("invpend",settings=opts)
empc = LinearMPC.ExplicitMPC(mpc;TH,build_tree=true)

x0 = [0.0,0.0,0.0,0.0,0]; 
rs = [zeros(2,20) repeat([20;0],1,730) repeat([19;0],1,80)];
dynamics = (x,u) -> mpc.F*x + mpc.G*u
@time xs,us,rs = LinearMPC.simulate(dynamics, empc,x0,1000;r=rs)

## Plotting
using Plots
plot(xs[1,:])
plot!(rs[1,:])
