## Inverted Pendulum
using LinearMPC

opts = LinearMPC.MPCSettings()
mpQP,P_theta,mpc = LinearMPC.mpc_examples("invpend",settings=opts)

x0 = [0.0,0.0,0.0,0.0,-1]; 
rs = [zeros(2,20) repeat([20;0],1,80)];
xs,us,rs = LinearMPC.simulate(mpc,x0,1000;r=rs)

