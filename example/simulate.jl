## Normal invpend
using LinearMPC

opts = LinearMPC.MPCSettings()
mpQP,P_theta,mpc = LinearMPC.mpc_examples("invpend",settings=opts)

x0 = [10.0,0.0,0.0,0.0,1]; 
rs = [zeros(2,50) repeat([10;0],1,250)];
xs,us,rs = LinearMPC.simulate(mpc,x0,3000;r=rs)

## Contact invpend 

mpQP,P_theta,mpc = LinearMPC.mpc_examples("invpend_contact",10,10)

x0 = [0.0, 0.0, 1.0, 0.0]; 
#rs = [zeros(2,50) repeat([10;0],1,250)];
#xs,us,rs = LinearMPC.simulate(mpc,x0,3000;r=rs)
xs,us,rs = LinearMPC.simulate(mpc,x0,30)

Nmax = 30
plt = lineplot(xs[1,1:Nmax] .- xs[2,1:Nmax] .- 0.5)
plt = lineplot(-xs[1,1:Nmax] .+ xs[2,1:Nmax] .- 0.5)
