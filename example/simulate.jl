## Normal invpend
using LinearMPC

opts = LinearMPC.MPCSettings()
mpQP,P_theta,mpc = LinearMPC.mpc_examples("invpend",settings=opts)

x0 = [0.0,0.0,0.0,0.0,-1]; 
rs = [zeros(2,20) repeat([20;0],1,80)];
xs,us,rs = LinearMPC.simulate(mpc,x0,1000;r=rs)

## Contact invpend 
using LinearMPC
using UnicodePlots
Nsim = 250 

#opts = LinearMPC.MPCSettings()
#opts.reference_tracking = false
#opts.explicit_soft = true 
#mpQP,P_theta,mpc = LinearMPC.mpc_examples("invpend_contact",5,5;settings=opts)
mpQP,P_theta,mpc = LinearMPC.mpc_examples("invpend_contact",10,10)

x0 = [0.0, 0.0, 1, 0.0]; 
#rs = [zeros(2,50) repeat([10;0],1,250)];
#xs,us,rs = LinearMPC.simulate(mpc,x0,3000;r=rs)
xs,us,rs,diffs = LinearMPC.simulate(mpc,x0,Nsim)

## Plot
plt = lineplot(xs[1,1:Nsim] .- xs[2,1:Nsim] .- 0.5)
plt = lineplot(-xs[1,1:Nsim] .+ xs[2,1:Nsim] .- 0.5)
plt = lineplot(us[1,1:Nsim])
