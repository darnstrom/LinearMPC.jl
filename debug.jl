using LinearMPC

F = [1 0.1; 0 1];
G = [0 0;1 1];
mpc = LinearMPC.MPC(F,G;C=[1 0;0 1],Np=10);
set_objective!(mpc, [1]; Q=[1,0], Rr=1e3);
set_objective!(mpc, [2]; Q=[0,1], Rr=1e3);
set_bounds!(mpc;umin=-ones(2),umax=ones(2));
move_block!(mpc,[1,1,8])
setup!(mpc)

sim_imp = Simulation(mpc;x0=10*ones(2), r = [10,0], N=500);
