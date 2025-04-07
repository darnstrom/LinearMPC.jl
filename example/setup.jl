using LinearMPC
# Dynamics (Continuous time)
A = [0 1 0 0 0; 0 -10 9.81 0 0 ; 0 0 0 1 0 ; 0 -20 39.24 0 2; 0 0 0 0 0];
B = [0;1.0;0;2.0;0];
C = [1.0 0 0 0 0 ; 0 0 1.0 0 0];
D = [0;0];

Bw = [0;0;0;2];
# Remove measured disturbance
A = A[1:4,1:4]
B = B[1:4]
C = C[:,1:4]

# Discretize  
Ts = 0.01;
F,G,Gw = LinearMPC.zoh(A,B,Ts;Bw);
G *= 100;

# Create basic MPC 
Np,Nc = 50,10;
mpc = LinearMPC.MPC(F,G,C,Np,Nc=Nc, Gw = Gw);

# Setup cost function
Q,R,Rr= [1.2^2,1], [0.0], [1.0]
set_weights!(mpc;Q,R,Rr)

# Setup control bounds
umin,umax = [-2.0], [2.0]
set_bounds!(mpc; umin,umax)
range = LinearMPC.ParameterRange(mpc)

range.xmin[:] .= -20*ones(4);
range.xmax[:] .= 20*ones(4);
range.rmin[:] .= -20*ones(2);
range.rmax[:] .= 20*ones(2);
range.dmin[:] .= -20*ones(1);
range.dmax[:] .= 20*ones(1);

mpQP = LinearMPC.mpc2mpqp(mpc);
# Add prestabilizing
empc = LinearMPC.ExplicitMPC(mpc;range);
fig = LinearMPC.plot_regions(empc)

set_prestabilizing_feedback!(mpc);
mpc.settings.move_block = :Zero
mpQP = LinearMPC.mpc2mpqp(mpc);
