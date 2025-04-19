using LinearMPC
# Example of setting up the inverted pendulum example
A = [0 1 0 0; 0 -10 9.81 0 ; 0 0 0 1 ; 0 -20 39.24 0];
B = 100*[0;1.0;0;2.0];
C = [1.0 0 0 0 ; 0 0 1.0 0];
Bw = [0;0;0;2];
# Discretize  
Ts = 0.01;
F,G,Gw = LinearMPC.zoh(A,B,Ts;Bw);

# Create basic MPC 
Np,Nc = 50,10;
mpc = LinearMPC.MPC(F,G,C,Np,Nc=Nc, Gw = Gw);

# Setup cost function
Q,R,Rr= [1.2^2,1], [0.0], [1.0]
set_weights!(mpc;Q,R,Rr)

# Setup control bounds
umin,umax = [-2.0], [2.0]
set_bounds!(mpc; umin,umax)

# Setup a ParameterRange for EMPC
range = LinearMPC.ParameterRange(mpc)

range.xmin[:] .= -20*ones(4);
range.xmax[:] .= 20*ones(4);
range.rmin[:] .= -20*ones(2);
range.rmax[:] .= 20*ones(2);
range.dmin[:] .= -20*ones(1);
range.dmax[:] .= 20*ones(1);

# Compute explicit
empc = LinearMPC.ExplicitMPC(mpc;range);
