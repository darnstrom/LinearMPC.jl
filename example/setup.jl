using LinearMPC
# Dynamics (Continuous time)
A = [0 1 0 0 0; 0 -10 9.81 0 0 ; 0 0 0 1 0 ; 0 -20 39.24 0 2; 0 0 0 0 0];
B = 100*[0;1.0;0;2.0;0];
C = [1.0 0 0 0 0 ; 0 0 1.0 0 0];
D = [0;0];
# Remove measured disturbance
A = A[1:4,1:4]
B = B[1:4]
C = C[:,1:4]

# Discretize  
Ts = 0.01;
F,G = LinearMPC.zoh(A,B,Ts);

# Create basic MPC 
Np,Nc = 50,5;
mpc = LinearMPC.MPC(F,G,C,Np,Nc=Nc);

# Setup cost function
Q,R,Rr= [1.2^2,1], [0.0], [1.0]
set_weights!(mpc;Q,R,Rr)

# Setup control bounds
umin,umax = [-2.0], [2.0]
set_bounds!(mpc; umin,umax)

# Setup cross term
#S = randn(4,1)
#set_weights!(mpc;S)
#mpQP = LinearMPC.mpc2mpqp(mpc);
set_prestabilizing_feedback!(mpc);
mpQP = LinearMPC.mpc2mpqp(mpc);
