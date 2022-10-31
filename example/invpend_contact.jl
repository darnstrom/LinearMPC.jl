## Init 
using LinearMPC
using UnicodePlots
using LinearAlgebra
using MatrixEquations
Nsim = 250 
Np,Nc = 10,10

## Formulate MPC problem
mc,mp,g,l,d=1,1,10,1,0.5;
κ,ν = 100,10;

A = [0 0 1 0;
     0 0 0 1; 
     0 (mp*g/mc) 0 0; 
     0 (mc+mp)*g/(mc*l) 0 0];
B = [0 0;
     0 0;
     1/mc 0;
     1/(mc*l) -1/(mp*l);]
B = [B zeros(4,2)]; 

C = I(4) 
D = zeros(4,1);
Ts = 0.05;

F,G = LinearMPC.zoh(A,B,Ts);

# State constraints 
uby = [d;pi/10;2;1];
lby = -uby

# MPC
mpc = LinearMPC.MPC(F,G,C,Np);
mpc.Nc = Nc;

mpc.weights.Q = diagm(1*[1.0,1,1,1]); 
mpc.weights.R = diagm([1.0;1e-4*ones(3)])
mpc.weights.Rr = diagm(1e-2*ones(4))
Qf,~ = ared(mpc.F,mpc.G[:,1],mpc.weights.R[1:1,1:1],mpc.weights.Q)
mpc.weights.Qf= Qf 

mpc.constraints.lb = [-1.0;0;zeros(2)];
mpc.constraints.ub = [1.0;1e30;ones(2)];
mpc.constraints.Ncc = mpc.Nc;
mpc.constraints.binary_controls = collect(3:4);

mpc.settings.QP_double_sided= true;
mpc.settings.reference_tracking=false;
mpc.settings.soft_constraints=false;

mpc.constraints.Cy = [C];
mpc.constraints.lby = [lby];
mpc.constraints.uby = [uby];
mpc.constraints.Ncy = [1:mpc.Nc]

δ2l, δ2u = -uby[1]+l*lby[2]-d, -lby[1]+l*uby[2]-d
dotδ2l, dotδ2u = -uby[3]+l*lby[4], -lby[3]+l*uby[4] 

u2l,u2u = κ*δ2l+ν*dotδ2l, κ*δ2u+ν*dotδ2u

Ax = [-1 l 0 0;
      1 -l 0 0;
      -κ κ*l -ν ν*l;
      κ -κ*l ν -ν*l;
      zeros(2,4);
      κ -κ*l ν -ν*l;
      -κ κ*l -ν ν*l;
     ]
Au = [0 0 -δ2u 0;
       0 0 -δ2l 0;
       0 0 0 -u2u;
       0 0 0 -u2l;
       0 1 -u2u 0;
       0 1 0 -u2u;
       0 1 0 -u2l;
       0 -1 u2u 0;
      ]
bg = [d; 
       -δ2l-d;
       κ*d;
       -κ*d-u2l;
       0;
       0;
       -u2l-κ*d;
       u2u+κ*d]

mpc.constraints.Au = Au;
mpc.constraints.Ax = Ax;
mpc.constraints.bg = bg;
mpc.constraints.Ncg = mpc.Nc;

mpQP = LinearMPC.mpc2mpqp(mpc);

## Simulate
x0 = [0.0, 0.0, -1, 0.0]; 
xs,us,rs,diffs = LinearMPC.simulate(mpc,x0,Nsim)

## Plot
plt = lineplot(-xs[1,1:Nsim] .+ xs[2,1:Nsim] .- 0.5)
plt = lineplot(us[1,1:Nsim])
