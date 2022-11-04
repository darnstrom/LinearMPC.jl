## Init 
using LinearMPC
using UnicodePlots
using LinearAlgebra
using MatrixEquations
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
Ts = 0.1;

F,G = LinearMPC.zoh(A,B,Ts);

# State constraints 
uby = [d;pi/10;2;1];
lby = -uby

# MPC
mpc = LinearMPC.MPC(F,G,C,Np);
mpc.Nc = Nc;

mpc.weights.Q = diagm(1*[1.0,1,1,1]); 
#mpc.weights.R = diagm([1.0;1e-4*ones(3)])
#mpc.weights.Rr = diagm(1e-2*ones(4))
mpc.weights.R = diagm([1.0;1e-2;zeros(2)])
mpc.weights.Rr = diagm([1e-2;zeros(3)])
Qf,~ = ared(mpc.F,mpc.G[:,1],mpc.weights.R[1:1,1:1],mpc.weights.Q)
mpc.weights.Qf= Qf 

mpc.constraints.lb = [-1.0;0;zeros(2)];
mpc.constraints.ub = [1.0;1e30;ones(2)];
mpc.constraints.Ncc = mpc.Nc;
mpc.constraints.binary_controls = collect(3:4);

mpc.settings.QP_double_sided= true;
mpc.settings.reference_tracking=false;
mpc.settings.soft_constraints=false;

Cy = C
#Cy = [1.0 0 0 0; 0 1 0 0];
mpc.constraints.Cy = [Cy];
mpc.constraints.lby = [lby];
mpc.constraints.uby = [uby];
#mpc.constraints.lby = [lby[1:2]];
#mpc.constraints.uby = [uby[1:2]];
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
mpLDP = LinearMPC.dualize(mpQP,1)

## Simulate
Nsim = 50 
x0 = [0.0, pi/40, -1, 0.0]; 
xs,us,rs,diffs = LinearMPC.simulate(mpc,x0,Nsim)

## Plot
plt = lineplot(-xs[1,1:Nsim] .+ xs[2,1:Nsim] .- 0.5)
plt = lineplot(us[1,1:Nsim])
plt = lineplot(us[2,1:Nsim])
## Code generation
#opts = Dict(:primal_tol => 1e-4)
opts = Dict(:zero_tol => 1e-7,
            :primal_tol => 1e-4,
            :dual_tol => 1e-4,
            :cycle_tol => 25,
            :iter_limit => 1e4)
LinearMPC.codegen(mpc,opt_settings=opts)
## Connect to Nucleo
using LibSerialPort
nucleo = LibSerialPort.open("/dev/ttyACM1",115200) 
u = zeros(4);
x = [0.0, pi/40, -1, 0];
Nsteps = 50;
exitflags = zeros(Nsteps);
cycles = zeros(Nsteps);
n_nodes = zeros(Nsteps);
n_iters = zeros(Nsteps);
us = zeros(mpc.nu,Nsteps)
xs = zeros(mpc.nx,Nsteps)
data = read(nucleo)
#close(nucleo)
## Step on Nucleo
for k = 1:Nsteps
    xs[:,k] = x;
    display(k)
    write(nucleo,Cfloat.(x));
    #raw_data=read(nucleo)
    wait_count = 0
    while (bytesavailable(nucleo) < 32)
        if(wait_count > 10)
            println("trying to resend!")
            data_flush = read(nucleo)
            write(nucleo,Cfloat.(x));
            wait_count = 0
        end
        sleep(0.25)
        wait_count +=1
    end
    for i = 1:4
        u[i] = read(nucleo,Cfloat)
    end
    cycles[k] = read(nucleo,Cfloat)
    exitflags[k] = read(nucleo,Cfloat)
    n_nodes[k] = read(nucleo,Cfloat)
    n_iters[k] = read(nucleo,Cfloat)
    us[:,k] = u
    @assert(exitflags[k]==1)
    #display(u)
    # step
    x = mpc.F*x+mpc.G*u
end

