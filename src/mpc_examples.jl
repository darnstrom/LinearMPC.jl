function mpc_examples(s, Np, Nc;nx=0,settings=nothing)
    if(s=="inv_pend"||s=="invpend")
        # Inverted pendulum
        A = [0 1 0 0 0; 0 -10 9.81 0 0 ; 0 0 0 1 0 ; 0 -20 39.24 0 2; 0 0 0 0 0];
        B = [0;1.0;0;2.0;0];
        C = [1.0 0 0 0 0 ; 0 0 1.0 0 0];
        D = [0;0];
        Ts = 0.01;

        F,G = zoh(A,B,Ts);
        G = 100*G; # TODO: Make scaling more systematic...

        # MPC

        mpc = MPC(F,G,C,Np);
        mpc.Nc = Nc;

        mpc.weights.Q = diagm([1.2^2, 1]); 
        mpc.weights.R = diagm([0.0])
        mpc.weights.Rr = diagm([1.0])

        mpc.constraints.lb = [-2.0];
        mpc.constraints.ub = [2.0];
        mpc.constraints.Ncc = mpc.Nc;

        if(isnothing(settings))
            mpc.settings.QP_double_sided = false;
        else
            mpc.settings=settings
        end

        mpQP = mpc2mpqp(mpc);
        P_theta =(A =zeros(8,0),
                  b = zeros(0),
                  lb =-[20*ones(5);20*ones(2);2],
                  ub = [20*ones(5);20*ones(2);2])

    elseif(s=="dc_motor"||s=="dcmotor")
        A = [0 1.0 0 0; -51.21 -1 2.56 0; 0 0 0 1; 128 0 -6.401 -10.2];
        B = [0;0;0;1];
        C = [1 0 0 0;1280 0 -64.01 0];
        D = [0;0];
        Ts = 0.1;
        tau = 78.5398;
        C = C./[2*pi;2*tau]; # Scaling 

        F,G= zoh(A,B,Ts);
        G = 440*G;

        mpc = MPC(F,G,C,Np);

        mpc.Nc = Nc;

        mpc.weights.Q = diagm([0.1^2, 0]); 
        mpc.weights.R = diagm([0.0])
        mpc.weights.Rr = diagm([0.1^2])

        mpc.constraints.lb = [-0.5];
        mpc.constraints.ub = [0.5];
        mpc.constraints.Ncc = mpc.Nc;

        mpc.constraints.Cy = [C[2:2,:]];
        mpc.constraints.lby = [[-0.5]];
        mpc.constraints.uby = [[0.5]];
        mpc.constraints.Ncy = [1:3];

        if(isnothing(settings))
            mpc.settings.QP_double_sided = false;
            mpc.settings.reference_tracking=true;
            mpc.settings.explicit_soft=true;
        else
            mpc.settings=settings
        end

        mpQP = mpc2mpqp(mpc);
        P_theta =(A =zeros(6,0),
                  b = zeros(0),
                  lb =-[100*ones(4);0.79577;0.5023],
                  ub = [100*ones(4);0.79577;0.5023;])
    elseif(s=="aircraft")
        A = [-0.0151 -60.5651  0       -32.174 0 0;
             -0.0001   -1.3411  0.9929   0     0 0;
             0.00018  43.2541 -0.86939  0     0 0;
             0         0       1        0     0 0;
             zeros(2,6)];
        B = [-2.516 -13.136;-0.1689 -0.2514;-17.251 -1.5766;0 0; 0 0; 0 0];
        C = [0 1 0 0 1 0;0 0 0 1 0 200];
        D = [0 0;0 0];

        Ts = 0.05;
        F,G = zoh(A,B,Ts);
        C = C./[1;200];

        mpc = MPC(F,50*G,C,Np);

        mpc.Nc = Nc;

        mpc.weights.Q = diagm([10,10].^2); 
        mpc.weights.R = diagm([0.0, 0.0])
        mpc.weights.Rr = diagm([0.1, 0.1].^2)

        mpc.constraints.lb = [-0.5;-0.5];
        mpc.constraints.ub = [0.5;0.5];
        mpc.constraints.Ncc = mpc.Nc;

        mpc.constraints.Cy = [C];
        mpc.constraints.lby = [[-0.5;-0.5]];
        mpc.constraints.uby = [[0.5;0.5]];
        mpc.constraints.Ncy = [1:1];

        if(isnothing(settings))
            mpc.settings.QP_double_sided= false;
            mpc.settings.reference_tracking=true;
            mpc.settings.explicit_soft=true;
        else
            mpc.settings=settings
        end

        mpQP = mpc2mpqp(mpc);
        P_theta =(A =zeros(10,0),
                  b = zeros(0),
                  lb =-[20*ones(6);1;0.05;0.5;0.5],
                  ub = [20*ones(6);1;0.05;0.5;0.5])
        return mpQP,P_theta,mpc
    elseif(s=="chained-firstorder" || s=="chained")
        A = -Matrix(I,nx,nx)+diagm(-1=>ones(nx-1));
        B = [1;zeros(nx-1,1)];
        C = Matrix(I,nx,nx);
        D = zeros(nx,1);
        Ts = 1;
        F,G=zoh(A,B,Ts);

        mpc = MPC(F,G,C,Np)
        mpc.Nc = Nc;

        mpc.weights.Q = diagm(ones(nx));
        #mpc.weights.Q = diagm([1.0;zeros(nx-1)]);
        mpc.weights.R = diagm([0.0])
        mpc.weights.Rr = diagm([1.0])

        mpc.constraints.lb = [-1.0]
        mpc.constraints.ub = [1.0]
        mpc.constraints.Ncc = mpc.Nc

        mpc.constraints.Cy = [C];
        mpc.constraints.lby = [-10.0*ones(nx)];
        mpc.constraints.uby = [10.0*ones(nx)];
        mpc.constraints.Ncy = [1:mpc.Nc]

        if(isnothing(settings))
            mpc.settings.QP_double_sided= false;
            mpc.settings.reference_tracking=true;
            mpc.settings.explicit_soft=true;
        else
            mpc.settings=settings
        end

        mpQP = mpc2mpqp(mpc);
        P_theta =(A =zeros(2*nx+1,0),
                  b = zeros(0),
                  lb =[-10*ones(2*nx);-1],
                  ub =[10*ones(2*nx);1])
    elseif(s=="mass-spring" || s=="mass" || s=="spring")
        κ=1; # spring
        λ=0; # damping
        nx = iseven(nx) ? nx : nx-1 # position+velocity=>2nm states
        nm = Int64(nx/2); # number of masses
        eye = Matrix(I,nm,nm)
        Fx = diagm(1=>κ*ones(nm-1), -1=>κ*ones(nm-1), 0=>-2κ*ones(nm))
        Fv = diagm(1=>λ*ones(nm-1), -1=>λ*ones(nm-1), 0=>-2λ*ones(nm))
        A = [zeros(nm,nm) Matrix(I,nm,nm);
             Fx Fv]
        B = [zeros(nm,1);
             1;
             zeros(nm-1,1)]
        display(A)
        display(B)
        C = Matrix(I,2*nm,2*nm);
        D = zeros(2*nm,1);
        Ts = 0.5;
        F,G=zoh(A,B,Ts);

        mpc = MPC(F,G,C,Np)
        mpc.Nc = Nc;

        mpc.weights.Q = diagm(100*ones(nx));
        #mpc.weights.Q = diagm([1.0;zeros(nx-1)]);
        mpc.weights.R = diagm([1.0])
        mpc.weights.Rr = diagm([0.0])

        mpc.constraints.lb = [-0.5]
        mpc.constraints.ub = [0.5]
        mpc.constraints.Ncc = mpc.Nc

        mpc.constraints.Cy = [Matrix(I,nm,2*nm)];
        mpc.constraints.lby = [-4.0*ones(nm)];
        mpc.constraints.uby = [4.0*ones(nm)];
        mpc.constraints.Ncy = [1:mpc.Nc]

        if(isnothing(settings))
            mpc.settings.QP_double_sided= false;
        else
            mpc.settings=settings
        end

        mpQP = mpc2mpqp(mpc);
        # Remove references and u_{-1} since regulation problem 
        mpQP.W = mpQP.W[:,1:nx]
        mpQP.f_theta = mpQP.f_theta[:,1:nx]
        # Remove slack
        mpQP.H=mpQP.H[1:end-1,1:end-1];
        mpQP.f_theta=mpQP.f_theta[1:end-1,:];
        mpQP.f=mpQP.f[1:end-1,:];
        mpQP.A=mpQP.A[:,1:end-1];
        # 
        P_theta =(A =zeros(nx,0),
                  b = zeros(0),
                  lb =-4.0*ones(nx),
                  ub = 4.0*ones(nx))
    elseif(s=="invpend_contact")
        # Inverted pendulum with contact forces
        mc,mp,g,l,d=1,1,10,1,0.5;
        κ,ν = 100,10;

        A = [0 0 1 0;
             0 0 0 1; 
             0 (mp*g/mc) 0 0; 
             0 (mc+mp)*g/(mc*l) 0 0];
        B = [0 0 0;
             0 0 0;
             1/mc 0 0;
             1/(mc*l) -1/(mp*l) 1/(mp*l);]
        B = [B zeros(4,4)] # Add binary constraints...

        C = I(4) 
        D = zeros(4,1);
        Ts = 0.05;

        F,G = zoh(A,B,Ts);

        # State constraints 
        uby = [d;pi/10;1;1];
        lby = -uby

        # MPC
        mpc = MPC(F,G,C,Np);
        mpc.Nc = Nc;

        mpc.weights.Q = diagm(1*[1.0,1,1,1]); 
        mpc.weights.R = diagm([1.0;1e-4*ones(6)])
        mpc.weights.Rr = diagm(zeros(7))
        Qf,~ = ared(mpc.F,mpc.G[:,1],mpc.weights.R[1:1,1:1],mpc.weights.Q)
        mpc.weights.Qf= Qf 

        mpc.constraints.lb = [-1.0;0;0;zeros(4)];
        mpc.constraints.ub = [1.0;1e30;1e30;ones(4)];
        mpc.constraints.Ncc = mpc.Nc;
        mpc.constraints.binary_controls = collect(4:7);

        if(isnothing(settings))
            mpc.settings.QP_double_sided= true;
            mpc.settings.reference_tracking=false;
            mpc.settings.soft_constraints=false;
        else
            mpc.settings=settings
        end

        mpc.constraints.Cy = [C];
        mpc.constraints.lby = [lby];
        mpc.constraints.uby = [uby];
        mpc.constraints.Ncy = [1:mpc.Nc]

        δ2l, δ2u = -uby[1]+l*lby[2]-d, -lby[1]+l*uby[2]-d
        dotδ2l, dotδ2u = -uby[3]+l*lby[4], -lby[3]+l*uby[4] 
        δ3l, δ3u = lby[1]-l*uby[2]-d, uby[1]-l*lby[2]-d
        dotδ3l, dotδ3u = lby[3]-l*uby[4], uby[3]-l*lby[4]

        u2l,u2u = κ*δ2l+ν*dotδ2l, κ*δ2u+ν*dotδ2u
        u3l,u3u = κ*δ3l+ν*dotδ3l, κ*δ3u+ν*dotδ3u

        println("δ2: $δ2l, $δ2u")
        println("dotδ2: $dotδ2l, $dotδ2u")
        println("u2: $u2l, $u2u")
        
        println("δ3: $δ3l, $δ3u")
        println("dotδ3: $dotδ3l, $dotδ3u")
        println("u3: $u3l, $u3u")

        Ax = [-1 l 0 0;
              1 -l 0 0;
              -κ κ*l -ν ν*l;
              κ -κ*l ν -ν*l;
              zeros(2,4);
              κ -κ*l ν -ν*l;
              -κ κ*l -ν ν*l;
             ]
        Au2 = [0 0 0 -δ2u 0 0 0;
               0 0 0 -δ2l 0 0 0;
               0 0 0 0 0 -u2u 0;
               0 0 0 0 0 -u2l 0;
               0 1 0 -u2u 0 0 0;
               0 1 0 0 0 -u2u 0;
               0 1 0 0 0 -u2l 0;
               0 -1 0 u2u 0 0 0;
              ]
        Au3 = [0 0 0 0 -δ3u 0 0;
               0 0 0 0 -δ3l 0 0;
               0 0 0 0 0 0 -u3u;
               0 0 0 0 0 0 -u3l;
               0 0 1 0 -u3u 0 0;
               0 0 1 0 0 0 -u3u;
               0 0 1 0 0 0 -u3l;
               0 0 -1 0 u3u 0 0;
              ]
        bg2 = [d; 
               -δ2l-d;
               κ*d;
               -κ*d-u2l;
               0;
               0;
               -u2l-κ*d;
               u2u+κ*d]
        bg3 = [d; 
               -δ3l-d;
               κ*d;
               -κ*d-u3l;
               0;
               0;
               -u3l-κ*d;
               u3u+κ*d]

        mpc.constraints.Au = [Au2;Au3]
        mpc.constraints.Ax = [Ax;-Ax]
        mpc.constraints.bg = [bg2;bg3];
        mpc.constraints.Ncg = mpc.Nc

        #ids = [1;2;3;4;7;8]
        #ids_tot = [ids;8 .+ ids];
        #mpc.constraints.Au  = mpc.constraints.Au[ids_tot,:];
        #mpc.constraints.Ax  = mpc.constraints.Ax[ids_tot,:];
        #mpc.constraints.bg  = mpc.constraints.bg[ids_tot];

        mpQP = mpc2mpqp(mpc);
        P_theta =(A =zeros(8,0),
                  b = zeros(0),
                  lb =-[20*ones(5);20*ones(2);2],
                  ub = [20*ones(5);20*ones(2);2])
    end
    return mpQP,P_theta,mpc
end

function mpc_examples(s;settings=nothing)
    if(s=="inv_pend"||s=="invpend")
        mpc_examples(s,50,5;settings)
    elseif(s=="dc_motor"||s=="dcmotor")
        mpc_examples(s,10,2;settings)
    elseif(s=="aircraft")
        mpc_examples(s,10,2;settings)
    end
end
