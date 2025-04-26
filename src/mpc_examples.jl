function mpc_examples(s, Np, Nc;params = Dict(),settings=nothing)
    if(s=="inv_pend"||s=="invpend")
        # Inverted pendulum
        A = [0 1 0 0; 
             0 -10 9.81 0; 
             0 0 0 1; 
             0 -20 39.24 0]; 
        B = 100*[0;1.0;0;2.0;;];
        C = [1.0 0 0 0; 0 0 1.0 0];
        D = [0;0];
        Ts = 0.01;

        mpc = MPC(A,B,Ts;C,Nc,Np);

        Q,R,Rr= [1.2^2,1], [0.0], [1.0]
        set_weights!(mpc;Q,R,Rr)

        umin,umax = [-2.0], [2.0]
        set_bounds!(mpc; umin,umax)

        if(isnothing(settings))
            mpc.settings.QP_double_sided = false;
        else
            mpc.settings=settings
        end

        range = ParameterRange(mpc);
        range.xmax[:] .=  20*ones(4)
        range.xmin[:] .= -20*ones(4)
        range.rmax[:] .= 20*ones(2)
        range.rmin[:] .= -20*ones(2)
        range.dmax[:] .= 20*ones(1)
        range.dmin[:] .= -20*ones(1)

    elseif(s=="dc_motor"||s=="dcmotor")
        A = [0 1.0 0 0; -51.21 -1 2.56 0; 0 0 0 1; 128 0 -6.401 -10.2];
        B = 440*[0;0;0;1.0;;];
        C = [1 0 0 0;1280 0 -64.01 0];
        Ts = 0.1;
        tau = 78.5398;
        C = C./[2*pi;2*tau]; # Scaling 


        mpc = MPC(A,B,Ts;C,Np,Nc);

        Q = [0.1^2, 0]; 
        R = [0.0]
        Rr = [0.1^2]
        set_weights!(mpc;Q,R,Rr)

        umin,umax = [-0.5],[0.5]
        set_bounds!(mpc;umin,umax)

        add_constraint!(mpc,Ax=C[2:2,:],lb = [-0.5] ,ub = [0.5], ks = 2:mpc.Nc+2, soft=true)

        if(isnothing(settings))
            mpc.settings.QP_double_sided = false;
            mpc.settings.reference_tracking=true;
            mpc.settings.explicit_soft=true;
        else
            mpc.settings=settings
        end

        range = ParameterRange(mpc);
        range.xmax[:] = [ 4*pi  4*pi/Ts  4*pi*20  4*pi*20/Ts]
        range.xmin[:] = -[ 4*pi  4*pi/Ts  4*pi*20  4*pi*20/Ts]
        range.rmax[:] .= [5;0.5]
        range.rmin[:] .=-[5;0.5]
        range.umax[:] .= [0.5023]
        range.umin[:] .= -[0.5023] 

    elseif(s=="aircraft")
        A = [-0.0151 -60.5651  0       -32.174;
             -0.0001   -1.3411  0.9929   0    ;
             0.00018  43.2541 -0.86939  0     ;
             0         0       1        0     ;]
        B = [-2.516 -13.136;-0.1689 -0.2514;-17.251 -1.5766;0 0];
        C = [0 1.0 0 0; 0 0 0 1];

        Ts = 0.05;
        F,G = zoh(A,B,Ts);
        C = C./[1;200];
        Dd = [1.0 0; 0 200]./[1;200]

        mpc = MPC(F,50*G;C,Np,Nc,Ts,Dd);


        Q = [10,10].^2; 
        R = [0.0, 0.0]
        Rr = [0.1, 0.1].^2
        set_weights!(mpc;Q,R,Rr)

        set_bounds!(mpc,umin=[-0.5;-0.5],umax=[0.5;0.5])

        set_output_bounds!(mpc, ymin = [-0.5;-0.5], ymax=[0.5;0.5],ks=2:2)

        if(isnothing(settings))
            mpc.settings.QP_double_sided= false;
            mpc.settings.reference_tracking=true;
            mpc.settings.explicit_soft=true;
        else
            mpc.settings=settings
        end

        range = ParameterRange(mpc);
        range.xmax[:] .= 20*ones(4); 
        range.xmin[:] .= -20*ones(4)
        range.dmax[:] .= 20*ones(2); 
        range.dmin[:] .= -20*ones(2)
        range.rmax[:] .= [1;0.05]; 
        range.rmin[:] .= -[1;0.05]

    elseif(s=="chained-firstorder" || s=="chained")
        nx = haskey(params,:nx) ? params[:nx] : 1
        A = -Matrix(I,nx,nx)+diagm(-1=>ones(nx-1));
        B = [1;zeros(nx-1,1);;];
        C = Matrix(I,nx,nx);
        Ts = 1;
        F,G=zoh(A,B,Ts);

        mpc = MPC(F,G;C,Np,Nc,Ts)

        Q = ones(nx);
        R = [0.0]
        Rr = [1.0]
        set_weights!(mpc;Q,R,Rr)

        set_bounds!(mpc,umin = [-1.0], umax=[1.0])
        set_output_bounds!(mpc,ymin = -10*ones(nx), ymax=10*ones(nx), ks = 2:mpc.Nc)

        if(isnothing(settings))
            mpc.settings.QP_double_sided= false;
            mpc.settings.reference_tracking=true;
            mpc.settings.explicit_soft=true;
        else
            mpc.settings=settings
        end

        range = ParameterRange(mpc);
        range.xmax[:] .= 10*ones(nx)
        range.xmin[:] .= -10*ones(nx)
        range.rmax[:] .= 10*ones(nx)
        range.rmin[:] .= -10*ones(nx)

    elseif(s=="mass-spring" || s=="mass" || s=="spring")
        κ=1; # spring
        λ=0; # damping
        nx = haskey(params,:nx) ? params[:nx] : 1
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
        C = Matrix(I,2*nm,2*nm);
        Ts = 0.5;
        F,G=zoh(A,B,Ts);

        mpc = MPC(F,G;C,Np,Nc,Ts)

        Q = 100*ones(nx);
        R = [1.0]
        Rr = [0.0]
        set_weights!(mpc;Q,R,Rr)
        set_bounds!(mpc,umin=[-0.5],umax=[0.5])
        add_constraint!(mpc, Ax = Matrix(I,nm,2*nm), 
                        lb = -4*ones(nm), ub = 4*ones(nm),
                        ks = 2:mpc.Nc)

        if(isnothing(settings))
            mpc.settings.QP_double_sided= false;
        else
            mpc.settings=settings
        end

        mpc.settings.reference_tracking= false
        range = ParameterRange(mpc);
        range.xmax[:] .= 4*ones(nx)
        range.xmin[:] .= -4*ones(nx)

    elseif(s=="nonlinear" || s == "nonlin")
        Ts = 0.2;
        F =[0.8187 zeros(1,4);
              0.1474 0.6550 -0.1637 0.0489 0.4878;
              0.01637 0.1637 0.9825 3.43*1e-3 0.0523;
              zeros(1,3) 0.8013 -0.1801;
              zeros(1,3) 0.1801 0.9813]; 
        G = [0.1813 0 0;
             0.0163 0.1637 3.43*1e-3
             1.14*1e-3 0.0175 1.77*1e-4;
             0 0 0.1801;
             0 0 0.0186;
            ]
        C  = [1.0 0 0 0 0; 0 1 2 0 0]

        mpc = MPC(F,G;C,Np,Nc,Ts)

        Q = [1.0, 1.0]
        R = zeros(3)
        Rr = (1e-1*[1, 1, 1]).^2
        set_weights!(mpc;Q,R,Rr)
        set_bounds!(mpc, umin = [-3.0,2,2], umax = [3.0,2,2])

        if(isnothing(settings))
            mpc.settings.QP_double_sided= false;
            mpc.settings.reference_tracking=true;
        else
            mpc.settings=settings
        end
        
        range = ParameterRange(mpc);
        range.xmax[:] .= [2;ones(4)]
        range.xmin[:] .= -0.5*ones(5)
        range.rmax[:] .= 10*ones(2)
        range.rmin[:] .=  -10*ones(2)

    elseif(s=="invpend_contact")
        # Inverted pendulum with contact forces
        nwalls = haskey(params,:nwalls) ? params[:nwalls] : 2
        nwalls = min(nwalls,2)

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

        C, D = Matrix{Float64}(I,4,4), zeros(4,1);
        Ts = 0.05;

        F,G = zoh(A,B,Ts);


        # MPC
        mpc = MPC(F,G;C,Np,Nc);

        Q = [1.0,1,1,1]; 
        R = [1.0;1e-4*ones(6)]
        Rr = zeros(7)
        Qf,~ = ared(mpc.model.F,mpc.model.G[:,1],mpc.weights.R[1:1,1:1],mpc.weights.Q)
        set_weights!(mpc;Q,R,Rr,Qf)

        # Control constraints
        set_bounds!(mpc,umin = [-1.0;0;zeros(4)], umax=[1.0;1e30;1e30;ones(4)])
        mpc.binary_controls = collect(4:7);

        if(isnothing(settings))
            mpc.settings.QP_double_sided= true;
            mpc.settings.reference_tracking=false;
            mpc.settings.soft_constraints=false;
        else
            mpc.settings=settings
        end

        # State constraints 
        uby = [d;pi/10;1;1];
        lby = -uby

        set_output_bounds!(mpc, ymin = lby, ymax=uby,ks = 2:mpc.Nc)
        #mpc.constraints.Ncy = [1:mpc.Nc]

        δ2l, δ2u = -uby[1]+l*lby[2]-d, -lby[1]+l*uby[2]-d
        dotδ2l, dotδ2u = -uby[3]+l*lby[4], -lby[3]+l*uby[4] 
        δ3l, δ3u = lby[1]-l*uby[2]-d, uby[1]-l*lby[2]-d
        dotδ3l, dotδ3u = lby[3]-l*uby[4], uby[3]-l*lby[4]

        u2l,u2u = κ*δ2l+ν*dotδ2l, κ*δ2u+ν*dotδ2u
        u3l,u3u = κ*δ3l+ν*dotδ3l, κ*δ3u+ν*dotδ3u

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

        add_constraint!(mpc, Au = Au2, Ax =  Ax, ub = bg2, ks = 2:mpc.Nc)
        if(nwalls == 2)
            add_constraint!(mpc, Au = Au3, Ax = -Ax, ub = bg3, ks = 2:mpc.Nc)
        end

        range = ParameterRange(mpc);
        range.xmax[:] .= 20*ones(4)
        range.xmin[:] .= -20*ones(4)
    end
    return mpc,range
end

function mpc_examples(s;settings=nothing)
    if(s=="inv_pend"||s=="invpend")
        mpc_examples(s,50,5;settings)
    elseif(s=="dc_motor"||s=="dcmotor")
        mpc_examples(s,10,2;settings)
    elseif(s=="aircraft")
        mpc_examples(s,10,2;settings)
    elseif(s=="nonlinear" || s=="nonlin")
        mpc_examples(s,5,2;settings)
    end
end
