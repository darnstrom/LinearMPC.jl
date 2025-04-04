# Create Φ and Γ such that 
# X = Φ x0 + Γ U  (where X and U contains xi and ui stacked..)
function state_predictor(F,G,Np,Nc;move_block=:Hold)
    if G isa AbstractMatrix
        nx,nu = size(G);
    else
        nx,nu = size(G[1]);
    end
    Γ = zeros((Np+1)*nx,Nc*nu);
    Φ = zeros((Np+1)*nx,nx);
    Φ[1:nx,:] = I(nx);

    Gtot = G; 
    Ftot = F;
    for i in 1:Nc 
        for j in 0:(Nc-i)
            Γ[((i+j)*nx+1):(i+j+1)*nx,(j*nu+1):(j+1)*nu] = Gtot; 
        end
        Φ[i*nx+1:(i+1)*nx,:] = Ftot;
        i == Nc && break

        Ftot*=F;
        Gtot=F*Gtot
    end
    # Set ui = u_Nc for i>Nc 
    # TODO add infinite LQR gain alternative
    for i in Nc+1:Np
        Γ[(nx*i+1):nx*(i+1),:] .= F*Γ[(nx*(i-1)+1):nx*i,:];
        if(move_block==:Hold)
            Γ[(nx*i+1):nx*(i+1),end-nu+1:end] +=G;
        #else Zero
        end

        Φ[i*nx+1:(i+1)*nx,:] = F*Φ[nx*(i-1)+1:nx*i,:]
    end
    return Φ,Γ
end

# Create A u <= b+W*theta 
# correspoding to  lb<=u_i<=ub for i ∈ {1,2,...,Nc}
function create_controlbounds(mpc::MPC;n_extra_states=0)
    nu,nx = mpc.nu, mpc.nx

    if !iszero(mpc.K) || !mpc.settings.QP_double_sided
        A = Matrix{Float64}(kron(I(mpc.Nb),I(nu)))
    else
        A = nothing;
    end
    ub = repeat(mpc.umax,mpc.Nb,1)
    lb = repeat(mpc.umin,mpc.Nb,1)
    return A,ub,lb, zeros(length(ub),nx+n_extra_states) 
end

# Create A u <= b+W*theta 
# correspoding to lb <= Au*uk+ Ax*xk <=ub for k ∈ ks 
function create_general_constraints(mpc::MPC,Γ,Φ;n_extra_states=0)
    # extract data
    Np, Nc= mpc.Np, mpc.Nc
    m= length(mpc.constraints)
    nu,nx = mpc.nu, mpc.nx 

    ubtot,lbtot = zeros(0,1),zeros(0,1);
    Axtot,Autot = zeros(0,nx*(Np+1)), zeros(0,nu*Nc);
    issoft,isbinary = falses(0),falses(0)
    prios = zeros(Int,0)

    eyeX, eyeU = I(Np+1), I(Nc);
    eyeU = [eyeU;zeros(Bool,1+Np-Nc,Nc)] # Zeros address that Nc < Np (terminal state)

    for c in mpc.constraints 
        mi = size(c.Au,1);
        Ni = length(c.ks);
        Autot = [Autot; kron(eyeU[c.ks,:],c.Au)]
        Axtot = [Axtot; kron(eyeX[c.ks,:],c.Ax)]

        ubtot = [ubtot;repeat(c.ub,Ni,1)]
        lbtot = [lbtot;repeat(c.lb,Ni,1)]

        issoft = [isbinary;repeat([c.binary],mi*Ni)]
        isbinary = [isbinary;repeat([c.soft],mi*Ni)]
        prios = [prios;repeat([c.prio],mi*Ni)]
    end
    A=Axtot*Γ+Autot;
    W = -Axtot*Φ;
    W = [W zeros(size(W,1),n_extra_states)]

    return A,ubtot,lbtot,W,issoft,isbinary,prios
end

# Compute A,b,W such that the constraints for a given MPC structure
# are on the form A U<=b W th

function create_constraints(mpc,Φ,Γ;n_extra_states=0)
    n = size(Γ,2);
    A = zeros(0,n);
    bu,bl,W = zeros(0),zeros(0),zeros(0,mpc.nx);
    issoft,isbinary = falses(0),falses(0)

    # Control bounds
    if(!isempty(mpc.umax))
        Ac,bu,bl,W = create_controlbounds(mpc;n_extra_states)
        if !isnothing(Ac) A = Ac end
        issoft= falses(n);
        prios = zeros(Int,n)
        isbinary_single = falses(mpc.nu) 
        isbinary_single[mpc.binary_controls] .= true;
        isbinary = repeat(isbinary_single,mpc.Nb)
    end

    # General constraints
    if(!isempty(mpc.constraints))
        Ag,bug,blg,Wg,softg,binaryg,priog = create_general_constraints(mpc,Γ,Φ;n_extra_states);
        my = Int(size(Ag,1));
        prios = [prios;priog]
        issoft = [issoft; softg]; # Start with sofetning all general constraints
        isbinary = [isbinary; binaryg]
        bu = [bu;bug];
        bl = [bl;blg];
        A = [A;Ag];
        W = [W;Wg];
    end
    # TODO remove inf bounds...

    return A,bu,bl,W,issoft,isbinary,prios
end


# Return V and W such that
# Δu = V*u+W u_minus
# u = inv(V)*Δu - inv(V) W u_minus
# Δu' D Δu = u'V'*D*Vu+2*(V' D W u_minus)' u + u_minus W' D W u_minus
# I.e. : Hd = V'DV, f_theta = 2V'DW, Hth = W'DW
function rate_to_abs(nu,N)
    V = diagm(0=>ones(nu*N),-nu=>-ones(nu*(N-1)));
    W = [-I(nu);zeros(nu*(N-1),nu)];
    return V,W
end

# Create H, f_theta, H_theta such that the objective function for a given
# MPC problem is formulated as 0.5 U' H U+th'F_theta' U + 0.5 th' H_theta th
# (where th contains x0, r and u(k-1))
function objective(Φ,Γ,C,Q,R,S,N,Nc,nu,nx,Qf,mpc)
    pos_ids= findall(diag(Q).>0); # Ignore zero indices... (and negative)
    Q = Q[pos_ids,pos_ids];
    C = C[pos_ids,:];
    nr = size(Q,1);

    f_theta = zeros(Nc*nu,nx)

    # ==== From u' R u ====
    H = kron(I(Nc),R);


    # ==== From (Cx)'Q(Cx) ====
    Qtot = kron(I(N+1),Q);
    if(!isnothing(Qf))
        Qtot[end-nx+1:end, end-nx+1:end] .= Qf
    end
    Ctot  = kron(I(N+1),C);
    CQCtot  = kron(I(N+1),C'*Q*C);

    H += Γ'*CQCtot*Γ; 
    f_theta += Γ'*CQCtot*Φ; # from x0
    H_theta = Φ'*CQCtot*Φ

    # ==== From x' S u ====
    if(!iszero(S))
        Stot = [kron(I(Nc),S);zeros((N-Nc+1)*nx,Nc*nu)]
        if(mpc.settings.move_block==:Hold)
            Stot[Nc*nx+1:N*nx,end-nu+1:end-nu+1] = repeat(S,N-Nc,1)
        end
        GS = Γ'*Stot
        H += GS + GS'
        f_theta += Stot'*Φ
    end


    # Add regularization for binaries (won't change the solution)
    f = zeros(size(H,1),1); 
    fbin_part = zeros(mpc.nu)
    fbin_part[mpc.binary_controls] .= 1 
    fbin = repeat(fbin_part,Nc)
    f -= 0.5*fbin
    H += diagm(fbin)

    return H,f,f_theta,H_theta
end

"""
    mpQP = mpc2mpqp(mpc)

For a given MPC structure `mpc`, form the multi-parametric QP `mpQP` in the form 
```
min			0.5 U' H U+(f+f_theta*θ)' U + 0.5 th' H_theta θ 
subject to 	A U <= b + W*θ
```
"""
function mpc2mpqp(mpc::MPC)

    F,G,C = mpc.F-mpc.G*mpc.K, mpc.G, mpc.C
    Q,R,Rr,S = mpc.weights.Q, mpc.weights.R, mpc.weights.Rr, mpc.weights.S 

    nr,nx = size(mpc.C)
    nu = size(mpc.G,2)
    n_extra_states = 0

    Np,Nc = mpc.Np, mpc.Nc

    if(mpc.settings.reference_tracking) # Reference tracking -> add reference to states
        F = cat(F,I(nr),dims=(1,2))
        G = [G;zeros(nr,nu)]
        C = [C -I(nr)] 
        S = [S;zeros(nr,nu)]
        n_extra_states+= nr;
    end
    if(!iszero(Rr)) # Penalizing u -> add uold to states 
        F = cat(F,zeros(nu,nu),dims=(1,2))
        G = [G;I(nu)]
        C = cat(C,I(nu), dims=(1,2))
        Q = cat(Q,Rr,dims=(1,2))
        S = [S;-Rr];
        R+=Rr
        n_extra_states+= nu;
    end

    Φ,Γ=state_predictor(F,G,Np,Nc; move_block = mpc.settings.move_block);

    # Create objective function 
    nxe,nue= size(G) 
    H,f,f_theta,H_theta = objective(Φ,Γ,C,
                                    Q,R,S,
                                    Np,Nc,nue,nxe,mpc.weights.Qf,mpc)
    H = (H+H')/2
    # Create Constraints 
    # TODO: correctly handle extension... 
    A, bu, bl, W, issoft, isbinary, prio = create_constraints(mpc,Φ,Γ;n_extra_states)

    senses = zeros(Cint,length(bu)); 

    if(mpc.settings.soft_constraints)
        if(mpc.settings.explicit_soft && any(issoft))
            A = [A zeros(size(A,1))];
            A[issoft[size(H,1)+1:end],end].=-1;

            H = cat(H,mpc.weights.rho,dims=(1,2));
            f_theta =[f_theta;zeros(1,size(f_theta,2))];
            f = [f;0];
        end
        if(!mpc.settings.explicit_soft)
            senses[issoft[:]].+=DAQP.SOFT
        end
    end
    senses[isbinary[:]].+=DAQP.BINARY
    if(mpc.settings.QP_double_sided)
        mpQP = (H=H,f=f, H_theta = H_theta, f_theta=f_theta,
                A=Matrix{Float64}(A), bu=bu, bl=bl, W=W, senses=senses, prio=prio)
    else # Transform bl + W θ ≤ A U ≤ bu + W θ → A U ≤ b + W
        ncstr = length(bu);
        n_bounds = ncstr-size(A,1);
        bounds_table=[collect(ncstr+1:2*ncstr);collect(1:ncstr)]
        A = [A;-A]
        if(mpc.settings.explicit_soft && any(issoft))# Correct sign for slack
            A[:,end].= -abs.(A[:,end])
        end
        b = [bu;-bl]
        W = [W;-W]
        senses = [senses;senses]
        prio = [prio;prio]
        mpQP = (H=H,f=f, H_theta = H_theta, f_theta=f_theta,
                A=Matrix{Float64}(A), b=b, W=W, senses=senses,
                bounds_table=bounds_table, prio=prio)
    end
    mpc.mpQP = mpQP
    mpc.opt_model = DAQP.Model()
    if(mpc.settings.QP_double_sided)
        DAQP.setup(mpc.opt_model,H[:,:],f[:],A,bu[:],bl[:],senses[:])
    else
        DAQP.setup(mpc.opt_model,H[:,:],f[:],A,mpQP.b[:],-1e30*ones(length(mpQP.b)),senses[:])
    end
    return mpQP
end
