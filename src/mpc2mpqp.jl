# Create Φ and Γ such that 
# X = Φ x0 + Γ U  (where X and U contains xi and ui stacked..)
function state_predictor(F,G,Np,Nc)
    nx,nu = size(G);
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
    for i in Nc+1:Np
        Γ[(nx*i+1):nx*(i+1),:] .= F*Γ[(nx*(i-1)+1):nx*i,:];
        Γ[(nx*i+1):nx*(i+1),end-nu+1:end] +=G;

        Φ[i*nx+1:(i+1)*nx,:] = F*Φ[nx*(i-1)+1:nx*i,:]
    end
    return Φ,Γ
end

function get_parameter_dims(mpc::MPC)
    nr = mpc.settings.reference_tracking ?  mpc.model.ny : 0
    nuprev = iszero(mpc.weights.Rr) ? 0 : mpc.model.nu
    return mpc.model.nx,nr,mpc.model.nd,nuprev
end

# Create A u <= b+W*theta 
# correspoding to  lb<=u_i<=ub for i ∈ {1,2,...,Nc}
function create_controlbounds(mpc::MPC, Γ, Φ)
    nu,nx,Nb = mpc.model.nu, mpc.model.nx, mpc.Nc
    nth = sum(get_parameter_dims(mpc))
    # Create bounds
    ub = repeat(mpc.umax,Nb,1)
    lb = repeat(mpc.umin,Nb,1)

    #-K*X + V = -K*(Γ V + Φ x0) + V 
    #         = (I-K*Γ)V -K Φ x0
    if(!iszero(mpc.K))
        K = kron(I(Nb),[mpc.K zeros(nu,nth-nx)])
        A = (I-K*Γ[1:Nb*nth, 1:Nb*nu])
        W = K*Φ[1:Nb*nth,:]
    else
        A = mpc.settings.QP_double_sided ? zeros(0,mpc.Nc*nu) : Matrix{Float64}(kron(I(Nb),I(nu)))
        W = zeros(length(ub),nth)
    end
    ub = repeat(mpc.umax,Nb,1)
    lb = repeat(mpc.umin,Nb,1)
    return A,ub,lb, W 
end

# Create A u <= b+W*theta 
# correspoding to lb <= Au*uk+ Ax*xk <=ub for k ∈ ks 
function create_general_constraints(mpc::MPC,Γ,Φ)
    # extract data
    Np, Nc= mpc.Np, mpc.Nc
    m= length(mpc.constraints)
    nu,nx = mpc.model.nu, mpc.model.nx 
    nth = sum(get_parameter_dims(mpc))

    ubtot,lbtot = zeros(0,1),zeros(0,1);
    Axtot,Autot = zeros(0,nth*(Np+1)), zeros(0,nu*Nc);
    issoft,isbinary = falses(0),falses(0)
    prios = zeros(Int,0)

    eyeX, eyeU = I(Np+1), I(Nc);
    eyeU = [eyeU;zeros(Bool,1+Np-Nc,Nc)] # Zeros address that Nc < Np (terminal state)

    for c in mpc.constraints 
        mi = size(c.Au,1);
        Ni = length(c.ks);

        Ar = isempty(c.Ar) ? zeros(mi,mpc.nr) : c.Ar
        Ad = isempty(c.Ad) ? zeros(mi,mpc.model.nd) : c.Ad
        Aup = isempty(c.Aup) ? zeros(mi,mpc.nuprev) : c.Auip

        Autot = [Autot; kron(eyeU[c.ks,:],c.Au)]
        Axtot = [Axtot; kron(eyeX[c.ks,:],[c.Ax-c.Au*mpc.K Ar Ad Aup])]

        ubtot = [ubtot;repeat(c.ub,Ni,1)]
        lbtot = [lbtot;repeat(c.lb,Ni,1)]

        issoft = [issoft;repeat([c.soft],mi*Ni)]
        isbinary = [isbinary;repeat([c.binary],mi*Ni)]
        prios = [prios;repeat([c.prio],mi*Ni)]
    end
    A=Axtot*Γ+Autot;
    W = -Axtot*Φ;

    return A,ubtot,lbtot,W,issoft,isbinary,prios
end

# Compute A,b,W such that the constraints for a given MPC structure
# are on the form A U<=b W th
function create_constraints(mpc,Φ,Γ)

    n = size(Γ,2);
    nth = sum(get_parameter_dims(mpc))
    # Control bounds
    if(!isempty(mpc.umax))
        A,bu,bl,W = create_controlbounds(mpc,Γ,Φ)
        issoft= falses(n);
        prios = zeros(Int,n)
        isbinary_single = falses(mpc.model.nu) 
        isbinary_single[mpc.binary_controls] .= true;
        isbinary = repeat(isbinary_single,mpc.Nc)
    else
        A = zeros(0,n);
        bu,bl,W = zeros(0),zeros(0),zeros(0,nth);
        issoft,isbinary = falses(0),falses(0)
        prios = zeros(0)
    end

    # General constraints
    if(!isempty(mpc.constraints))
        Ag,bug,blg,Wg,softg,binaryg,priog = create_general_constraints(mpc,Γ,Φ);
        my = Int(size(Ag,1));
        prios = [prios;priog]
        issoft = [issoft; softg];
        isbinary = [isbinary; binaryg]
        bu = [bu;bug];
        bl = [bl;blg];
        A = [A;Ag];
        W = [W;Wg];
    end
    # TODO remove inf bounds...

    return A,bu,bl,W,issoft,isbinary,prios
end

# Create H, f_theta, H_theta such that the objective function for a given
# MPC problem is formulated as 0.5 U' H U+th'F_theta' U + 0.5 th' H_theta th
# (where th contains x0, r and u(k-1))
function objective(Φ,Γ,C,Q,R,S,Qf,N,Nc,nu,nx,mpc)

    pos_ids= findall(diag(Q).>0); # Ignore zero indices... (and negative)
    Q = Q[pos_ids,pos_ids];
    Cp = C[pos_ids,:];

    pos_ids= findall(diag(Q).>0); # Ignore zero indices... (and negative)
    Qf = Qf[pos_ids,pos_ids];
    Cf = C[pos_ids,:];


    f_theta = zeros(Nc*nu,nx)

    # ==== From u' R u ====
    H = kron(I(Nc),R);
    H[end-nu+1:end,end-nu+1:end] .+= (N-Nc)*R # To accound for Nc < N...


    # ==== From (Cx)'Q(Cx) ====
    CQCtot  = kron(I(N),Cp'*Q*Cp);
    CQCtot = cat(CQCtot,Cf'*Qf*Cf,dims=(1,2))

    H += Γ'*CQCtot*Γ; 
    f_theta += Γ'*CQCtot*Φ; # from x0
    H_theta = Φ'*CQCtot*Φ

    # ==== From x' S u ====
    if(!iszero(S))
        Stot = [kron(I(Nc),S);zeros((N-Nc+1)*nx,Nc*nu)]
        Stot[Nc*nx+1:N*nx,end-nu+1:end] = repeat(S,N-Nc,1) # Due to control horizon
        GS = Γ'*Stot
        H += (GS + GS')
        f_theta += Stot'*Φ
    end


    # Add regularization for binary variables (won't change the solution)
    f = zeros(size(H,1),1); 
    fbin_part = zeros(mpc.model.nu)
    fbin_part[mpc.binary_controls] .= 1 
    fbin = repeat(fbin_part,Nc)
    f -= 0.5*fbin
    H += diagm(fbin)

    return (H+H')/2,f,f_theta,H_theta
end

"""
    mpc2mpqp(mpc)

For a given MPC structure `mpc`, form the multi-parametric QP `mpQP`. 

"""
function mpc2mpqp(mpc::MPC)

    F,G,C = mpc.model.F-mpc.model.G*mpc.K, mpc.model.G, mpc.model.C
    Q,R,Rr,S = mpc.weights.Q, mpc.weights.R, mpc.weights.Rr, mpc.weights.S 
    Qf = iszero(mpc.weights.Qf) ? Q : mpc.weights.Qf

    nx,nr,nd,nuprev = get_parameter_dims(mpc)
    mpc.nr, mpc.nuprev =  nr,nuprev
    nu = mpc.model.nu 

    Np,Nc = mpc.Np, mpc.Nc

    if(nr > 0) # Reference tracking -> add reference to states
        F = cat(F,I(nr),dims=(1,2))
        G = [G;zeros(nr,nu)]
        C = [C -I(nr)] 
        S = [S;zeros(nr,nu)]
        #Qf = cat(Qf,zeros(nr,nr),dims=(1,2))
    end

    if(nd > 0) # add measureable disturbnace
        F = cat(F,I(nd),dims=(1,2))
        F[1:nx,end-nd+1:end] .= mpc.model.Gd
        G = [G;zeros(nd,nu)]
        S = [S;zeros(nd,nu)]
        C = [C mpc.model.Dd]
    end

    if(nuprev > 0) # Penalizing Δu -> add uold to states 
        F = cat(F,zeros(nu,nu),dims=(1,2))
        F[end-nu+1:end,1:nx] .= -mpc.K
        G = [G;I(nu)]
        C = [C zeros(nr,nu); mpc.K zeros(nu,nr+nd) I(nu)]
        Q = cat(Q,Rr,dims=(1,2))
        Qf = cat(Qf,zeros(nu,nu),dims=(1,2))
        S = [S;-Rr];
        S[1:nx,:] -=mpc.K'*Rr
        R+=Rr
    end

    if(!iszero(mpc.weights.R) && !iszero(mpc.K)) # terms from prestabilizing feedback
        Q = cat(Q,mpc.weights.R,dims=(1,2))
        Qf = cat(Qf,zeros(nu,nu),dims=(1,2))
        C = [C; mpc.K zeros(nu,nr+nd+nuprev)]
        S[1:nx,:] -=mpc.K'*mpc.weights.R
    end

    Φ,Γ=state_predictor(F,G,Np,Nc);

    # Create objective function 
    nxe,nue= size(G) 
    H,f,f_theta,H_theta = objective(Φ,Γ,C,
                                    Q,R,S,Qf,
                                    Np,Nc,nue,nxe,mpc)
    # Create Constraints 
    A, bu, bl, W, issoft, isbinary, prio = create_constraints(mpc,Φ,Γ)

    # Apply move blocking
    if(!isempty(mpc.move_blocks))
        T,id = zeros(0,0), 0
        keep_bounds = Int[] 
        for mb in mpc.move_blocks
            T = cat(T,repeat(I(mpc.model.nu),mb,1),dims=(1,2))
            keep_bounds = keep_bounds ∪ collect(id+1:id+mpc.model.nu)
            id += mb*mpc.model.nu
        end
        A, H, f, f_theta = A*T, T'*H*T, T'*f, T'*f_theta

        # remove superfluous control bounds
        keep = keep_bounds ∪ collect(mpc.model.nu*mpc.Nc+1:length(bu))
        bu,bl,W = bu[keep],bl[keep], W[keep,:]
        issoft,isbinary,prio = issoft[keep],isbinary[keep],prio[keep]
        if (!iszero(mpc.K) || !mpc.settings.QP_double_sided) # prestab feedback -> A rows for bounds
            A = A[keep,:] 
        end
    end

    # Setup sense
    senses = zeros(Cint,length(bu)); 

    # ignore constraints which have inf bounds
    for i in 1:length(bu)
        if(bu[i] > 1e20 && bl[i] < -1e20)
            senses[i] += DAQP.IMMUTABLE
        end
    end

    # Handle soft constrints  
    if(mpc.settings.soft_constraints)
        if(mpc.settings.explicit_soft && any(issoft))
            A = [A zeros(size(A,1))];
            ns = length(issoft)-size(A,1)
            A[issoft[ns+1:end],end].=-1;

            H = cat(H,mpc.settings.soft_weight,dims=(1,2));
            f_theta =[f_theta;zeros(1,size(f_theta,2))];
            f = [f;0];
        end
        if(!mpc.settings.explicit_soft)
            senses[issoft[:]].+=DAQP.SOFT
        end
    end
    # Mark binary constraints
    senses[isbinary[:]].+=DAQP.BINARY

    # Stack constraints in case of QP is assumed to be single sided. 
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
    return mpQP
end
