# Create Φ and Γ such that 
# X = Φ x0 + Γ U  (where X and U contains xi and ui stacked..)
function state_predictor(F,G,Np,Nc)
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
    for i in Nc+1:Np
        Γ[(nx*i+1):nx*(i+1),:] .= F*Γ[(nx*(i-1)+1):nx*i,:];
        Γ[(nx*i+1):nx*(i+1),end-nu+1:end] +=G;

        Φ[i*nx+1:(i+1)*nx,:] = F*Φ[nx*(i-1)+1:nx*i,:]
    end
    return Φ,Γ
end

function get_parameter_dims(mpc::MPC)
    nr,nx = size(mpc.C)
    nw = iszero(mpc.Gw) ? 0 : size(mpc.Gw,2)
    nd = iszero(mpc.Dd) ? 0 : size(mpc.Dd,2)
    nuprev = iszero(mpc.weights.Rr) ? 0 : mpc.nu
    return nx,nr,nw,nd,nuprev
end

# Create A u <= b+W*theta 
# correspoding to  lb<=u_i<=ub for i ∈ {1,2,...,Nc}
function create_controlbounds(mpc::MPC, Γ, Φ)
    nu,nx,Nb = mpc.nu, mpc.nx, mpc.Nb
    n_extra_states = mpc.nr + mpc.nw + mpc.nd + mpc.nuprev
    # Create bounds
    ub = repeat(mpc.umax,Nb,1)
    lb = repeat(mpc.umin,Nb,1)

    #-K*X + V = -K*(Γ V + Φ x0) + V 
    #         = (I-K*Γ)V -K Φ x0
    if(!iszero(mpc.K))
        K = kron(I(Nb),[mpc.K zeros(nu,n_extra_states)])
        A = (I-K*Γ[1:Nb*(nx+n_extra_states), 1:Nb*nu])
        W = K*Φ[1:Nb*(nx+n_extra_states),:]
    else
        A = mpc.settings.QP_double_sided ? nothing : Matrix{Float64}(kron(I(mpc.Nb),I(nu)))
        W = zeros(length(ub),nx+n_extra_states)
    end
    ub = repeat(mpc.umax,mpc.Nb,1)
    lb = repeat(mpc.umin,mpc.Nb,1)
    return A,ub,lb, W 
end

# Create A u <= b+W*theta 
# correspoding to lb <= Au*uk+ Ax*xk <=ub for k ∈ ks 
function create_general_constraints(mpc::MPC,Γ,Φ)
    # extract data
    Np, Nc= mpc.Np, mpc.Nc
    m= length(mpc.constraints)
    nu,nx = mpc.nu, mpc.nx 
    n_extra_states = mpc.nr + mpc.nw + mpc.nd + mpc.nuprev

    ubtot,lbtot = zeros(0,1),zeros(0,1);
    Axtot,Autot = zeros(0,(nx+n_extra_states)*(Np+1)), zeros(0,nu*Nc);
    issoft,isbinary = falses(0),falses(0)
    prios = zeros(Int,0)

    eyeX, eyeU = I(Np+1), I(Nc);
    eyeU = [eyeU;zeros(Bool,1+Np-Nc,Nc)] # Zeros address that Nc < Np (terminal state)

    for c in mpc.constraints 
        mi = size(c.Au,1);
        Ni = length(c.ks);

        Ar = isempty(c.Ar) ? zeros(mi,mpc.nr) : c.Ar
        Aw = isempty(c.Aw) ? zeros(mi,mpc.nw) : c.Aw
        Ad = isempty(c.Ad) ? zeros(mi,mpc.nd) : c.Ad
        Aup = isempty(c.Aup) ? zeros(mi,mpc.nuprev) : c.Auip

        Autot = [Autot; kron(eyeU[c.ks,:],c.Au)]
        Axtot = [Axtot; kron(eyeX[c.ks,:],[c.Ax-c.Au*mpc.K Ar Aw Ad Aup])]

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
    A = zeros(0,n);
    bu,bl,W = zeros(0),zeros(0),zeros(0,mpc.nx);
    issoft,isbinary = falses(0),falses(0)

    # Control bounds
    if(!isempty(mpc.umax))
        Ac,bu,bl,W = create_controlbounds(mpc,Γ,Φ)
        if !isnothing(Ac) A = Ac end
        issoft= falses(n);
        prios = zeros(Int,n)
        isbinary_single = falses(mpc.nu) 
        isbinary_single[mpc.binary_controls] .= true;
        isbinary = repeat(isbinary_single,mpc.Nb)
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
    fbin_part = zeros(mpc.nu)
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

    F,G,C = mpc.F-mpc.G*mpc.K, mpc.G, mpc.C
    Q,R,Rr,S = mpc.weights.Q, mpc.weights.R, mpc.weights.Rr, mpc.weights.S 
    Qf = isempty(mpc.weights.Qf) ? Q : mpc.weights.Qf

    nx,nr,nw,nd,nuprev = get_parameter_dims(mpc)
    mpc.ny, mpc.nr, mpc.nw, mpc.nd, mpc.nuprev =  size(C,1),nr,nw,nd,nuprev
    nu = size(mpc.G,2)

    Np,Nc = mpc.Np, mpc.Nc

    if(nr > 0) # Reference tracking -> add reference to states
        F = cat(F,I(nr),dims=(1,2))
        G = [G;zeros(nr,nu)]
        C = [C -I(nr)] 
        S = [S;zeros(nr,nu)]
        #Qf = cat(Qf,zeros(nr,nr),dims=(1,2))
    end

    if(nw > 0) # add measureable input disturbance
        F = cat(F,I(nw),dims=(1,2))
        F[1:nx,end-nw+1:end] .= mpc.Gw
        G = [G;zeros(nw,nu)]
        #Qf = cat(Qf,zeros(nw,nw),dims=(1,2))
        S = [S;zeros(nw,nu)]
        C = [C zeros(size(C,1),nw)]
    end

    if(nd > 0) # add measureable output disturbnace
        F = cat(F,I(nd),dims=(1,2))
        G = [G;zeros(nd,nu)]
        S = [S;zeros(nd,nu)]
        C = [C mpc.Dd]
    end

    if(nuprev > 0) # Penalizing Δu -> add uold to states 
        F = cat(F,zeros(nu,nu),dims=(1,2))
        F[end-nu+1:end,1:nx] .= -mpc.K
        G = [G;I(nu)]
        C = [C zeros(nr,nu); mpc.K zeros(nu,nr+nd+nw) I(nu)]
        Q = cat(Q,Rr,dims=(1,2))
        Qf = cat(Qf,zeros(nu,nu),dims=(1,2))
        S = [S;-Rr];
        S[1:nx,:] -=mpc.K'*Rr
        R+=Rr
    end

    if(!iszero(mpc.weights.R) && !iszero(mpc.K)) # terms from prestabilizing feedback
        Q = cat(Q,mpc.weights.R,dims=(1,2))
        Qf = cat(Qf,zeros(nu,nu),dims=(1,2))
        C = [C; mpc.K zeros(nu,nr+nw+nd+nuprev)]
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
            T = cat(T,repeat(I(mpc.nu),mb,1),dims=(1,2))
            keep_bounds = keep_bounds ∪ collect(id+1:id+mpc.nu)
            id += mb*mpc.nu
        end
        A, H, f, f_theta = A*T, T'*H*T, T'*f, T'*f_theta

        # remove superfluous control bounds
        keep = keep_bounds ∪ collect(mpc.nu*mpc.Nc+1:length(bu))
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

            H = cat(H,mpc.weights.rho,dims=(1,2));
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
