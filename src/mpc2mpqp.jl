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
    if mpc.settings.reference_preview && !mpc.settings.reference_condensation && nr > 0
        nr = nr * mpc.Np  # Reference preview uses Np time steps
    end
    nuprev = iszero(mpc.weights.Rr) ? 0 : mpc.model.nu
    return mpc.model.nx,nr,mpc.model.nd,nuprev
end

function get_parameter_names(mpc::MPC)
    nx,nr,nd,nuprev = get_parameter_dims(mpc)
    names = copy(mpc.model.labels.x)
    if nr > 0
        if mpc.settings.reference_preview && !mpc.settings.reference_condensation
            # For reference preview, create names for each time step
            for k in 0:mpc.Np-1
                push!(names, Symbol.(string.(mpc.model.labels.y).*"r_$k")...)
            end
        else
            push!(names,Symbol.(string.(mpc.model.labels.y).*"r")...)
        end
    end
    nd>0 && push!(names,mpc.model.labels.d...)
    nuprev>0 && push!(names,Symbol.(string.(mpc.model.labels.u).*"p")...)
    return names
end

# Create A u <= b+W*theta 
# correspoding to  lb<=u_i<=ub for i ∈ {1,2,...,Nc}
function create_controlbounds(mpc::MPC, Γ, Φ)
    nu,nx,Nb = mpc.model.nu, mpc.model.nx, mpc.Nc
    nth = sum(get_parameter_dims(mpc))

    #-K*X + V = -K*(Γ V + Φ x0) + V 
    #         = (I-K*Γ)V -K Φ x0
    if(!iszero(mpc.K))
        K = kron(I(Nb),[mpc.K zeros(nu,nth-nx)])
        A = (I-K*Γ[1:Nb*nth, 1:Nb*nu])
        W = K*Φ[1:Nb*nth,:]
    else
        A = zeros(0,mpc.Nc*nu)
        W = zeros(Nb*nu,nth)
    end
    # Create bounds
    ub = repeat(mpc.umax,Nb,1)
    lb = repeat(mpc.umin,Nb,1)

    # Tighten constraint
    if(!iszero(mpc.K) && !iszero(mpc.model.wmin) && !iszero(mpc.model.wmax))
        FK= mpc.model.F-mpc.model.G*mpc.K
        ut,lt= constraint_tightening(-mpc.K,FK,1:Nb,mpc.model.wmin,mpc.model.wmax)
        ub -= ut
        lb += lt
    end
    return A,ub,lb, W 
end

# Create A u <= b+W*theta 
# correspoding to lb <= Au*uk+ Ax*xk <=ub for k ∈ ks 
function create_general_constraints(mpc::MPC,Γ,Φ)
    # extract data
    Np, Nc= mpc.Np, mpc.Nc
    m= length(mpc.constraints)
    nu,nx = mpc.model.nu, mpc.model.nx 

    if mpc.settings.reference_preview
        nxe = sum(get_parameter_dims(mpc))-mpc.nr
        nrx = 0
    else
        nxe = sum(get_parameter_dims(mpc))
        nrx = mpc.nr
    end

    ubtot,lbtot = zeros(0,1),zeros(0,1);
    Axtot,Autot = zeros(0,nxe*(Np+1)), zeros(0,nu*Nc);
    issoft,isbinary = falses(0),falses(0)
    prios = zeros(Int,0)

    eyeX, eyeU = I(Np+1), I(Nc);
    eyeU = [eyeU;zeros(Bool,1+Np-Nc,Nc)] # Zeros address that Nc < Np (terminal state)

    tighten_constraints = !iszero(mpc.model.wmin) && !iszero(mpc.model.wmax)
    FK = tighten_constraints ? mpc.model.F-mpc.model.G*mpc.K : zeros(0,0)

    for c in mpc.constraints 
        mi = size(c.Au,1);
        ks = [k for k in c.ks if k<= Np]
        Ni = length(ks);

        Ax = c.Ax-c.Au*mpc.K
        Ar = isempty(c.Ar) ? zeros(mi,nrx) : c.Ar
        Ad = isempty(c.Ad) ? zeros(mi,mpc.model.nd) : c.Ad
        Aup = isempty(c.Aup) ? zeros(mi,mpc.nuprev) : c.Auip

        Autot = [Autot; kron(eyeU[ks,:],c.Au)]
        Axtot = [Axtot; kron(eyeX[ks,:],[Ax Ar Ad Aup])]

        ubtot = [ubtot;repeat(c.ub,Ni,1)]
        lbtot = [lbtot;repeat(c.lb,Ni,1)]

        if(tighten_constraints)
            ut,lt= constraint_tightening(Ax,FK,ks,mpc.model.wmin,mpc.model.wmax)
            ubtot -= ut
            lbtot += lt
        end

        issoft = [issoft;repeat([c.soft],mi*Ni)]
        isbinary = [isbinary;repeat([c.binary],mi*Ni)]
        prios = [prios;repeat([c.prio],mi*Ni)]
    end
    A=Axtot*Γ+Autot;
    W = -Axtot*Φ;

    # Correct for r not being state when using reference preview
    if mpc.settings.reference_preview
        W = [W[:,1:nx] zeros(size(W,1),mpc.nr) W[:,nx+1:end]]
    end

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

    return A,bu,bl,W,issoft,isbinary,prios
end

# Create H, f_theta, H_theta such that the objective function for a given
# MPC problem is formulated as 0.5 U' H U+th'F_theta' U + 0.5 th' H_theta th
# (where th contains x0, r and u(k-1))
function objective(Φ,Γ,C,Q,R,S,Qf,N,Nc,nu,nx,mpc)

    ny = mpc.model.ny
    Q_full,Qf_full = Q[1:ny,1:ny],Qf[1:ny,1:ny]
    C_full = C[1:ny,:]

    pos_ids_Q = findall(diag(Q).>0); # Ignore zero indices... (and negative)
    Q = Q[pos_ids_Q,pos_ids_Q];
    Cp = C[pos_ids_Q,:];

    pos_ids_Qf = findall(diag(Qf).>0); # Ignore zero indices... (and negative)
    Qf = Qf[pos_ids_Qf,pos_ids_Qf];
    Cf = C[pos_ids_Qf,:];

    # Get parameter dimensions
    nxp, nrp, ndp, nup = get_parameter_dims(mpc)

    # ==== From u' R u ====
    H = kron(I(Nc),R);
    H[end-nu+1:end,end-nu+1:end] .+= (N-Nc)*R # To accound for Nc < N...

    # ==== From (Cx)'Q(Cx) ====
    CQCtot  = kron(I(N),Cp'*Q*Cp);
    CQCf = Cf'*Qf*Cf + cat(mpc.weights.Qfx,zeros(nx-nxp,nx-nxp), dims=(1,2))
    CQCtot = cat(CQCtot,CQCf,dims=(1,2))

    H += Γ'*CQCtot*Γ; 
    # f_theta & H_theta for state parameters
    f_theta  = Γ'*CQCtot*Φ; # from x0
    H_theta  = Φ'*CQCtot*Φ

    # ==== From x' S u ====
    if(!iszero(S))
        Stot = [kron(I(Nc),S);zeros((N-Nc+1)*nx,Nc*nu)]
        Stot[Nc*nx+1:N*nx,end-nu+1:end] = repeat(S,N-Nc,1) # Due to control horizon
        GS = Γ'*Stot
        H += (GS + GS')
        f_theta += Stot'*Φ
    end

    # ==== Reference tracking terms ====
    if mpc.settings.reference_tracking && nrp > 0
        if mpc.settings.reference_preview
            # Reference preview mode: handle time-varying references
            # Needs to add terms for r to f_theta and H_theta
            # Recall that θ = [x0 r nd uprev]
            if nrp > 0
                Fr = -Γ'*cat(kron(I(N),C_full'*Q_full), C_full'*Qf_full,dims=(1,2))
                Fr = Fr[:,ny+1:end] # First reference superfluous
                Hr = cat(kron(I(N-1),Q_full),Qf_full,dims=(1,2))
                if mpc.settings.reference_condensation
                    Is = repeat(I(ny),mpc.Np)
                    mpc.traj2setpoint = (inv(H)*Fr*Is)\(inv(H)*Fr)
                    #mpc.traj2setpoint = (Fr*Is)\(Fr)
                    Fr *= Is
                    Hr = Is'*Hr*Is
                end
                f_theta = [f_theta[:,1:nxp] Fr f_theta[:,nxp+1:end]]

                H_theta = [H_theta[1:nxp, 1:nxp] zeros(nxp,nrp) H_theta[1:nxp,nxp+1:end];
                           zeros(nrp,nxp) Hr zeros(nrp,ndp+nup);
                           H_theta[nxp+1:end,1:nxp] zeros(ndp+nup,nrp) H_theta[nxp+1:end,nxp+1:end]]
            end
        else
            # Standard mode: references are handled as augmented states
            # No additional terms needed here since references are in the state vector
        end
    end


    # Add regularization for binary variables (won't change the solution)
    f = zeros(size(H,1),1); 
    fbin_part = zeros(mpc.model.nu)
    fbin_part[mpc.binary_controls] = (mpc.umax[mpc.binary_controls] + mpc.umin[mpc.binary_controls])/2
    fbin = repeat(fbin_part,Nc)
    f -= fbin
    H += diagm(fbin .!= 0)

    return (H+H')/2,f,f_theta,H_theta
end

"""
    mpc2mpqp(mpc)

For a given MPC structure `mpc`, form the multi-parametric QP `mpQP`. 

"""
function mpc2mpqp(mpc::MPC; singlesided=false, single_soft=false)

    F,G,C = mpc.model.F-mpc.model.G*mpc.K, mpc.model.G, mpc.model.C
    Q,R,Rr,S = mpc.weights.Q, mpc.weights.R, mpc.weights.Rr, mpc.weights.S 
    Qf = iszero(mpc.weights.Qf) && iszero(mpc.weights.Qfx) ? Q : mpc.weights.Qf

    nx,nr,nd,nuprev = get_parameter_dims(mpc)
    mpc.nr, mpc.nuprev =  nr,nuprev
    nu = mpc.model.nu 

    Np,Nc = mpc.Np, mpc.Nc

    if(nr > 0) # Reference tracking -> add reference to states
        if mpc.settings.reference_preview
            # Reference preview mode: references are part of parameter vector
            # No need to add reference states to F, G matrices
            # References will be handled directly in the objective function
        else
            # Standard mode: add reference as constant states
            F = cat(F,I(mpc.model.ny),dims=(1,2))
            G = [G;zeros(mpc.model.ny,nu)]
            C = [C -I(mpc.model.ny)] 
            S = [S;zeros(mpc.model.ny,nu)]
        end
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
        nye,nxe = size(C)
        C = [C zeros(nye,nu); mpc.K zeros(nu,nxe-nx) I(nu)]
        Q = cat(Q,Rr,dims=(1,2))
        Qf = cat(Qf,zeros(nu,nu),dims=(1,2))
        S = [S;-Rr];
        S[1:nx,:] -=mpc.K'*Rr
        R+=Rr
    end

    if(!iszero(mpc.weights.R) && !iszero(mpc.K)) # terms from prestabilizing feedback
        Q = cat(Q,mpc.weights.R,dims=(1,2))
        Qf = cat(Qf,zeros(nu,nu),dims=(1,2))
        C = [C; mpc.K zeros(nu,size(C,2)-nx)]
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
        for mb in [mpc.move_blocks[1:end-1];1]
            T = cat(T,repeat(I(mpc.model.nu),mb,1),dims=(1,2))
            keep_bounds = keep_bounds ∪ collect(id+1:id+mpc.model.nu)
            id += mb*mpc.model.nu
        end
        A, H, f, f_theta = A*T, T'*H*T, T'*f, T'*f_theta

        # remove superfluous control bounds
        keep = keep_bounds ∪ collect(mpc.model.nu*mpc.Nc+1:length(bu))
        bu,bl,W = bu[keep],bl[keep], W[keep,:]
        issoft,isbinary,prio = issoft[keep],isbinary[keep],prio[keep]
        if (!iszero(mpc.K)) # prestab feedback -> A rows for bounds
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

    # Mark soft constrints  
    senses[issoft[:]].+=DAQP.SOFT

    # Mark binary constraints
    senses[isbinary[:]].+=DAQP.BINARY

    # Replace Infs for codegen
    clamp!(bu,-1e30,1e30)
    clamp!(bl,-1e30,1e30)

    # Stack constraints in case of QP is assumed to be single sided. 
    mpQP = (H=H,f=f, H_theta = H_theta, f_theta=f_theta,
            A=Matrix{Float64}(A), bu=bu, bl=bl, W=W, senses=senses,
            prio=prio, has_binaries=any(isbinary))

    if(singlesided)
        mpQP = make_singlesided(mpQP;single_soft,soft_weight=mpc.settings.soft_weight)
    end
    return mpQP
end
