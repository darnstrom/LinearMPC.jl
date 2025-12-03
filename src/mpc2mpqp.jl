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
    # Use stored values if QP is set up, otherwise compute from settings
    # This ensures consistency between the QP and parameter vector at runtime
    if mpc.mpqp_issetup
        return mpc.model.nx, mpc.nr, mpc.model.nd, mpc.nuprev, mpc.nl
    end
    nr = mpc.settings.reference_tracking ?  mpc.model.ny : 0
    if mpc.settings.reference_preview && !mpc.settings.reference_condensation && nr > 0
        nr = nr * mpc.Np  # Reference preview uses Np time steps
    end
    nuprev = iszero(mpc.weights.Rr) ? 0 : mpc.model.nu
    nl = mpc.settings.linear_cost ? mpc.model.nu * mpc.Nc : 0
    return mpc.model.nx,nr,mpc.model.nd,nuprev,nl
end

function get_parameter_names(mpc::MPC)
    nx,nr,nd,nuprev,nl = get_parameter_dims(mpc)
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
    if nl > 0
        for k in 0:mpc.Nc-1
            push!(names, Symbol.(string.(mpc.model.labels.u).*"l_$k")...)
        end
    end
    return names
end

# Create A u <= b+W*theta
# correspoding to  lb<=u_i<=ub for i ∈ {1,2,...,Nc}
function create_controlbounds(mpc::MPC, Γ, Φ)
    nu,nx,Nb = mpc.model.nu, mpc.model.nx, mpc.Nc
    _,_,_,_,nl = get_parameter_dims(mpc)
    nth = sum(get_parameter_dims(mpc)) - nl  # Exclude linear cost from state-related dimensions
    !iszero(mpc.model.offset) && (nth+=1); # constant offset in dynamics

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
    # Append zero columns for linear cost parameters (they don't affect constraints)
    nl > 0 && (W = [W zeros(size(W,1), nl)])
    # Create bounds
    ub = repeat(mpc.umax,Nb,1)
    lb = repeat(mpc.umin,Nb,1)

    # Tighten constraint
    if(!iszero(mpc.K) && (!iszero(mpc.model.wmin) || !iszero(mpc.model.wmax)))
        FK= mpc.model.F-mpc.model.G*mpc.K
        ut,lt= constraint_tightening(-mpc.K,FK,1:Nb,mpc.model.wmin,mpc.model.wmax,mpc.Δx0)
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
    _,_,_,_,nl = get_parameter_dims(mpc)

    if mpc.settings.reference_preview
        nxe = sum(get_parameter_dims(mpc))-mpc.nr-nl
        nrx = 0
    else
        nxe = sum(get_parameter_dims(mpc))-nl
        nrx = mpc.nr
    end
    !iszero(mpc.model.offset) && (nxe+=1); # constant offset in dynamics

    ubtot,lbtot = zeros(0,1),zeros(0,1);
    Axtot,Autot = zeros(0,nxe*(Np+1)), zeros(0,nu*Nc);
    issoft,isbinary = falses(0),falses(0)
    prios = zeros(Int,0)

    eyeX, eyeU = I(Np+1), I(Nc);
    eyeU = [eyeU;zeros(Bool,1+Np-Nc,Nc)] # Zeros address that Nc < Np (terminal state)

    tighten_constraints = !iszero(mpc.model.wmin) || !iszero(mpc.model.wmax) || !iszero(mpc.Δx0)

    for c in mpc.constraints 
        mi = size(c.Au,1);
        ks = [k for k in c.ks if k<= Np]
        Ni = length(ks);

        Ax = c.Ax-c.Au*mpc.K
        Ar = isempty(c.Ar) || nrx == 0 ? zeros(mi,nrx) : c.Ar
        Ad = isempty(c.Ad) ? zeros(mi,mpc.model.nd) : c.Ad
        Aup = isempty(c.Aup) ? zeros(mi,mpc.nuprev) : c.Auip
        Ah = iszero(mpc.model.offset) ? zeros(mi,0) : zeros(mi,1)

        Autot = [Autot; kron(eyeU[ks,:],c.Au)]
        Axtot = [Axtot; kron(eyeX[ks,:],[Ax Ar Ad Aup Ah])]

        ubtot = [ubtot;repeat(c.ub,Ni,1)]
        lbtot = [lbtot;repeat(c.lb,Ni,1)]

        if(tighten_constraints)
            FK = mpc.model.F-mpc.model.G*mpc.K
            ut,lt= constraint_tightening(Ax,FK,ks,mpc.model.wmin,mpc.model.wmax,mpc.Δx0)
            ubtot -= ut
            lbtot += lt
        end

        issoft = [issoft;repeat([c.soft],mi*Ni)]
        isbinary = [isbinary;repeat([c.binary],mi*Ni)]
        prios = [prios;repeat([c.prio],mi*Ni)]
    end
    A=Axtot*Γ+Autot;
    W = -Axtot*Φ;

    # Add extra row due to reference preview
    if mpc.settings.reference_tracking && mpc.settings.reference_preview
        Wr = zeros(0,mpc.nr)
        if mpc.settings.reference_condensation
            for c in mpc.constraints
                mi = size(c.Au,1);
                ks = [k for k in c.ks if k<= Np]

                Ar = isempty(c.Ar) ? zeros(mi,mpc.nr) : c.Ar
                Wr = [Wr;repeat(-Ar,length(ks))] 
            end
        else
            eye_r = I(mpc.Np)
            for c in mpc.constraints 
                mi,Ni = size(c.Au,1),sum(c.ks .<=Np);
                if isempty(c.Ar)
                    Wrn = zeros(mi*Ni,mpc.nr)
                else
                    ks = [k-1 for k in c.ks if k<= Np && k>=1] # first ref is at k=2, not k=1 
                    Wrn = [zeros(mi*(Ni-length(ks)),mpc.nr); 
                           kron(eye_r[ks,:],-c.Ar) zeros(mi*length(ks),mpc.model.ny*(Ni-length(ks)))]
                end
                Wr = [Wr; Wrn] 
            end
        end
        W = [W[:,1:nx] Wr W[:,nx+1:end]]
    end

    # Append zero columns for linear cost parameters (they don't affect constraints)
    nl > 0 && (W = [W zeros(size(W,1), nl)])

    return A,ubtot,lbtot,W,issoft,isbinary,prios
end

# Compute A,b,W such that the constraints for a given MPC structure
# are on the form A U<=b W th
function create_constraints(mpc,Φ,Γ)

    n = size(Γ,2);
    nth = sum(get_parameter_dims(mpc))
    !iszero(mpc.model.offset) && (nth +=1); # constant offset in dynamics
    # Control bounds
    if(!isempty(mpc.umax))
        A,bu,bl,W = create_controlbounds(mpc,Γ,Φ)
        issoft= falses(n);
        prios = zeros(Int,n)
        isbinary_single = falses(mpc.model.nu) 
        isbinary_single[mpc.binary_controls] .= true;
        isbinary = repeat(isbinary_single,mpc.Nc)
        if mpc.Nc_binary >= 0
            isbinary[mpc.Nc_binary*mpc.model.nu+1:end].=false
        end
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
    nxp, nrp, ndp, nup, nlp = get_parameter_dims(mpc)

    # ==== From u' R u ====
    H = kron(I(Nc),R);
    f = zeros(size(H,1),1);
    H[end-nu+1:end,end-nu+1:end] .+= (N-Nc)*R # To accound for Nc < N...

    # Compensate for nonzero operating point
    if(!mpc.settings.reference_tracking && !iszero(mpc.model.uo))
        Uo = repeat(mpc.model.uo,Nc)
        f-=H*Uo
        if(!iszero(mpc.K) && !iszero(R)) # contribution from prestabilizing feedback
            KR = [-mpc.K'*R;zeros(nx-size(mpc.K,2),nu)]
            KRtot = [kron(I(Nc),KR);zeros((N-Nc+1)*nx,Nc*nu)]
            KRtot[Nc*nx+1:N*nx,end-nu+1:end] = repeat(KR,N-Nc,1) # Due to control horizon
            GKR= Γ'*KRtot
            f-=(GKR+GKR')*Uo
        end
    end

    # ==== From (Cx)'Q(Cx) ====
    CQCtot  = kron(I(N),Cp'*Q*Cp);
    CQCf = Cf'*Qf*Cf + cat(mpc.weights.Qfx,zeros(nx-nxp,nx-nxp), dims=(1,2))
    CQCtot = cat(CQCtot,CQCf,dims=(1,2))

    H += Γ'*CQCtot*Γ; 
    # f_theta & H_theta for state parameters
    f_theta  = Γ'*CQCtot*Φ; # from x0
    H_theta  = Φ'*CQCtot*Φ
    if(!mpc.settings.reference_tracking && !iszero(mpc.model.xo))
        f -= Γ'*CQCtot*repeat([mpc.model.xo;zeros(nx-nxp)],N+1)
    end

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
                    if(isempty(mpc.settings.traj2setpoint))
                        W = cat(1e6*I(nu),I(nu*(Nc-1)),dims=(1,2))
                        mpc.traj2setpoint = (W*inv(H)*Fr*Is)\(W*inv(H)*Fr)
                    else
                        mpc.traj2setpoint = mpc.settings.traj2setpoint
                    end
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

    # ==== Linear cost on control inputs ====
    if mpc.settings.linear_cost && nlp > 0
        # Linear cost l'U adds directly to the linear term
        # f_theta gets an identity block for the linear cost parameters
        # This maps l = [l_0; l_1; ...; l_{Nc-1}] directly to U
        Fl = I(nlp)
        f_theta = [f_theta Fl]

        # Extend H_theta with zero blocks (linear cost doesn't contribute to quadratic terms)
        nth_current = size(H_theta, 1)
        H_theta = [H_theta          zeros(nth_current, nlp);
                   zeros(nlp, nth_current)  zeros(nlp, nlp)]
    end

    # Add regularization for binary variables (won't change the solution)
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
function mpc2mpqp(mpc::MPC)
    if(mpc.settings.reference_tracking && !iszero(mpc.model.uo) && !iszero(mpc.weights.R))
        @warn("Tracking and a direct penalty on u can yield steady-state error. 
Consider instead to:
1) Setting R=0 and Rr≠0               →  penalize change in u instead
2) Setting reference_tracking = false →  regulation to the operating point.")
    end

    F,G,C,Q,R,S,Qf = create_extended_system_and_cost(mpc) 

    Φ,Γ=state_predictor(F,G,mpc.Np,mpc.Nc);

    # Create objective function 
    nxe,nue= size(G) 
    H,f,f_theta,H_theta = objective(Φ,Γ,C,
                                    Q,R,S,Qf,
                                    mpc.Np,mpc.Nc,nue,nxe,mpc)
    # Create Constraints 
    A, bu, bl, W, issoft, isbinary, prio = create_constraints(mpc,Φ,Γ)

    # Collapse constant term in dynamics
    if(!iszero(mpc.model.offset))
        f += f_theta[:,end]
        f_theta = f_theta[:,1:end-1]
        H_theta = H_theta[1:end-1,1:end-1]
        bu += W[:,end]
        bl += W[:,end]
        W = W[:,1:end-1]
    end

    # Apply move blocking
    if(!isempty(mpc.move_blocks))
        keep_bounds = Int[] 
        nu = mpc.model.nu
        if(!mpc.settings.move_block_foh)
            T= zeros(0,0)
            for (k,mb) in enumerate([mpc.move_blocks[1:end-1];1])
                T = cat(T,repeat(I(nu),mb,1),dims=(1,2))
                keep_bounds = keep_bounds ∪ collect((k-1)*nu+1:k*nu)
            end
        else
            T = zeros(length(f),length(mpc.move_blocks)*nu)
            offset = 0
            for (k,mb) in enumerate(mpc.move_blocks[1:end-1])
                block = [[(mb-i)/mb for i in 0:mb-1] [i/mb for i in 0:mb-1]]
                T[offset*nu+1:(offset+mb)*nu,(k-1)*nu+1:(k+1)*nu] = kron(block,I(nu))
                keep_bounds = keep_bounds ∪ collect((k-1)*nu+1:k*nu)
                offset += mb
            end
            T[end-nu+1:end,end-nu+1:end] = I(nu)
        end
        A, H, f, f_theta = A*T, T'*H*T, T'*f, T'*f_theta

        # remove superfluous control bounds
        keep = keep_bounds ∪ collect(nu*mpc.Nc+1:length(bu))
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

    # Resort based on priorities
    ns = length(prio)-size(A,1)
    prio_order = sortperm(prio[ns+1:end])
    A = A[prio_order,:]
    prio_order = [1:ns; prio_order .+ ns] # Offset correctly
    bu,bl,W = bu[prio_order], bl[prio_order], W[prio_order,:]
    issoft,isbinary,prio = issoft[prio_order],isbinary[prio_order], prio[prio_order]
    break_points = unique(i->prio[i], eachindex(prio))[2:end].-1;
    break_points = Cint.(break_points)
    isempty(break_points) || push!(break_points,length(prio))
                                                
    # Mark soft constrints  
    senses[issoft[:]].+=DAQP.SOFT

    # Mark binary constraints
    senses[isbinary[:]].+=DAQP.BINARY

    # Replace Infs for codegen
    clamp!(bu,-1e30,1e30)
    clamp!(bl,-1e30,1e30)

    mpQP = MPQP(H,f[:],H_theta,f_theta,
                A,bu,bl,W,senses,Cint.(prio),break_points,
                any(isbinary))

    return mpQP
end

function create_extended_system_and_cost(mpc::MPC)
    F,G,C = mpc.model.F-mpc.model.G*mpc.K, mpc.model.G, mpc.model.C
    Q,R,Rr,S = copy(mpc.weights.Q), copy(mpc.weights.R), copy(mpc.weights.Rr), copy(mpc.weights.S)
    Qf = iszero(mpc.weights.Qf) && iszero(mpc.weights.Qfx) ? Q : copy(mpc.weights.Qf)

    nx,nr,nd,nuprev,nl = get_parameter_dims(mpc)
    mpc.nr, mpc.nuprev, mpc.nl = nr, nuprev, nl
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

    if(!iszero(mpc.model.offset))
        F = cat(F,1,dims=(1,2))
        F[1:nx,end] .= mpc.model.offset
        G = [G;zeros(1,nu)]
        S = [S;zeros(1,nu)]
        C = [C zeros(size(C,1),1)]
    end

    return F,G,C,Q,R,S,Qf 
end
