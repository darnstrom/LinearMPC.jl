struct DenseObjective
    H::Matrix{Float64}
    f::Vector{Float64}
    f_theta::Matrix{Float64}
    H_theta::Matrix{Float64}
end

struct DenseConstraints
    A::Matrix{Float64}
    bu::Vector{Float64}
    bl::Vector{Float64}
    W::Matrix{Float64}
    issoft::BitVector
    isbinary::BitVector
    prio::Vector{Cint}
end

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
    return Φ::Matrix{Float64},Γ::Matrix{Float64}
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
    nuprev = !iszero(mpc.weights.Rr) || any(!iszero(c.Aup) for c in mpc.constraints) ? mpc.model.nu : 0
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
    nxe = sum(get_parameter_dims(mpc))-nl
    mpc.settings.reference_preview && (nxe-=mpc.nr) # reference not part of extended state
    !iszero(mpc.model.f_offset) && (nxe+=1); # constant offset in dynamics

    #-K*X + V = -K*(Γ V + Φ x0) + V 
    #         = (I-K*Γ)V -K Φ x0
    if(!iszero(mpc.K))
        K = kron(I(Nb),[mpc.K zeros(nu,nxe-nx)])
        A = (I-K*Γ[1:Nb*nxe, 1:Nb*nu])
        W = K*Φ[1:Nb*nxe,:]
        if mpc.settings.reference_tracking && mpc.settings.reference_preview
            # append rows to W from reference that are not part of extended state
            W = [W[:,1:nx] zeros(size(W,1),mpc.nr) W[:,nx+1:end]]
        end
        nl > 0 && (W = [W zeros(size(W,1), nl)]) # append zero rows for linear cost
    else
        nth =  mpc.settings.reference_preview ? nxe+nl+mpc.nr : nxe+nl
        A = zeros(0,mpc.Nc*nu)
        W = zeros(Nb*nu,nth)
    end
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

    nxe = sum(get_parameter_dims(mpc))-nl
    if mpc.settings.reference_preview
        nxe-=mpc.nr # reference not part of extended state
        nrx = 0
    else
        nrx = mpc.nr
    end
    !iszero(mpc.model.f_offset) && (nxe+=1); # constant offset in dynamics

    ubtot,lbtot = zeros(0,1),zeros(0,1);
    Axtot,Autot = zeros(0,nxe*(Np+1)), zeros(0,nu*Nc);
    issoft,isbinary = falses(0),falses(0)
    prios = zeros(Int,0)

    eyeX, eyeU = I(Np+1), I(Nc);
    eyeU = [eyeU;zeros(Bool,1+Np-Nc,Nc)] # Zeros address that Nc < Np (terminal state)

    tighten_constraints = !iszero(mpc.model.wmin) || !iszero(mpc.model.wmax) || !iszero(mpc.Δx0)

    for c in mpc.constraints 
        mi = size(c.Au,1);
        kmax = iszero(c.Au) ? Np+1 : Np
        ks = [k for k in c.ks if k<= kmax]
        Ni = length(ks);

        Ax = c.Ax-c.Au*mpc.K
        Ar = isempty(c.Ar) || nrx == 0 ? zeros(mi,nrx) : c.Ar
        Ad = isempty(c.Ad) ? zeros(mi,mpc.model.nd) : c.Ad
        Aup = isempty(c.Aup) ? zeros(mi,mpc.nuprev) : c.Aup
        Ah = iszero(mpc.model.f_offset) ? zeros(mi,0) : zeros(mi,1)

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
function create_constraints(mpc::MPC,Φ,Γ)

    n = size(Γ,2);
    nth = sum(get_parameter_dims(mpc))
    !iszero(mpc.model.f_offset) && (nth +=1); # constant offset in dynamics
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

    # Collapse constant term in constraints 
    if(!iszero(mpc.model.f_offset))
        bu += W[:,end]
        bl += W[:,end]
        W = W[:,1:end-1]
    end

    return DenseConstraints(A,bu[:],bl[:],W,issoft,isbinary,prios)
end

# Create H, f_theta, H_theta such that the objective function for a given
# MPC problem is formulated as 0.5 U' H U+th'F_theta' U + 0.5 th' H_theta th
# (where th contains x0, r and u(k-1))
function create_objective(mpc::MPC,Φ,Γ,C,w::MPCWeights,nu::Int,nx::Int)

    Q,R,S,Qf = w.Q, w.R, w.S, w.Qf
    N,Nc = mpc.Np,mpc.Nc
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

    # Collapse constant term in objective
    if(!iszero(mpc.model.f_offset))
        f += f_theta[:,end]
        f_theta = f_theta[:,1:end-1]
        H_theta = H_theta[1:end-1,1:end-1]
    end
    # Contant offset in h enters similar to as reference (can be seen as reference r - h_offset)
    if(nrp > 0 && !iszero(mpc.model.h_offset))
        if mpc.settings.reference_preview && !mpc.settings.reference_condensation
            f .-= f_theta[:,nxp+1:nxp+nrp]*repeat(mpc.model.h_offset,mpc.Np)
        else
            f .-= f_theta[:,nxp+1:nxp+nrp]*mpc.model.h_offset
        end
    end

    return DenseObjective((H+H')/2,f[:],f_theta,H_theta)
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

    F,G,C = create_extended_system(mpc) 
    weights = create_extended_cost(mpc,mpc.weights) 
    nxe,nue= size(G) 

    Φ,Γ=state_predictor(F,G,mpc.Np,mpc.Nc);

    objective  = create_objective(mpc,Φ,Γ,C,weights,nue,nxe)
    constraints = create_constraints(mpc,Φ,Γ)

    if(!isempty(mpc.move_blocks))
        objective,constraints = apply_move_block(mpc,objective,constraints)
    end

    # Resort based on priorities
    constraints = sort_constraints(constraints)
    
    # Remove redundant constraints
    if(mpc.settings.preprocess_mpqp)
        constraints = remove_redundant(constraints)
        constraints = remove_duplicate(constraints)
    end

    return MPQP(objective,constraints)
end

function create_extended_system(mpc::MPC)
    F,G = mpc.model.F-mpc.model.G*mpc.K, mpc.model.G
    C = mpc.model.C
    nx,nr,nd,nuprev,nl = get_parameter_dims(mpc)
    mpc.nr, mpc.nuprev, mpc.nl = nr, nuprev, nl
    nu = mpc.model.nu 

    if(nr > 0) # Reference tracking -> add reference to states
        if !mpc.settings.reference_preview
            F = cat(F,I(mpc.model.ny),dims=(1,2))
            G = [G;zeros(mpc.model.ny,nu)]
            C = [C -I(mpc.model.ny)] 
        end
    end

    if(nd > 0) # add measureable disturbance
        F = cat(F,I(nd),dims=(1,2))
        F[1:nx,end-nd+1:end] .= mpc.model.Gd
        G = [G;zeros(nd,nu)]
        C = [C mpc.model.Dd]
    end

    if(nuprev > 0) # Penalizing Δu -> add uold to states 
        F = cat(F,zeros(nu,nu),dims=(1,2))
        F[end-nu+1:end,1:nx] .= -mpc.K
        G = [G;I(nu)]
        nye,nxe = size(C)
        C = [C zeros(nye,nu); mpc.K zeros(nu,nxe-nx) I(nu)]
    end

    if(!iszero(mpc.weights.R) && !iszero(mpc.K)) # terms from prestabilizing feedback
        C = [C; mpc.K zeros(nu,size(C,2)-nx)]
    end

    if(!iszero(mpc.model.f_offset)) # Add offset to states
        F = cat(F,1,dims=(1,2))
        F[1:nx,end] .= mpc.model.f_offset
        G = [G;zeros(1,nu)]
        C = [C zeros(size(C,1),1)]
    end
    return F::Matrix{Float64},G::Matrix{Float64}, C::Matrix{Float64}
end

function create_extended_cost(mpc::MPC, weights::MPCWeights)
    Q,R,Rr,S = copy(weights.Q), copy(weights.R), copy(weights.Rr), copy(weights.S)
    Qf = iszero(weights.Qf) && iszero(weights.Qfx) ? Q : copy(weights.Qf)
    nx,nr,nd,nuprev,nl = get_parameter_dims(mpc)
    nu = mpc.model.nu 

    if nr > 0 && !mpc.settings.reference_preview # Reference tracking -> add reference to states
        S = [S;zeros(mpc.model.ny,nu)]
    end

    if(nd > 0) # add measureable disturbance 
        S = [S;zeros(nd,nu)]
    end

    if(nuprev > 0) # Penalizing Δu -> add uold to states 
        Q = cat(Q,Rr,dims=(1,2))
        Qf = cat(Qf,zeros(nu,nu),dims=(1,2))
        S = [S;-Rr];
        S[1:nx,:] -=mpc.K'*Rr
        R+=Rr
    end

    if(!iszero(mpc.weights.R) && !iszero(mpc.K)) # terms from prestabilizing feedback
        Q = cat(Q,mpc.weights.R,dims=(1,2))
        Qf = cat(Qf,zeros(nu,nu),dims=(1,2))
        S[1:nx,:] -=mpc.K'*mpc.weights.R
    end

    if(!iszero(mpc.model.f_offset))
        S = [S;zeros(1,nu)]
    end

    return MPCWeights(Q,R,zeros(0,0),S,Qf,zeros(0,0))
end

function remove_redundant(c::DenseConstraints)
    A,bu,bl,W = c.A, c.bu, c.bl, c.W
    issoft,isbinary,prio = c.issoft, c.isbinary, c.prio
    nsimple = length(bu) - size(A,1)
    keep = collect(1:nsimple) # Don't remove simple bounds
    norm_factors = ones(nsimple)
    for (i,a) in enumerate(eachrow(A))
        id = nsimple+i
        norm_factor = norm(a)
        if norm_factor > 1e-10
            nz_id = findfirst(x -> abs(x) > 1e-12, a)
            if sign(a[nz_id]) < 0 # To ensure unique half planes, make sure first nonzero is zero 
                a.=-a.+0.0 #+0.0 to avoid -0.0
                bui = bu[id]
                bu[id] = -bl[id]
                bl[id] = -bui
                W[id,:] .=-W[id,:] .+ 0.0
            end
            if count(x -> abs(x) > 1e-12, a) == 1 # check if bound is a simple bound
                if nz_id <= nsimple  && prio[nz_id] == prio[id]
                    if issoft[nz_id] == issoft[id] && isbinary[nz_id] == isbinary[id]
                        if iszero(W[id,:]-W[nz_id,:])
                            bu[nz_id] = min(bu[nz_id],bu[id]/norm_factor)
                            bl[nz_id] = max(bl[nz_id],bl[id]/norm_factor)
                            continue
                        end
                    end
                end
            end
            push!(keep,nsimple+i) # Only add nonzero rows
            push!(norm_factors,1/norm_factor)
        end
    end
    if length(keep) < length(bu)
        keepA = keep[nsimple+1:end].-nsimple
        A = A[keepA,:].*norm_factors[nsimple+1:end]
        bu,bl,W=bu[keep].*norm_factors,bl[keep].*norm_factors,W[keep,:].*norm_factors
        issoft,isbinary,prio = issoft[keep],isbinary[keep],prio[keep]
    end
    return DenseConstraints(A,bu,bl,W,issoft,isbinary,prio)
end

function find_duplicate_rows(A::AbstractMatrix; digits=6)
    # Group by rounded values
    row_groups = Dict{Vector{Float64}, Vector{Int}}()
    order = Vector{Float64}[]
    for i in 1:size(A, 1)
        # Round the row to handle floating point jitter
        row = round.(A[i, :], digits=digits)

        if haskey(row_groups, row)
            push!(row_groups[row], i)
        else
            row_groups[row] = [i]
            push!(order,row)
        end
    end
    return [row_groups[row] for row in order]
end

function remove_duplicate(c::DenseConstraints)
    A,bu,bl,W = c.A, c.bu, c.bl, c.W
    issoft,isbinary,prio = c.issoft, c.isbinary, c.prio

    nsimple = length(bu) - size(A,1)
    idsA = nsimple+1:length(bu)
    Aext = [A W[idsA,:] issoft[idsA] isbinary[idsA] prio[idsA]]
    duplicate_map = find_duplicate_rows(Aext)

    if(length(duplicate_map) == size(A,1)) # No duplicates
        return c 
    end

    A_new = zeros(length(duplicate_map),size(A,2))
    bu_new = [bu[1:nsimple];zeros(length(duplicate_map))]
    bl_new = [bl[1:nsimple];zeros(length(duplicate_map))]
    W_new = [W[1:nsimple,:];zeros(length(duplicate_map),size(W,2))]
    issoft_new = [issoft[1:nsimple];falses(length(duplicate_map))]
    isbinary_new = [isbinary[1:nsimple];falses(length(duplicate_map))]
    prio_new = [prio[1:nsimple];zeros(Cint,length(duplicate_map))]
    for i in 1:length(duplicate_map)
        # TODO, correctly handle prios
        id =  duplicate_map[i][1]
        A_new[i,:] = A[id,:]
        ids = duplicate_map[i].+nsimple
        id+=nsimple # To go account for simple bounds
        bu_new[nsimple+i]= minimum(bu[ids])
        bl_new[nsimple+i] = maximum(bl[ids])
        W_new[nsimple+i,:] = W[id,:]
        issoft_new[nsimple+i] = issoft[id]
        isbinary_new[nsimple+i] = isbinary[id]
        prio_new[nsimple+i] = prio[id]
    end

    return DenseConstraints(A_new,bu_new,bl_new,W_new,issoft_new,isbinary_new,prio_new)
end

function apply_move_block(mpc::MPC, obj::DenseObjective, c::DenseConstraints)
    nu = mpc.model.nu
    nu_bounds = length(mpc.umax)

    nUold = nu*mpc.Nc 
    nUnew = sum(length(mb) for mb in mpc.move_blocks)

    new_id,T,counter = 1,zeros(nUold,nUnew),collect(1:nu)
    keep = Int[]
    for pass in 1:maximum(length,mpc.move_blocks)
        for (iu,mb) in enumerate(mpc.move_blocks)
            length(mb) < pass  && continue # No more blocks for control iu
            block = length(mb) != pass ? mb[pass] : 1 # clipping since the end will be superfluous...
            T[counter[iu]:nu:counter[iu]+nu*(block-1),new_id] .= 1
            counter[iu] <= nu_bounds*mpc.Nc && append!(keep,counter[iu])
            counter[iu]+=nu*block
            new_id +=1
        end
    end
    new_obj = DenseObjective(T'*obj.H*T, T'*obj.f, T'*obj.f_theta, obj.H_theta)

    append!(keep,nu_bounds*mpc.Nc+1:length(c.bu))

    # Remove superfluous control bounds
    Anew = !iszero(mpc.K) ? c.A[keep,:]*T : c.A*T # control bounds are in A if prestabilizing feedback
    new_c = DenseConstraints(Anew, c.bu[keep], c.bl[keep], c.W[keep,:], c.issoft[keep], c.isbinary[keep], c.prio[keep])
    return new_obj,new_c
end

function sort_constraints(c::DenseConstraints)
    ns = length(c.prio)-size(c.A,1)
    order = sortperm(c.prio[ns+1:end])
    Anew = c.A[order,:]
    order = [1:ns; order .+ ns] # Offset correctly
    return DenseConstraints(Anew, c.bu[order], c.bl[order], c.W[order,:], 
                            c.issoft[order], c.isbinary[order], c.prio[order])
end

function MPQP(obj::DenseObjective, c::DenseConstraints)
    # Setup sense
    senses = zeros(Cint,length(c.bu)); 

    # ignore constraints which have inf bounds
    for i in 1:length(c.bu)
        if(c.bu[i] > 1e20 && c.bl[i] < -1e20)
            senses[i] = DAQP.IMMUTABLE
        elseif abs(c.bu[i]-c.bl[i]) < 1e-12
            senses[i] = DAQP.EQUALITY
        end
    end
    # Mark soft constrints  
    senses[c.issoft[:]].+=DAQP.SOFT

    # Mark binary constraints
    senses[c.isbinary[:]].+=DAQP.BINARY

    # Replace Infs for codegen
    clamp!(c.bu,-1e30,1e30)
    clamp!(c.bl,-1e30,1e30)

    break_points = unique(i->c.prio[i], eachindex(c.prio))[2:end].-1;
    break_points = Cint.(break_points)
    isempty(break_points) || push!(break_points,length(c.prio))

    m,n = length(c.bu),length(obj.f)
    mpQP = MPQP(obj.H,obj.f,obj.H_theta,obj.f_theta,
                c.A,c.bu,c.bl,c.W, senses, c.prio, break_points,
                any(c.isbinary),
                zeros(m),zeros(m),zeros(m),zeros(n))
end
