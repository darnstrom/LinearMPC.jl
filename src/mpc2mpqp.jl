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

    G = G isa AbstractMatrix ? [G for i in 1:Nc] : G ∪ [G[end] for i in length(G)+1:Nc]
    F = F isa AbstractMatrix ? [F for i in 1:Nc] : F ∪ [F[end] for i in length(F)+1:Nc]

    Gtot = copy(G);
    Ftot = copy(F[1]);
    for i in 1:Nc 
        for j in 0:(Nc-i)
            Γ[((i+j)*nx+1):(i+j+1)*nx,(j*nu+1):(j+1)*nu] = Gtot[j+1]; 
        end
        Φ[i*nx+1:(i+1)*nx,:] = Ftot;
        i == Nc && break
        Ftot*=F[i+1];
        Gtot=[F[i+1]*Gtot[j] for j in 1:(Nc-i)];
    end
    # Set ui = u_Nc for i>Nc 
    # TODO add infinite LQR gain alternative
    F,G = F[end],G[end] # For now, just keep this constant... #TODO do this
    for i in Nc+1:Np
        Γ[(nx*i+1):nx*(i+1),:] .= F*Γ[(nx*(i-1)+1):nx*i,:];
        Γ[(nx*i+1):nx*(i+1),end-nu+1:end] +=G;

        Φ[i*nx+1:(i+1)*nx,:] = F*Φ[nx*(i-1)+1:nx*i,:]
    end
    return Φ,Γ
end

# Create A u <= b+W*theta 
# correspoding to  lb<=u_i<=ub for i ∈ {1,2,...,Nc}
function create_controlbounds(mpc::MPC)
    nu,nx = mpc.nu, mpc.nx

    if !isempty(mpc.K) || !mpc.settings.QP_double_sided
        A = Matrix{Float64}(kron(I(mpc.Nb),I(nu)))
    else
        A = nothing;
    end
    ub = repeat(mpc.umax,mpc.Nb,1)
    lb = repeat(mpc.umin,mpc.Nb,1)
    return A,ub,lb, zeros(length(ub),nx) 
end

# Create A u <= b+W*theta 
# correspoding to lb <= Au*uk+ Ax*xk <=ub for k ∈ ks 
function create_general_constraints(mpc::MPC,Γ,Φ)
    # extract data
    Np, Nc= mpc.Np, mpc.Nc
    m= length(mpc.constraints)
    nu,nx = mpc.nu, mpc.nx 

    ubtot,lbtot = zeros(0,1),zeros(0,1);
    Axtot,Autot = zeros(0,nx*(Np+1)), zeros(0,nu*Nc);
    issoft,isbinary = falses(0),falses(0)

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
    end
    A=Axtot*Γ+Autot;
    W = -Axtot*Φ;

    return A,ubtot,lbtot,W,issoft,isbinary
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
        Ac,bu,bl,W = create_controlbounds(mpc)
        if !isnothing(Ac) A = Ac end
        issoft= falses(n);
        isbinary_single = falses(mpc.nu) 
        isbinary_single[mpc.binary_controls] .= true;
        isbinary = repeat(isbinary_single,mpc.Nb)
    end

    # General constraints
    if(!isempty(mpc.constraints))
        Ag,bug,blg,Wg,softg,binaryg = create_general_constraints(mpc,Γ,Φ);
        my = Int(size(Ag,1));
        issoft = [issoft; softg]; # Start with sofetning all general constraints
        isbinary = [isbinary; binaryg]
        bu = [bu;bug];
        bl = [bl;blg];
        A = [A;Ag];
        W = [W;Wg];
    end
    # TODO remove inf bounds...

    return A,bu,bl,W,issoft,isbinary
end


# Return T and S such that
# Δu = T*u+S u_minus
# u = inv(T)*Δu - inv(T) S u_minus
# Δu' D Δu = u'T'*D*Tu+2*(T' D S u_minus)' u + u_minus S' D S u_minus
# I.e. : Hd = T'DT, f_theta = 2T'DS, Hth = S'DS
function rate_to_abs(nu,N)
    T = diagm(0=>ones(nu*N),-nu=>-ones(nu*(N-1))); 
    S = [-I(nu);zeros(nu*(N-1),nu)];
    return T,S
end

# Create H, f_theta, H_theta such that the objective function for a given
# MPC problem is formulated as 0.5 U' H U+th'F_theta' U + 0.5 th' H_theta th
# (where th contains x0, r and u(k-1))
function objective(mpc,Φ,Γ)
    return objective(Φ,Γ,mpc.C,mpc.weights.Q,mpc.weights.R,mpc.weights.Rr,mpc.Np,mpc.Nc,mpc.nu,mpc.weights.Qf,mpc.nx,mpc)
end
function objective(Φ,Γ,C,Q,R,Rr,N,Nc,nu,Qf,nx,mpc)
    pos_ids= findall(diag(Q).>0); # Ignore zero indices... (and negative)
    Q = Q[pos_ids,pos_ids];
    C = C[pos_ids,:];
    nr = size(Q,1);
    T,S= rate_to_abs(nu,Nc);

    # Quadratic term 
    H = kron(I(Nc),R);  # from u'R u
    Rrtot = kron(I(Nc),Rr)
    H += T'*Rrtot*T; # from Δu'Rr Δu

    Qtot = kron(I(N+1),Q);
    if(!isnothing(Qf))
        Qtot[end-nx+1:end, end-nx+1:end] .= Qf
    end
    Ctot  = kron(I(N+1),C);

    H += Γ'*Ctot'*Qtot*Ctot*Γ; # from (Cx-r)Q(Cx-r)

    Reftot = repeat(I(nr),N+1); # Assume r constant over the horizon

    # Linear term
    f_theta = Γ'*Ctot'*Qtot*Ctot*Φ; # from x0
    f_theta = [f_theta -Γ'*Ctot'*Qtot*Reftot]; # from r
    f_theta = [f_theta T'*Rrtot*S]; # from u[k-1]

    # Quadratic parameter term 
    # (relevant to maintain a positive definite value function)
    Hth_xx=Φ'*Ctot'*Qtot*Ctot*Φ; 
    Hth_rx = Reftot'*Qtot*Ctot*Φ; 
    Hth_rr= Reftot'*Qtot*Reftot
    Hth_uu=  S'*Rrtot*S

    H_theta = cat([Hth_xx Hth_rx';Hth_rx Hth_rr],Hth_uu,dims=(1,2)); 

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
    Φ,Γ=state_predictor(mpc.F,mpc.G,mpc.Np,mpc.Nc);

    # Create objective function 
    H,f, f_theta, H_theta = objective(mpc,Φ,Γ)
    H = (H+H')/2
    # Create Constraints 
    A, bu, bl, W, issoft, isbinary = create_constraints(mpc,Φ,Γ)

    if(mpc.settings.reference_tracking)
        W = [W zeros(size(W,1),size(f_theta,2)-mpc.nx)]; # Add zeros for r and u-
    else
        f_theta = f_theta[:,1:mpc.nx] # Remove terms for r and u-
    end

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
                A=Matrix{Float64}(A), bu=bu, bl=bl, W=W, senses=senses)
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
        mpQP = (H=H,f=f, H_theta = H_theta, f_theta=f_theta,
                A=Matrix{Float64}(A), b=b, W=W, senses=senses,
                bounds_table=bounds_table)
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
