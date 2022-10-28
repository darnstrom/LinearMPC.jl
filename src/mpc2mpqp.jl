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
        Ftot*=F;
        Gtot=F*Gtot;
    end
    # Set ui = u_Nc for i>Nc 
    # TODO add infinite LQR gain alternative
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
    c = mpc.constraints
    Ncc= min(c.Ncc,mpc.Nc); # Cannot have constraint over N...
    nu,nx = mpc.nu, mpc.nx 
    lb,ub = c.lb, c.ub


    A = kron(I(Ncc),I(nu)); 
    ubs = repeat(ub,Ncc,1)
    lbs = repeat(lb,Ncc,1)
    return A,ubs,lbs, zeros(length(ubs),nx) 
end

# Create A u <= b+W*theta 
# correspoding to  lby<=Cy xi<=uby for i ∈ {1,2,...,Ncy}
function create_statebounds(mpc::MPC,Γ,Φ)
    c = mpc.constraints
    N = mpc.Np
    Cy,lby,uby,yb_inds = c.Cy, c.lby, c.uby, c.Ncy

    m = length(Cy);
    nx = size(Φ,2);
    Ctot = zeros(0,nx*(N+1));
    ubtot = zeros(0,1);
    lbtot = zeros(0,1);
    eye = I(N);

    for i in 1:m
        mi = size(Cy[i],1);
        Ni = length(yb_inds[i]);
        Ctot = [Ctot;
                zeros(Ni*mi,nx) kron(eye[yb_inds[i],:],Cy[i]) ];
        ubtot = [ubtot;repeat(uby[i],Ni,1)];
        lbtot= [lbtot;repeat(lby[i],Ni,1)];
    end
    A=Ctot*Γ
    return A,ubtot,lbtot,-Ctot*Φ
end

# Create A u <= b+W*theta 
# correspoding to Au * u_i+ Ax *x_i <=bg for i ∈ {1,2,...,Ncg} 
# (and Au  u_i <=bg for i = 0)
function create_generalconstraints(mpc::MPC,Γ,Φ)
    # extract data
    c = mpc.constraints
    Au,Ax,bg = c.Au, c.Ax, c.bg
    N = mpc.Nc
    Ncg= min(c.Ncg,N); # Cannot have constraint over N...
    m = length(bg)
    nx,nu = mpc.nx, mpc.nu

    Axtot = [zeros(m*Ncg,nx) kron(I(Ncg),Ax) zeros(m*Ncg,nx*(N-Ncg))];
    Autot = [kron(I(Ncg),Au) zeros(m*Ncg,nu*(N-Ncg))];
    b= repeat(bg,Ncg,1);
    A=Axtot*Γ+Autot;
    W = - Axtot*Φ;

    return A,b,repeat([-1e30],length(b)),W
end

# Compute A,b,W such that the constraints for a given MPC structure
# are on the form A U<=b W th

function create_constraints(mpc,Φ,Γ)
    n = size(Γ,2);
    A = zeros(0,n);
    Ax = mpc.constraints.Ax
    bu = zeros(0);
    bl = zeros(0);
    W = zeros(0,mpc.nx);
    issoft = falses(0)
    isbinary= falses(0)

    # Control bounds
    if(!isempty(mpc.constraints.ub))
        ~,bu,bl,W = create_controlbounds(mpc)
        issoft= falses(n);
        isbinary_single = falses(mpc.nu) 
        isbinary_single[mpc.constraints.binary_controls] .= true;
        isbinary = repeat(isbinary_single,mpc.constraints.Ncc)
    end

    # State bounds
    if(!isempty(mpc.constraints.lby))
        m = size(A,1);
        Ay,buy,bly,Wy = create_statebounds(mpc,Γ,Φ);
        my = Int(size(Ay,1));
        issoft = [issoft; trues(my)];
        isbinary = [isbinary; falses(my)]
        bu = [bu;buy];
        bl = [bl;bly];
        A = [A;Ay];
        W = [W;Wy];
    end

    # General constraints 
    if(!isempty(mpc.constraints.bg))
        m = size(A,1);
        Ncg = mpc.constraints.Ncg
        mg = Ncg*length(mpc.constraints.bg);
        Ag,bug,blg,Wg = create_generalconstraints(mpc,Γ,Φ);
        bu = [bu;bug];
        bl = [bl;blg];
        issoft = [issoft; repeat(sum(abs.(Ax),dims=2).>0,Ncg)];
        isbinary = [isbinary; falses(length(bug))];
        A = [A;Ag];
        W = [W;Wg];
    end
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
    return objective(Φ,Γ,mpc.C,mpc.weights.Q,mpc.weights.R,mpc.weights.Rr,mpc.Np,mpc.Nc,mpc.nu,mpc.weights.Qf,mpc.nx)
end
function objective(Φ,Γ,C,Q,R,Rr,N,Nc,nu,Qf,nx)
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
    f_theta = [f_theta T'*Rrtot *S]; # from u[k-1]

    # Quadratic parameter term (relevant to maintain a positive definite value function)
    Hth_xx=Φ'*Ctot'*Qtot*Ctot*Φ; 
    Hth_rx = Reftot'*Qtot*Ctot*Φ; 
    Hth_rr= Reftot'*Qtot*Reftot
    Hth_uu=  S'*Rrtot*S

    H_theta = cat([Hth_xx Hth_rx';Hth_rx Hth_rr],Hth_uu,dims=(1,2)); 

    return H,f_theta,H_theta
end

"""
    mpQP = mpc2mpqp(mpc)

For a given MPC structure `mpc`, form the multi-parametric QP `mpQP` in the form 
```
min			0.5 U' H U+th'F_theta' U + 0.5 th' H_theta th
subject to 	A U <= b + Wth
```
"""
function mpc2mpqp(mpc::MPC)
    mpQP = LinearMPC.MPQP();

    Φ,Γ=state_predictor(mpc.F,mpc.G,mpc.Np,mpc.Nc);

    # Create objective function 
    H, f_theta, H_theta = objective(mpc,Φ,Γ)
    H = (H+H')/2
    # Create Constraints 
    A, bu, bl, W, issoft, isbinary = create_constraints(mpc,Φ,Γ)

    if(mpc.settings.reference_tracking)
        W = [W zeros(size(W,1),size(f_theta,2)-mpc.nx)]; # Add zeros for r and u-
    else
        f_theta = f_theta[:,1:mpc.nx] # Remove terms for r and u-
    end

    senses = zeros(Cint,size(W,1)); 

    if(mpc.settings.soft_constraints)
        if(mpc.settings.explicit_soft && any(issoft))
            A = [A zeros(size(A,1))];
            A[issoft[size(H,1)+1:end],end].=-1;

            H = cat(H,mpc.weights.rho,dims=(1,2));
            f_theta =[f_theta;zeros(1,size(f_theta,2))];
        end
        if(!mpc.settings.explicit_soft)
            senses[issoft[:]].+=DAQP.SOFT
        end
    end
    f = zeros(size(H,1),1); 
    senses[isbinary[:]].+=DAQP.BINARY
    if(mpc.settings.QP_double_sided)
        mpQP = (H=H,f=f, H_theta = H_theta, f_theta=f_theta,
                A=Matrix{Float64}(A), bu=bu, bl=bl, W=W, senses=senses)
    else # Transform bl + W θ ≤ A U ≤ bu + W θ → A U ≤ b + W
        ncstr = length(bu);
        n_bounds = ncstr-size(A,1);
        bounds_table=[collect(ncstr+1:2*ncstr);collect(1:ncstr)]
        A = [I(n_bounds) zeros(n_bounds,size(A,2)-n_bounds);A]
        A = [A;-A]
        if(mpc.settings.explicit_soft && any(issoft))# Correct sign for slack
            A[:,end].= -abs.(A[:,end])
        end
        b = [bu;-bl]
        W = [W;-W]
        senses = [senses;senses]
        mpQP = (H=H,f=f, H_theta = H_theta, f_theta=f_theta,
                A=Matrix{Float64}(A), b=b, W=W, senses=senses, bounds_table=bounds_table)
    end
    mpc.mpQP = mpQP
    return mpQP
end

function dualize(mpQP,n_control;normalize=true)
    n = size(mpQP.H,1)
    R = cholesky((mpQP.H+mpQP.H')/2)
    M = mpQP.A/R.U
    Vth = (R.L)\mpQP.f_theta
    Dth = mpQP.W + [Matrix{Float64}(I(n))/R.U; M]*Vth
    du = mpQP.bu
    dl = mpQP.bl

    Dth = Dth'[:,:]# Col major...
    if(normalize)
        # Normalize
        norm_factor = 0
        for i in 1:size(M,1)
            norm_factor = norm(M[i,:],2)
            M[i,:]./=norm_factor
            Dth[:,i]./=norm_factor
            du[i]/= norm_factor
            dl[i]/= norm_factor
        end
    end
    Xth = -mpQP.H\mpQP.f_theta;
    Xth = Xth[1:n_control,:];
    Xth = Xth'[:,:] # col major => row major

    return (M=M, Dth=Dth, du=du, dl=dl, Xth=Xth)
end

