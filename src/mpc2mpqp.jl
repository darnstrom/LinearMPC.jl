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
function create_controlbounds(lb,ub,N,Ncc,nx; double_sided=false)
    Ncc= min(Ncc,N); # Cannot have constraint over N...
    nu = length(lb);
    A = kron(I(Ncc),I(nu)); 
    if(double_sided)
        ubs = repeat(ub,Ncc,1)
        lbs = repeat(lb,Ncc,1)
        return A,ubs,lbs, zeros(length(ubs),nx) 
    else
        A = [A;-A];
        b= [repeat(ub,Ncc,1);-repeat(lb,Ncc,1)];
        W = zeros(length(b),nx); 
        return A,b,W
    end
end

# Create A u <= b+W*theta 
# correspoding to  lby<=Cy xi<=uby for i ∈ {1,2,...,Ncy}
function create_statebounds(Cy,lby,uby,Γ,Φ,N,yb_inds;double_sided=false)
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
    if(double_sided)
        return A,ubtot,lbtot,-Ctot*Φ
    else
        A=[A;-A] # TODO soft constraints
        b= [ubtot;-lbtot];
        W = [-Ctot;Ctot]*Φ;
        return A,b,W
    end
end


# Create A u <= b+W*theta 
# correspoding to Au * u_i+ Ax *x_i <=bg for i ∈ {1,2,...,Ncg} 
# (and Au  u_i <=bg for i = 0)
function create_generalconstraints(Au,Ax,bg,Γ,Φ,N,Ncg;double_sided=false)
    Ncg= min(Ncg,N); # Cannot have constraint over N...
    m,nx = size(Ax);
    nu = size(Au,2);
    Axtot = [zeros(m*Ncg,nx) kron(I(Ncg),Ax) zeros(m*Ncg,nx*(N-Ncg))];
    Autot = [kron(I(Ncg),Au) zeros(m*Ncg,nu*(N-Ncg))];
    b= repeat(bg,Ncg,1);
    A=Axtot*Γ+Autot;
    W = - Axtot*Φ;
    if(double_sided)
        return A,b,repeat(-Inf,length(b)),W
    else
        return A,b,W
    end
end

# Compute A,b,W such that the constraints for a given MPC structure
# are on the form A U<=b W th
function create_constraints(mpc,Φ,Γ;double_sided=false)
    c=mpc.constraints
    return create_constraints(c.lb,c.ub,c.Ncc,c.Cy,c.lby,c.uby,c.Ncy,c.Au,c.Ax,c.bg,c.Ncg,Γ,Φ,mpc.Np,mpc.nx,c.double_sided)
end

function create_constraints(lb,ub,Ncc,Cy,lby,uby,Ncy,Au,Ax,bg,Ncg,Γ,Φ,N,nx,double_sided)
    n = size(Γ,2);
    A = zeros(0,n);
    if(double_sided)
        bu = zeros(0);
        bl = zeros(0);
    else
        b = zeros(0);
    end
    W = zeros(0,nx);
    bounds_table=zeros(Int64,0);
    issoft = falses(0)
    # Control bounds
    if(!isempty(ub))
        if(double_sided)
            ~,bu,bl,W = create_controlbounds(lb,ub,N,Ncc,nx;double_sided=true)
            issoft= falses(n);
        else
            A,b,W = create_controlbounds(lb,ub,N,Ncc,nx);
            bounds_table= [collect(n+1:2*n); collect(1:n)];
            issoft= falses(2*n);
        end
    end

    # State bounds
    if(!isempty(lby))
        m = size(A,1);
        if(double_sided)
            Ay,buy,bly,Wy = create_statebounds(Cy,lby,uby,Γ,Φ,N,Ncy;double_sided=true);
            my = Int(size(Ay,1));
            issoft = [issoft; trues(my)];
            bu = [bu;buy];
            bl = [bl;bly];
        else
            Ay,by,Wy = create_statebounds(Cy,lby,uby,Γ,Φ,N,Ncy);
            b = [b;by];
            my = Int(size(Ay,1)/2);
            bounds_table = [bounds_table; m.+collect(my+1:2*my); m.+collect(1:my)];
            issoft = [issoft; trues(2*my)];
        end
        A = [A;Ay];
        W = [W;Wy];
    end

    # General constraints 
    if(!isempty(bg))
        m = size(A,1);
        mg = Ncg*length(bg);
        if(double_sided)
            Ag,bug,blg,Wg = create_generalconstraints(Au,Ax,bg,Γ,Φ,N,Ncg);
            bu = [bu;bug];
            bl = [bl;blg];
            bounds_table = [bounds_table; m.+collect(1:mg)];
            issoft = [issoft; repeat(sum(abs.(Ax),dims=2).>0,Ncg)];
        else
            Ag,bg,Wg = create_generalconstraints(Au,Ax,bg,Γ,Φ,N,Ncg);
            b = [b;bg];
            bounds_table = [bounds_table; m.+collect(1:mg)];
            issoft = [issoft; repeat(sum(abs.(Ax),dims=2).>0,Ncg)];
        end
        A = [A;Ag];
        W = [W;Wg];
    end
    if(double_sided)
        return A,bu,bl,W,bounds_table,issoft
    else
        return A,b,W,bounds_table,issoft
    end
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
    return objective(Φ,Γ,mpc.C,mpc.weights.Q,mpc.weights.R,mpc.weights.Rr,mpc.Np,mpc.Nc,mpc.nu)
end
function objective(Φ,Γ,C,Q,R,Rr,N,Nc,nu)
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
function mpc2mpqp(mpc::MPC;explicit_soft=true)
    mpQP = LinearMPC.MPQP();

    Φ,Γ=state_predictor(mpc.F,mpc.G,mpc.Np,mpc.Nc);

    # Create objective function 
    H, f_theta, H_theta = objective(mpc,Φ,Γ)
    # Create Constraints 
    if(mpc.constraints.double_sided)
        A, bu, bl, W, bounds_table, issoft = create_constraints(mpc,Φ,Γ)
    else
        A, b, W, bounds_table, issoft = create_constraints(mpc,Φ,Γ)
    end

    W = [W zeros(size(W,1),size(f_theta,2)-mpc.nx)]; # Add zeros for r and u-

    if(explicit_soft && any(issoft))
        A = [A zeros(size(A,1))];
        A[issoft[size(H,1)+1:end],end].=-1;

        H = cat(H,mpc.weights.rho,dims=(1,2));
        f_theta =[f_theta;zeros(1,size(f_theta,2))];
    end
    f = zeros(size(H,1),1); 
    senses = zeros(Cint,size(W,1)); 
    if(mpc.constraints.double_sided)
        return (H=H,f=f, H_theta = H_theta, f_theta=f_theta,
                A=Matrix{Float64}(A), bu=bu, bl=bl, W=W, senses=senses, bounds_table=bounds_table)
    else
        return (H=H,f=f, H_theta = H_theta, f_theta=f_theta,
                A=Matrix{Float64}(A), b=b, W=W, senses=senses, bounds_table=bounds_table)
    end
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

