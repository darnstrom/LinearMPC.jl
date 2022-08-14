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
function create_controlbounds(lb,ub,N,Ncc,nx)
    Ncc= min(Ncc,N); # Cannot have constraint over N...
    nu = length(lb);
    A = kron(I(Ncc),I(nu)); 
    A = [A;-A];
    b= [repeat(ub,Ncc,1);-repeat(lb,Ncc,1)];
    W = zeros(length(b),nx); 
    return A,b,W
end

# Create A u <= b+W*theta 
# correspoding to  lby<=Cy xi<=uby for i ∈ {1,2,...,Ncy}
function create_statebounds(Cy,lby,uby,Γ,Φ,N,yb_inds)
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
    b= [ubtot;-lbtot];
    A=Ctot*Γ
    A=[A;-A] # TODO soft constraints
    W = [-Ctot;Ctot]*Φ;
    return A,b,W
end


# Create A u <= b+W*theta 
# correspoding to Au * u_i+ Ax *x_i <=bg for i ∈ {1,2,...,Ncg} 
# (and Au  u_i <=bg for i = 0)
function create_generalconstraints(Au,Ax,bg,Γ,Φ,N,Ncg)
    Ncg= min(Ncg,N); # Cannot have constraint over N...
    m,nx = size(Ax);
    nu = size(Au,2);
    Axtot = [zeros(m*Ncg,nx) kron(I(Ncg),Ax) zeros(m*Ncg,nx*(N-Ncg))];
    Autot = [kron(I(Ncg),Au) zeros(m*Ncg,nu*(N-Ncg))];
    b= repeat(bg,Ncg,1);
    A=Axtot*Γ+Autot;
    W = - Axtot*Φ;
    return A,b,W
end

# Compute A,b,W such that the constraints for a given MPC structure
# are on the form A U<=b W th
function create_constraints(mpc,Φ,Γ)
    c=mpc.constraints
    return create_constraints(c.lb,c.ub,c.Ncc,c.Cy,c.lby,c.uby,c.Ncy,c.Au,c.Ax,c.bg,c.Ncg,Γ,Φ,mpc.Np,mpc.nx)
end

function create_constraints(lb,ub,Ncc,Cy,lby,uby,Ncy,Au,Ax,bg,Ncg,Γ,Φ,N,nx)
    n = size(Γ,2);
    A = zeros(0,n);
    b = zeros(0);
    W = zeros(0,nx);
    bounds_table=zeros(Int64,0);
    issoft = falses(0)
    # Control bounds
    if(!isempty(ub))
        A,b,W = create_controlbounds(lb,ub,N,Ncc,nx);
        bounds_table= [collect(n+1:2*n); collect(1:n)];
        issoft= falses(2*n);
    end

    # State bounds
    if(!isempty(lby))
        m = size(A,1);
        Ay,by,Wy = create_statebounds(Cy,lby,uby,Γ,Φ,N,Ncy);
        A = [A;Ay];
        b = [b;by];
        W = [W;Wy];
        my = Int(size(Ay,1)/2);
        bounds_table = [bounds_table; m.+collect(my+1:2*my); m.+collect(1:my)];
        issoft = [issoft; trues(2*my)];
    end

    # General constraints 
    if(!isempty(bg))
        m = size(A,1);
        mg = Ncg*length(bg);
        Ag,bg,Wg = create_generalconstraints(Au,Ax,bg,Γ,Φ,N,Ncg);
        A = [A;Ag];
        b = [b;bg];
        W = [W;Wg];
        bounds_table = [bounds_table; m.+collect(1:mg)];
        issoft = [issoft; repeat(sum(abs.(Ax),dims=2).>0,Ncg)];
    end
    return A,b,W,bounds_table,issoft
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
function mpc2mpqp(mpc::MPC)
    mpQP = LinearMPC.MPQP();

    Φ,Γ=state_predictor(mpc.F,mpc.G,mpc.Np,mpc.Nc);

    # Create objective function 
    mpQP.H, mpQP.f_theta, mpQP.H_theta = objective(mpc,Φ,Γ)
    # Create Constraints 
    mpQP.A, mpQP.b,mpQP.W, mpQP.bounds_table, issoft = create_constraints(mpc,Φ,Γ)
    mpQP.W = [mpQP.W zeros(length(mpQP.b),size(mpQP.f_theta,2)-mpc.nx)]; # Add zeros for r and u-

    if(any(issoft))
        mpQP.A = [mpQP.A zeros(size(mpQP.A,1))];
        mpQP.A[issoft,end].=-1;

        mpQP.H = cat(mpQP.H,mpc.weights.rho,dims=(1,2));
        mpQP.f_theta =[mpQP.f_theta;zeros(1,size(mpQP.f_theta,2))];
    end
    mpQP.f = zeros(size(mpQP.H,1),1); 
    mpQP.senses = zeros(Cint,length(mpQP.b)); 
    return mpQP
end
