function cis(mpc::MPC;xmin=zeros(0),xmax=zeros(0),wmin=zeros(0),wmax=zeros(0),max_iter=500,eps_shrink=1e-3,verbose=true)
    F,G = mpc.model.F,mpc.model.G
    nx,nu = size(G)
    umin,umax = mpc.umin,mpc.umax
    Ax,bx = zeros(0,nx),zeros(0)
    Agu,Agx,bg = zeros(0,nu),zeros(0,nx),zeros(0)
    #  Extract constraint into general form
    # TODO : account for special constraints and ks
    # TODO : account for Ar,Aw,Ad,Aup...
    for c in mpc.constraints
        if(!isempty(c.Ax))
            if(!isempty(c.Au))
                Agu,Agx,bg = [Agu;c.Au;-c.Au],[Agx;c.Ax;-c.Ax],[bg;c.ub;-c.lb]
            else
                Ax,bx = [Ax;c.Ax;-c.Ax],[bx;c.ub;-c.lb]
            end
        else
            if(!isempty(c.Au))
                Agu,Agx,bg = [Agu;c.Au;-c.Au],[Agx;zeros(2*length(c.ub),nx)],[bg;c.ub;-c.lb]
            end
        end
    end
    # TODO prune inf bounds
    invariant_set(F,xmin,xmax;wmax,wmin,G,umin,umax,Ax,bx,Agu,Agx,bg,max_iter,eps_shrink)
end
function invariant_set(mpc::MPC,K=zeros(0,0);xmin=zeros(0),xmax=zeros(0),wmin=zeros(0),wmax=zeros(0),max_iter=500,eps_shrink=1e-3,verbose=true)
    F,G = mpc.model.F,mpc.model.G
    if(isempty(K))
        K = -mpc.K
    end
    nx,nu = size(G)
    umin,umax = mpc.umin,mpc.umax
    Ax,bx = zeros(0,nx),zeros(0)
    #  Extract constraint into general form
    # TODO : account for special constraints and ks
    # TODO : account for Ar,Aw,Ad,Aup...
    for c in mpc.constraints
        if(!isempty(c.Ax))
            if(!isempty(c.Au))
                Ax,bx = [Ax;c.Ax+c.Au*K;-c.Ax-c.Au*K],[bx;c.ub;-c.lb]
            else
                Ax,bx = [Ax;c.Ax;-c.Ax],[bx;c.ub;-c.lb]
            end
        else
            if(!isempty(c.Au))
                Ax,bx = [Ax;c.Au*K;-c.Au*K],[bx;c.ub;-c.lb]
            end
        end
    end
    if(!isempty(umin))
        Ax,bx = [Ax;I(nu)*K;-I(nu)*K],[bx;umax;-umin]
    end
    # TODO prune inf bounds
    invariant_set(F,xmin,xmax;wmax,wmin,Ax,bx,max_iter,eps_shrink)
end
function invariant_set(F::Matrix{<:Real},xmin=zeros(0),xmax=zeros(0);
              wmin=zeros(0),wmax=zeros(0),
              G=zeros(0),umin=zeros(0),umax=zeros(0),
              Ax=zeros(0),bx=zeros(0),
              Agu=zeros(0,0),Agx=zeros(0,0),bg=zeros(0),
              max_iter=500, eps_shrink=1e-3, verbose=true)

    # Setup the sets X = {x : H*x ≤ h }, U = {u : Hu*u ≤ hu}
    nx,nu = size(F,1),size(G,2)
    if(!isempty(xmax))
        H = [I(nx) -I(nx)]
        h = [xmax;-xmin]
    else
        H,h = zeros(nx,0),zeros(0)
    end
    if(!isempty(Ax))
        H,h = [H Ax'], [h;bx]
    end

    if(!isempty(umax))
        Hgu = [I(nu) -I(nu)]
        hg = [umax;-umin]
        Hgx = zeros(nx,2*nu)
    else
        Hgu,hg,Hgx = zeros(nu,0),zeros(0),zeros(nx,0)
    end

    if(!isempty(Agu))
        Hgu = [Hgu Agu']
        hg = [hg;bg]
        if(!isempty(Agx))
            Hgx = [Hgx Agx']
        end
    end

    # Start the iterations
    for iter in 1:max_iter
        verbose && println("\r>> #$iter | Constraints: $(length(h))        ");
        hadd = copy(h)
        if(!isempty(wmax))
            # compute hi - max_{w∈W} hi'*w (where W is a box => explicit solution)
            for (i,hi) in enumerate(eachcol(H))
                hadd[i] -= sum(hi[j] * ( hi[j] > 0 ? wmax[j] : wmin[j]) for j in 1:nx)
            end
        end

        if(nu > 0 && !isempty(Hgx))
            Ht = [F'*H Hgx; G'*H Hgu]
            Hadd,hadd = eliminate(Ht,[hadd;hg],collect(nx+1:nx+nu))
        else
            Hadd = F'*H
            normalize!(Hadd,hadd)
        end
        nold  = length(h)
        H,h = minrep([H Hadd], [h;hadd]; keep=1:nold, tol_weak=1e-6+1e-5)
        length(h) == nold && break # All new constraints redundant -> terminate

        h[nold+1:end].-=eps_shrink # shrink new constraints
        H,h = minrep(H, h; keep=nold+1:length(h))
    end
    verbose && println("")
    return H,h
end
