function invariant_set(F::Matrix{<:Real},xmin=zeros(0),xmax=zeros(0);
              wmin=zeros(0),wmax=zeros(0),
              G=zeros(0),umin=zeros(0),umax=zeros(0),
              Ax=zeros(0),bx=zeros(0),
              Agu=zeros(0,0),Agx=zeros(0,0),bg=zeros(0),
              max_iter=500, eps_shrink=1e-3)

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
    return H,h
end
