function zoh(A,B,Ts)
  dims = size(B);
  if(length(dims)==1)
	nx = dims[1]
	nu = 1
  else
	nx = dims[1]
	nu = dims[2]
  end
  M = exp([A*Ts  B*Ts;
		   zeros(nu, nx + nu)])
  F = M[1:nx, 1:nx]
  G = M[1:nx, nx+1:nx+nu]
  return F,G
end

function compute_control(mpc::MPC,θ)
    mpQP = mpc.mpQP
    if(mpc.settings.QP_double_sided)
        bth = mpQP.W*θ
        bu = mpQP.bu + bth 
        bl = mpQP.bl + bth 
        f = mpQP.f +mpQP.f_theta*θ
        d = DAQP.Model()
        DAQP.setup(d,mpQP.H,f[:],mpQP.A,bu[:],bl[:],mpQP.senses[:])
        DAQP.settings(d,Dict(:iter_limit=>1e4))
        DAQP.settings(d,Dict(:cycle_tol=>25))
        DAQP.settings(d,Dict(:primal_tol=>1e-6))
        DAQP.settings(d,Dict(:progress_tol=>1e-14))
        #DAQP.settings(d,Dict(:dual_tol=>1e-6))
        DAQP.settings(d,Dict(:pivot_tol=>1e-3))
        udaqp,~,exitflag,info = DAQP.solve(d)
        f_daqp = 0.5*udaqp'*mpQP.H*udaqp+(f'*udaqp)[1]
        slack = minimum([bu-[I(length(udaqp));mpQP.A]*udaqp; [I(length(udaqp));mpQP.A]*udaqp-bl])
        #println(info)

        # Test jump 
        model = create_jump(mpQP.H,f[:],mpQP.A,bu[:],bl[:],findall(mpQP.senses .== 16))
        optimize!(model)
        #display(solution_summary(model))
        u = value.(model[:x])
        #u = zeros(length(udaqp))
        f_jmp = 0.5*u'*mpQP.H*u+(f'*u)[1]
        #u,~,exitflag,info = DAQP.quadprog(mpQP.H,f[:],mpQP.A,bu[:],bl[:],mpQP.senses[:])
        println("flag: $exitflag, iterations: $(info.iterations), nodes: $(info.nodes)")
        #diff = norm(u-udaqp)
        diff = f_jmp-f_daqp 
        println("DIFF: $(diff)")
        display([u udaqp])
        display([f_daqp f_jmp])
        display(slack)
        #@assert(diff<1e0)
        @assert(exitflag==1)
        if(exitflag != 1)
            #@error("DAQP failed")
            #return
            #mpQP.senses[:].=0;
            #u,~,exitflag,info = DAQP.quadprog(mpQP.H,f[:],mpQP.A,bu[:],bl[:],mpQP.senses)
            #println("Completely relaxed: $exitflag") 
        end
    end
    return udaqp[1:mpc.nu],diff
end

function step(mpc::MPC,x;r=nothing,uprev=nothing)
    if(mpc.settings.reference_tracking)
        if(r==nothing) r = zeros(size(mpc.C,1)) end
        if(uprev==nothing) uprev = zeros(mpc.nu) end
        θ = [x;r;uprev]
    else
        θ = x
    end
    u,diff = compute_control(mpc,θ)
    x = mpc.F*x+mpc.G*u
    return x,u,diff
end

function simulate(mpc::MPC,x0,N_steps;r=nothing)
    x = x0
    u = zeros(mpc.nu)
    
    rs = zeros(size(mpc.C,1),N_steps);

    xs = zeros(mpc.nx,N_steps+1); xs[:,1] = x0
    us = zeros(mpc.nu,N_steps)
    diffs = zeros(N_steps)

    # Setup reference 
    if(!isnothing(r))
        rs[:,1:size(r,2)].= r
        rs[:,size(r,2)+1:end] .= r[:,end]
    end

    for i = 1:N_steps
        x,u,diff = step(mpc,x;r=rs[:,i],uprev=u)
        us[:,i] = u 
        xs[:,i+1] = x
        diffs[i] = diff
    end
    return xs,us,rs,diffs
end

function create_jump(H,f,A,bu,bl,bin_ids)
    n = length(f)
    nb = length(bu)-size(A,1)
    if(n > nb)
        bl = [bl[1:nb];-1e30*ones(n-nb);bl[nb+1:end]]
        bu = [bu[1:nb];1e30*ones(n-nb);bu[nb+1:end]]
    end
    model = direct_model(Gurobi.Optimizer())
    set_optimizer_attribute(model, "OutputFlag", 0)
    @variable(model, bl[i] <= x[i = 1:n] <= bu[i])
    for i in bin_ids
        set_binary(x[i])
    end
    @objective(model, Min, 0.5*x'*H*x+f'*x)
    @constraint(model, A*x .<= bu[n+1:end])
    @constraint(model, A*x .>= bl[n+1:end])
    return model
end
