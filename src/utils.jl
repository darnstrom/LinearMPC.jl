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
        u,~,exitflag,info = DAQP.quadprog(mpQP.H,f[:],mpQP.A,bu[:],bl[:],mpQP.senses)
        println(exitflag)
    end
    return u[1:mpc.nu]
end

function step(mpc::MPC,x;r=nothing,uprev=nothing)
    if(mpc.settings.reference_tracking)
        if(r==nothing) r = zeros(size(mpc.C,1)) end
        if(uprev==nothing) uprev = zeros(mpc.nu) end
        θ = [x;r;uprev]
    else
        θ = x
    end
    u = compute_control(mpc,θ)
    x = mpc.F*x+mpc.G*u
    return x,u
end

function simulate(mpc::MPC,x0,N_steps;r=nothing)
    x = x0
    u = zeros(mpc.nu)
    
    rs = zeros(size(mpc.C,1),N_steps);

    xs = zeros(mpc.nx,N_steps+1); xs[:,1] = x0
    us = zeros(mpc.nu,N_steps)

    # Setup reference 
    if(!isnothing(r))
        rs[:,1:size(r,2)].= r
        rs[:,size(r,2)+1:end] .= r[:,end]
    end

    for i = 1:N_steps
        x,u = step(mpc,x;r=rs[:,i],uprev=u)
        us[:,i] = u 
        xs[:,i+1] = x
    end
    return xs,us,rs
end
