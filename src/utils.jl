function zoh(A,B,Ts)
  dims = size(B);
  nx,nu = length(dims)==1 ? (dims[1],1) : dims
  M = exp([A*Ts  B*Ts; zeros(nu, nx + nu)])
  return M[1:nx, 1:nx], M[1:nx, nx+1:nx+nu]
end

function solve(mpc::MPC,θ)
    bth = mpc.mpQP.W*θ
    bu = mpc.settings.QP_double_sided ? mpc.mpQP.bu + bth : mpc.mpQP.b + bth
    bl = mpc.settings.QP_double_sided ? mpc.mpQP.bl + bth : nothing 
    f = mpc.mpQP.f +mpc.mpQP.f_theta*θ
    DAQP.update(mpc.opt_model,nothing,f,nothing,bu,bl,nothing)
    udaqp,fval,exitflag,info = DAQP.solve(mpc.opt_model)
    @assert(exitflag>=1)
    return udaqp[1:mpc.nu]
end

function solve(empc::ExplicitMPC,θ)
    if isnothing(empc.bst) 
        @error "Need to build a binary search tree to evaluate control law"
    else
        return ParametricDAQP.evaluate(empc.bst,θ) 
    end
end

function compute_control(mpc::Union{MPC,ExplicitMPC},x;r=zeros(mpc.ny),uprev=zeros(mpc.nu))
    θ = mpc.settings.reference_tracking ? [x;r;uprev] : x
    return solve(mpc,θ)
end

function simulate(dynamics,mpc::Union{MPC,ExplicitMPC},x0,N_steps;r=nothing, callback=(x,u,k)->nothing)
    x,u = x0,zeros(mpc.nu)
    
    rs = zeros(mpc.ny,N_steps);
    xs = zeros(mpc.nx,N_steps+1); xs[:,1] = x0
    us = zeros(mpc.nu,N_steps)

    # Setup reference 
    if(!isnothing(r))
        rs[:,1:size(r,2)].= r
        rs[:,size(r,2)+1:end] .= r[:,end] # hold last value of reference
    end

    # Start the simulation
    for k = 1:N_steps
        u = compute_control(mpc,x;r=rs[:,k],uprev=u)
        x = dynamics(x,u)
        callback(x,u,k)
        us[:,k], xs[:,k+1] = u, x
    end
    return xs,us,rs
end
