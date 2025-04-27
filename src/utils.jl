"""
    compute_control(mpc,x;r,uprev)

For a given MPC `mpc` and state `x`, compute the optimal control action. 
Optional arguments: 
* `r` - reference value
* `uprev` - previous control action
all of them defaults to zero.
"""
function compute_control(mpc::Union{MPC,ExplicitMPC},x;r=nothing,d=nothing,uprev=nothing)
    # Setup parameter vector θ
    nx,nr,nd,nuprev = get_parameter_dims(mpc)
    r = isnothing(r) ? zeros(nr) : r
    d = isnothing(d) ? zeros(nd) : d 
    uprev = isnothing(uprev) ? zeros(nuprev) : uprev
    return solve(mpc,[x;r;d;uprev])
end

function solve(mpc::MPC,θ)
    isnothing(mpc.mpQP) && setup!(mpc) # ensure mpQP is setup 
    bth = mpc.mpQP.W*θ
    bu = mpc.settings.QP_double_sided ? mpc.mpQP.bu + bth : mpc.mpQP.b + bth
    bl = mpc.settings.QP_double_sided ? mpc.mpQP.bl + bth : -1e30*ones(length(bu))
    f = mpc.mpQP.f +mpc.mpQP.f_theta*θ
    DAQP.update(mpc.opt_model,nothing,f,nothing,bu,bl,nothing)
    udaqp,fval,exitflag,info = DAQP.solve(mpc.opt_model)
    @assert(exitflag>=1)
    return udaqp[1:mpc.model.nu]-mpc.K*θ[1:mpc.model.nx]
end

function solve(empc::ExplicitMPC,θ)
    if isnothing(empc.bst) 
        @error "Need to build a binary search tree to evaluate control law"
    else
        return ParametricDAQP.evaluate(empc.bst,θ) 
    end
end


function simulate(dynamics,mpc::Union{MPC,ExplicitMPC},x0,N_steps;r=nothing, callback=(x,u,k)->nothing)
    x,u = x0,zeros(mpc.model.nu)
    
    rs = zeros(mpc.model.ny,N_steps);
    xs = zeros(mpc.model.nx,N_steps+1); xs[:,1] = x0
    us = zeros(mpc.model.nu,N_steps)

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

function range2region(range)
    lb = [range.xmin;range.rmin;range.dmin;range.umin]
    ub = [range.xmax;range.rmax;range.dmax;range.umax]
    return (A = zeros(length(ub), 0), b=zeros(0), lb=lb,ub=ub)
end

function zoh(A,B,Ts; Bd = zeros(0,0))
  dims = size(B);
  nx,nu = length(dims)==1 ? (dims[1],1) : dims
  if isempty(Bd)
      M = exp([A*Ts  B*Ts; zeros(nu, nx + nu)])
      return M[1:nx, 1:nx], M[1:nx, nx+1:nx+nu],zeros(0,0)
  else
      nd = size(Bd,2);
      M = exp([A*Ts  [B Bd]*Ts; zeros(nu+nd, nx + nu + nd)])
      return M[1:nx, 1:nx], M[1:nx, nx+1:nx+nu], M[1:nx, nx+nu+1:nx+nu+nd]
  end
end

matrixify(x::Number, n::Int) = diagm(fill(float(x),n))
matrixify(v::AbstractVector, n::Int) = diagm(float(v))
matrixify(M::AbstractMatrix, n::Int) = float(M)

# get paramer id of a label
function label2id(mpc, label::Symbol)

    id = findfirst(x->x==label,mpc.model.labels.x)
    isnothing(id)  || return id,string(label);

    if(mpc.nr > 0 && string(label)[end] == 'r')
        l = Symbol(string(label)[1:end-1])
        id = findfirst(x->x==l,mpc.model.labels.y)
        isnothing(id)  || return mpc.model.nx + id,string(l)*"^r";
    end

    id = findfirst(x->x==label,mpc.model.labels.d)
    isnothing(id)  || return mpc.model.nx+mpc.nr+id,string(label);

    if(mpc.nuprev > 0 && string(label)[end] == 'p')
        l = Symbol(string(label)[1:end-1])
        id = findfirst(x->x==l,mpc.model.labels.u)
        isnothing(id)  || return mpc.model.nx+mpc.nr+mpc.model.nd+id,string(l)*"^-";
    end

    return nothing,string(label)
end

function make_subscript(label::String)
    nid = findfirst(isdigit,collect(label))
    return isnothing(nid) ? label : label[1:nid-1]*"_"*label[nid:end]
end
