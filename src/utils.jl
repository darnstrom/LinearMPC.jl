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
    r = isnothing(r) || !mpc.settings.reference_tracking ? zeros(nr) : r
    d = isnothing(d) ? zeros(nd) : d 
    uprev = isnothing(uprev) || mpc.nuprev == 0 ? zeros(nuprev) : uprev
    return solve(mpc,[x;r;d;uprev])
end

function solve(mpc::MPC,θ)
    mpc.mpqp_issetup || setup!(mpc) # ensure mpQP is setup
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
function make_singlesided(mpQP;explicit_soft=true)
    ncstr = length(mpQP.bu);
    n_bounds = ncstr-size(mpQP.A,1);
    bounds_table=[collect(ncstr+1:2*ncstr);collect(1:ncstr)]
    A = [I(n_bounds) zeros(n_bounds,size(mpQP.A,2)-n_bounds);mpQP.A]
    A = [A;-A]
    if(explicit_soft && any(senses.==DAQP.SOFT))# Correct sign for slack
        A[:,end].= -abs.(A[:,end])
    end
    b = [mpQP.bu;-mpQP.bl]
    W = [mpQP.W;-mpQP.W]
    senses = [mpQP.senses;mpQP.senses]
    return (H=mpQP.H,f=mpQP.f, H_theta = mpQP.H_theta, f_theta=mpQP.f_theta,
            A=Matrix{Float64}(A), b=b, W=W, senses=senses,
            bounds_table=bounds_table)
end
