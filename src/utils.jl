"""
    compute_control(mpc,x;r,uprev)

For a given MPC `mpc` and state `x`, compute the optimal control action. 

Optional arguments: 
* `r` - reference value. Can be:
  - Vector of length `ny` for constant reference
  - Matrix of size `(ny, Np)` for reference preview (when `mpc.settings.reference_preview = true`)
* `uprev` - previous control action

All arguments default to zero.

# Examples
```julia
# Standard reference tracking
u = compute_control(mpc, x; r=[1.0, 0.0])

# Reference preview (requires mpc.settings.reference_preview = true)
r_trajectory = [1.0 1.5 2.0 2.0 2.0;   # ny × Np matrix
                0.0 0.0 0.5 1.0 1.0]
u = compute_control(mpc, x; r=r_trajectory)
```
"""
function compute_control(mpc::Union{MPC,ExplicitMPC},x;r=nothing,d=nothing,uprev=nothing)
    # Setup parameter vector θ
    nx,nr,nd,nuprev = get_parameter_dims(mpc)
    
    r = format_reference(mpc, r)
    d = isnothing(d) ? zeros(nd) : d 
    uprev = isnothing(uprev) ? mpc.uprev[1:nuprev] : uprev[1:nuprev]

    mpc.uprev = solve(mpc,[x;r;d;uprev])
    return mpc.uprev
end

"""
    format_reference(mpc, r)

Format reference input for MPC controller. Handles both single reference 
and reference preview scenarios.
"""
function format_reference(mpc::Union{MPC,ExplicitMPC}, r)
    if isnothing(r)
        r = mpc.settings.reference_tracking ?  zeros(mpc.model.ny) : zeros(0)
    end
    isempty(r) && return r
    
    if mpc.settings.reference_preview
        # Get ids
        if(isempty(mpc.ref_blocks))
            ids = 1:mpc.Np
        else
            id,ids = 1,Int[]
            for mb in mpc.ref_blocks
                push!(ids,id)
                id += mb
            end
        end


        # Reference preview mode
        ny = mpc.model.ny
        
        if r isa AbstractVector
            # Single reference - broadcast across prediction horizon
            if length(r) == ny
                return repeat(r, length(ids))
            else
                error("Reference vector length ($(length(r))) must match number of outputs ($(ny))")
            end
        elseif r isa AbstractMatrix
            # Reference trajectory matrix
            if size(r, 1) != ny
                error("Reference matrix must have $(ny) rows (number of outputs)")
            end
            
            # Flatten and pad/truncate to prediction horizon
            if size(r, 2) >= mpc.Np
                return vec(r[:, ids])
            else
                # Pad with last column if trajectory is shorter than prediction horizon
                r_extended = zeros(ny, mpc.Np)
                r_extended[:, 1:size(r, 2)] = r
                r_extended[:, (size(r, 2)+1):end] .= repeat(r[:, end], 1, mpc.Np - size(r, 2))
                return vec(r_extended[:,ids])
            end
        else
            error("Reference must be a vector or matrix for reference preview")
        end
    else
        # Standard mode - single reference
        if r isa AbstractVector
            if length(r) == mpc.model.ny
                return r
            else
                error("Reference vector length ($(length(r))) must match number of outputs ($(mpc.model.ny))")
            end
        elseif r isa AbstractMatrix
            # Use first column if matrix provided in non-preview mode
            if size(r, 1) == mpc.model.ny
                return r[:, 1]
            else
                error("Reference matrix must have $(mpc.model.ny) rows (number of outputs)")
            end
        else
            error("Reference must be a vector or matrix")
        end
    end
end

function solve(mpc::MPC,θ)
    mpc.mpqp_issetup || setup!(mpc) # ensure mpQP is setup
    mpc.mpqp_issetup || throw("Could not setup optimization problem")
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
        l = Symbol(split(string(label),'r')[1])
        id = findfirst(x->x==l,mpc.model.labels.y)
        isnothing(id)  || return mpc.model.nx + id,string(l)*"^r";
    end

    id = findfirst(x->x==label,mpc.model.labels.d)
    isnothing(id)  || return mpc.model.nx+mpc.nr+id,string(label);

    if(mpc.nuprev > 0 && string(label)[end] == 'p')
        l = Symbol(split(string(label),'p')[1])
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
