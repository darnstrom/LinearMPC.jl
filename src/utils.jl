"""
    compute_control(mpc,x;r,d,uprev,l)

For a given MPC `mpc` and state `x`, compute the optimal control action.

Optional arguments:
* `r` - reference value. Can be:
  - Vector of length `ny` for constant reference
  - Matrix of size `(ny, Np)` for reference preview (when `mpc.settings.reference_preview = true`)
* `d` - measured disturbance
* `uprev` - previous control action
* `l` - linear cost on control (requires `mpc.settings.linear_cost = true`). Can be:
  - Vector of length `nu` for constant linear cost across horizon
  - Matrix of size `(nu, Np)` for time-varying linear cost (mean computed per move block)
  - Matrix of size `(nu, Nc)` for pre-blocked linear cost
  - Nothing (default) for no linear cost

All arguments default to zero.

# Examples
```julia
# Standard reference tracking
u = compute_control(mpc, x; r=[1.0, 0.0])

# Reference preview (requires mpc.settings.reference_preview = true)
r_trajectory = [1.0 1.5 2.0 2.0 2.0;   # ny × Np matrix
                0.0 0.0 0.5 1.0 1.0]
u = compute_control(mpc, x; r=r_trajectory)

# With linear cost (requires mpc.settings.linear_cost = true)
u = compute_control(mpc, x; l=[1.0])  # Constant cost

# Time-varying linear cost
l_trajectory = [1.0 0.5 0.0 -0.5 -1.0]  # nu × Np matrix
u = compute_control(mpc, x; l=l_trajectory)
```
"""
function compute_control(mpc::MPC,x;r=nothing,d=nothing,uprev=nothing,l=nothing, check=true)
    θ = form_parameter(mpc,x,r,d,uprev,l)
    udaqp,fval,exitflag,info = solve(mpc,θ)
    check && @assert(exitflag>=1)
    # mpc.uprev = udaqp[1:mpc.model.nu]-mpc.K*θ[1:mpc.model.nx]
    mpc.uprev .= udaqp[1:mpc.model.nu]
    mul!(mpc.uprev, mpc.K, θ[1:mpc.model.nx], -1, 1)
    return copy(mpc.uprev)
end

function compute_control(empc::ExplicitMPC,x;r=nothing,d=nothing,uprev=nothing,l=nothing, check=true)
    if isnothing(empc.bst)
        @error "Need to build a binary search tree to evaluate control law"
    else
        return ParametricDAQP.evaluate(empc.bst,form_parameter(empc,x,r,d,uprev,l))
    end
end

function compute_control_trajectory(mpc::MPC,x;r=nothing,d=nothing,uprev=nothing,l=nothing, check=true)
    θ = form_parameter(mpc,x,r,d,uprev,l)
    udaqp,_,exitflag,_ = solve(mpc,θ)
    check && @assert(exitflag>=1)
    # mpc.uprev = udaqp[1:mpc.model.nu]-mpc.K*θ[1:mpc.model.nx]
    mpc.uprev .= udaqp[1:mpc.model.nu]
    mul!(mpc.uprev, mpc.K, θ[1:mpc.model.nx], -1, 1)
    return udaqp
end

"""
    format_reference(mpc, r)

Format reference input for MPC controller. Handles both single reference 
and reference preview scenarios.
"""
function format_reference(mpc::Union{MPC,ExplicitMPC}, r)
    !mpc.settings.reference_tracking && return zeros(0)
    if isnothing(r) 
        r = zeros(mpc.model.ny)
    end
    isempty(r) && return r
    
    if mpc.settings.reference_preview
        # Reference preview mode
        ny = mpc.model.ny
        
        if r isa AbstractVector
            # Single reference - broadcast across prediction horizon
            if length(r) == ny
                return condense_reference(mpc,vec(repeat(r, mpc.Np)))
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
                return condense_reference(mpc,vec(r[:, 1:mpc.Np]))
            else
                # Pad with last column if trajectory is shorter than prediction horizon
                r_extended = zeros(ny, mpc.Np)
                r_extended[:, 1:size(r, 2)] = r
                r_extended[:, (size(r, 2)+1):end] .= repeat(r[:, end], 1, mpc.Np - size(r, 2))
                return condense_reference(mpc,vec(r_extended))
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
"""
    condense_reference(mpc, r)

Condense reference trajectory to single setpoint.
"""
function condense_reference(mpc::Union{MPC,ExplicitMPC}, r)
    if mpc.settings.reference_condensation
        isempty(mpc.traj2setpoint) && setup!(mpc)
        return mpc.traj2setpoint*r
    else
        return r
    end
end

"""
    format_linear_cost(mpc, l)

Format linear cost input for MPC controller.

# Arguments
- `mpc`: MPC or ExplicitMPC controller
- `l`: Linear cost input. Can be:
  - `nothing`: Returns zero vector
  - Vector of length `nu`: Broadcast across control horizon
  - Vector of length `nl`: Used directly (already blocked)
  - Matrix of size `(nu, Np)`: Mean computed over each move block
  - Matrix of size `(nu, Nc)`: Used directly
"""
function format_linear_cost(mpc::Union{MPC,ExplicitMPC}, l)
    # Use stored mpc.nl (set during setup!) to ensure consistency with the QP
    mpc.nl == 0 && return zeros(0)
    isnothing(l) && return zeros(mpc.nl)

    nu = mpc.model.nu
    Nc_blocked = mpc.nl ÷ nu  # Number of blocked control moves
    move_blocks = mpc.move_blocks
    Np = mpc.Np

    # Handle single vector (constant linear cost)
    if l isa AbstractVector && length(l) == nu
        return repeat(l, Nc_blocked)
    end

    # Handle already-blocked input
    if l isa AbstractVector && length(l) == mpc.nl
        return l
    end

    # Matrix input
    if l isa AbstractMatrix
        size(l, 1) != nu && error("Linear cost matrix must have $nu rows")
        ncols = size(l, 2)

        if ncols == Nc_blocked
            # Already blocked size - use directly
            return vec(l)
        elseif ncols >= Np || !isempty(move_blocks)
            # Full Np horizon (or longer) - apply move blocking
            l_mat = ncols >= Np ? l[:, 1:Np] : begin
                # Pad with last column
                tmp = zeros(nu, Np)
                tmp[:, 1:ncols] = l
                tmp[:, ncols+1:end] .= l[:, end]
                tmp
            end

            if isempty(move_blocks)
                # No blocking - use first Nc_blocked columns
                return vec(l_mat[:, 1:Nc_blocked])
            end

            # Apply move blocking via direct mean computation
            l_blocked = zeros(mpc.nl)
            offset = 0
            for (k, mb) in enumerate(move_blocks)
                block_inds = (offset + 1):min(offset + mb, Np)
                l_blocked[(k-1)*nu+1:k*nu] = vec(mean(l_mat[:, block_inds], dims=2))
                offset += mb
            end
            return l_blocked
        else
            # Smaller than Np but not Nc_blocked - pad to Nc_blocked
            l_ext = zeros(nu, Nc_blocked)
            l_ext[:, 1:ncols] = l
            l_ext[:, ncols+1:end] .= l[:, end]
            return vec(l_ext)
        end
    else
        error("Linear cost must be a vector or matrix")
    end
end

"""
    xdaqp,fval,exitflag,info = solve(mpc,θ)

Solve corresponding QP given the parameter θ
"""
function solve(mpc::MPC,θ)
    mpc.mpqp_issetup || setup!(mpc) # ensure mpQP is setup
    mpc.mpqp_issetup || throw("Could not setup optimization problem")
    _inner_solve(mpc, mpc.mpQP, mpc.opt_model, θ)
end

# The function below is a "function barrier" to work around mpQP and opt_model being
# type unstable
function _inner_solve(mpc, mpQP, opt_model, θ)
    mul!(mpc._bth, mpQP.W, θ)
    mpc._bu .= mpQP.bu .+ mpc._bth
    mpc._bl .= mpQP.bl .+ mpc._bth
    mul!(mpc._f, mpQP.f_theta, θ)
    mpc._f .+= mpQP.f
    if mpQP.has_binaries # Make sure workspace is clean
        ccall(("node_cleanup_workspace", DAQP.libdaqp),Cvoid,(Cint,Ptr{DAQP.Workspace}),0, opt_model.work)
    end
    DAQP.update(opt_model,nothing,mpc._f,nothing,mpc._bu,mpc._bl,nothing)
    return DAQP.solve(opt_model)
end

function range2region(range)
    lb = [range.xmin;range.rmin;range.dmin;range.umin;range.lmin]
    ub = [range.xmax;range.rmax;range.dmax;range.umax;range.lmax]
    return (A = zeros(length(ub), 0), b=zeros(0), lb=lb,ub=ub)
end

function zoh(A,B,Ts)
  nx,nu= size(B);
  M = exp([A*Ts  B*Ts; zeros(nu, nx + nu)])
  return M[1:nx, 1:nx], M[1:nx, nx+1:nx+nu]
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

# 1. Transform bl + W θ ≤ A U ≤ bu + W θ → A U ≤ b + W
# 2. Make the soft constraints explicit
function make_singlesided(mpQP;single_soft=false, soft_weight=1e6)
    ncstr = length(mpQP.bu);
    n_bounds = ncstr-size(mpQP.A,1);
    bounds_table=[collect(ncstr+1:2*ncstr);collect(1:ncstr)]
    A0 = [I(n_bounds) zeros(n_bounds,size(mpQP.A,2)-n_bounds);mpQP.A]
    A = [A0;-A0]

    senses = repeat(mpQP.senses,2)
    prio = repeat(mpQP.prio,2)

    H,f,f_theta = mpQP.H, mpQP.f, mpQP.f_theta

    # Make soft constraints explicit
    soft_mask = (mpQP.senses .& DAQP.SOFT .== DAQP.SOFT)
    if(any(soft_mask))
        soft_ids = findall(soft_mask)

        R = cholesky((mpQP.H+mpQP.H')/2)
        Ms= A0[soft_mask,:]/R.U
        norm_factors = [norm(view(Ms,i,:),2) for i in 1:size(Ms,1)]

        if(single_soft)
            nsoft = 1
            A = [A zeros(size(A,1),1)]
            A[soft_ids,end] .= -norm_factors
            A[soft_ids .+ ncstr,end] .= -norm_factors
        else
            nsoft, n = length(soft_ids), size(A,2)
            A = [A zeros(2*ncstr,nsoft)]
            A[soft_ids,n+1:end] = -diagm(norm_factors)
            A[soft_ids .+ ncstr,n+1:end] = -diagm(norm_factors)

        end
        H = cat(H,soft_weight*I(nsoft),dims=(1,2))
        f = [f;zeros(nsoft)] 
        f_theta = [f_theta;zeros(nsoft,size(f_theta,2))]
    end

    b = [mpQP.bu;-mpQP.bl]
    W = [mpQP.W;-mpQP.W]

    # Prune possible Inf bounds
    rm_ids = findall(b[:] .>= 1e20)
    if(!isempty(rm_ids))
        bounds_table[bounds_table[rm_ids]] = bounds_table[rm_ids] # Make other bound point to itself
        # Correct bounds table 
        rm_offset, keep_ids = 1, Int[]
        for i in 1:2*ncstr
            if(i==rm_ids[rm_offset])
                rm_offset+=1
            else
                bounds_table[i] -= (rm_offset-1)
                push!(keep_ids,i)
            end
        end
        A,b,W = A[keep_ids,:],b[keep_ids],W[keep_ids,:]
        senses,prio,bounds_table = senses[keep_ids],prio[keep_ids],bounds_table[keep_ids]
    end

    return (H=H,f=f, H_theta = mpQP.H_theta, f_theta=f_theta,
            A=Matrix{Float64}(A), b=b, W=W, senses=senses,
            bounds_table=bounds_table, prio =prio, has_binaries=mpQP.has_binaries)
end

"""
    evaluate_cost(mpc,xs,us,rs;Q,Rr,S)
Compute the cost 0.5 ∑ x'*Q x + u' R u + Δu' Rr Δu + x' S u
"""
function evaluate_cost(mpc::MPC,xs,us,rs=zeros(0,0);
        Q=mpc.weights.Q, R = mpc.weights.R, Rr = mpc.weights.Rr, S = mpc.weights.S)
    nu,N = size(us)
    rs = isempty(rs) ? zeros(mpc.model.ny,N) : rs
    Δus = diff([zeros(nu) us],dims=2)
    cost = 0.0
    for i = 1:N
        err = mpc.model.C*xs[:,i]-rs[:,i]
        cost += dot(err,Q,err)
        cost += dot(us[:,i],R,us[:,i])
        cost += dot(Δus[:,i],Rr,Δus[:,i])
        cost += dot(xs[:,i],S,us[:,i])
    end
    return 0.5*cost
end

"""
    constraint_violation(c,xs,us)
evaluates the possible violation of constraint c at state x and control u
"""
function constraint_violation(c::Constraint,x::Vector{Float64},u::Vector{Float64})
    Axx_Aux = c.Ax*x+c.Au*u
    return maximum([c.lb-Axx_Aux;Axx_Aux-c.ub;0])
end

function constraint_violation(c::Constraint,xs::Matrix{Float64},us::Matrix{Float64})
    @assert(size(xs,2) == size(us,2))
    return [constraint_violation(c,xs[:,i],us[:,i]) for i = 1:size(xs,2)]
end
