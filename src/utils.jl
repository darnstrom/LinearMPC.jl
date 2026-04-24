"""
    compute_control(mpc,x;r,d,uprev,p)

For a given MPC `mpc` and state `x`, compute the optimal control action.

Optional arguments:
* `r` - reference value. Can be:
  - Vector of length `ny` for constant reference
  - Matrix of size `(ny, Np)` for reference preview (when `mpc.settings.reference_preview = true`)
* `d` - measured disturbance. Can be:
  - Vector of length `nd` for constant disturbance
  - Matrix of size `(nd, Np)` for disturbance preview (when `mpc.settings.disturbance_preview = true`)
* `uprev` - previous control action
* `p` - generalized parameter trajectory. Can be:
  - Vector of length `np` for a constant parameter across the horizon
  - Matrix of size `(np, Np)` for time-varying generalized parameters
  - Nothing (default) for no generalized-parameter contribution

All arguments default to zero.

# Examples
```julia
# Standard reference tracking
u = compute_control(mpc, x; r=[1.0, 0.0])

# Reference preview (requires mpc.settings.reference_preview = true)
r_trajectory = [1.0 1.5 2.0 2.0 2.0;   # ny × Np matrix
                0.0 0.0 0.5 1.0 1.0]
u = compute_control(mpc, x; r=r_trajectory)

# Disturbance preview (requires mpc.settings.disturbance_preview = true)
d_trajectory = [0.0 0.2 0.4 0.4 0.4]  # nd × Np matrix
u = compute_control(mpc, x; d=d_trajectory)

# Generalized parameter preview
u = compute_control(mpc, x; p=[0.5, -0.25])
```
"""
function compute_control(mpc::MPC,x;r=nothing,d=nothing,uprev=nothing,p=nothing, check=true)
    θ = form_parameter(mpc,x,r,d,uprev,p)
    udaqp,fval,exitflag,info = solve(mpc,θ)
    check && @assert(exitflag>=1)
    # mpc.uprev = udaqp[1:mpc.model.nu]-mpc.K*θ[1:mpc.model.nx]
    mpc.uprev .= udaqp[1:mpc.model.nu]
    mul!(mpc.uprev, mpc.K, θ[1:mpc.model.nx], -1, 1)
    return copy(mpc.uprev)
end

function compute_control(empc::ExplicitMPC,x;r=nothing,d=nothing,uprev=nothing,p=nothing, check=true)
    if isnothing(empc.bst)
        @error "Need to build a binary search tree to evaluate control law"
    else
        empc.uprev .= ParametricDAQP.evaluate(empc.bst,form_parameter(empc,x,r,d,uprev,p))
        return copy(empc.uprev) 
    end
end

function compute_control_trajectory(mpc::MPC,x;r=nothing,d=nothing,uprev=nothing,p=nothing, check=true)
    θ = form_parameter(mpc,x,r,d,uprev,p)
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
    format_disturbance(mpc, d)

Format disturbance input for MPC controller. Handles both single disturbance
and disturbance preview scenarios.
"""
function format_disturbance(mpc::Union{MPC,ExplicitMPC}, d)
    d = get_control_disturbance(mpc, d)
    nd_base = mpc.model.nd
    nd_base == 0 && return zeros(0)

    if isnothing(d)
        d = zeros(nd_base)
    end
    isempty(d) && return d

    if mpc.settings.disturbance_preview
        if d isa AbstractVector
            if length(d) == nd_base
                return vec(repeat(d, 1, mpc.Np))
            else
                error("Disturbance vector length ($(length(d))) must match number of disturbances ($(nd_base))")
            end
        elseif d isa AbstractMatrix
            if size(d, 1) != nd_base
                error("Disturbance matrix must have $(nd_base) rows (number of disturbances)")
            end

            if size(d, 2) >= mpc.Np
                return vec(d[:, 1:mpc.Np])
            else
                d_extended = zeros(nd_base, mpc.Np)
                d_extended[:, 1:size(d, 2)] = d
                d_extended[:, (size(d, 2)+1):end] .= repeat(d[:, end], 1, mpc.Np - size(d, 2))
                return vec(d_extended)
            end
        else
            error("Disturbance must be a vector or matrix for disturbance preview")
        end
    else
        if d isa AbstractVector
            if length(d) == nd_base
                return d
            else
                error("Disturbance vector length ($(length(d))) must match number of disturbances ($(nd_base))")
            end
        elseif d isa AbstractMatrix
            if size(d, 1) == nd_base
                return d[:, 1]
            else
                error("Disturbance matrix must have $(nd_base) rows (number of disturbances)")
            end
        else
            error("Disturbance must be a vector or matrix")
        end
    end
end

function get_affine_parameter_base_dim(mpc::MPC)
    mpc.mpqp_issetup && return mpc.np == 0 ? 0 : mpc.np ÷ mpc.Np
    dims = Int[size(mpc.weights.Ex, 2), size(mpc.weights.E, 2)]
    append!(dims, (size(c.Ap, 2) for c in mpc.constraints))
    append!(dims, (size(first(c).Ex, 2) for c in mpc.objectives))
    append!(dims, (size(first(c).E, 2) for c in mpc.objectives))
    return maximum(dims)
end

get_affine_parameter_base_dim(mpc::ExplicitMPC) = mpc.np == 0 ? 0 : mpc.np ÷ mpc.Np

"""
    format_affine_parameters(mpc, p)

Format generalized parameter input for MPC controller.
"""
function format_affine_parameters(mpc::Union{MPC,ExplicitMPC}, p)
    np_total = mpc isa MPC && !mpc.mpqp_issetup ? get_affine_parameter_base_dim(mpc) * mpc.Np : mpc.np
    np_total == 0 && return zeros(0)
    isnothing(p) && return zeros(np_total)

    np_base = get_affine_parameter_base_dim(mpc)
    Np = mpc.Np

    if p isa AbstractVector && length(p) == np_base
        return vec(repeat(p, 1, Np))
    end

    if p isa AbstractVector && length(p) == np_total
        return float(p)
    end

    if p isa AbstractMatrix
        size(p, 1) == np_base || error("Generalized parameter matrix must have $np_base rows")
        if size(p, 2) >= Np
            return vec(p[:, 1:Np])
        end

        p_extended = zeros(np_base, Np)
        p_extended[:, 1:size(p, 2)] = p
        p_extended[:, size(p, 2)+1:end] .= repeat(p[:, end], 1, Np - size(p, 2))
        return vec(p_extended)
    end

    error("Generalized parameters must be a vector or matrix")
end

"""
    xdaqp,fval,exitflag,info = solve(mpc,θ)

Solve corresponding QP given the parameter θ
"""
function solve(mpc::MPC,θ)
    mpc.mpqp_issetup || setup!(mpc) # ensure mpQP is setup
    mpc.mpqp_issetup || throw("Could not setup optimization problem")

    mul!(mpc.mpQP._bth, mpc.mpQP.W, θ)
    mpc.mpQP._bu .= mpc.mpQP.bu .+ mpc.mpQP._bth
    mpc.mpQP._bl .= mpc.mpQP.bl .+ mpc.mpQP._bth
    mul!(mpc.mpQP._f, mpc.mpQP.f_theta, θ)
    mpc.mpQP._f .+= mpc.mpQP.f
    if mpc.mpQP.has_binaries # Make sure workspace is clean
        ccall(("daqp_node_cleanup_workspace", DAQP.libdaqp),
              Cvoid,(Cint,Ptr{DAQP.Workspace}),0, mpc.opt_model.work)
    end
    DAQP.update(mpc.opt_model,nothing,mpc.mpQP._f,nothing,mpc.mpQP._bu,mpc.mpQP._bl,nothing)
    return DAQP.solve(mpc.opt_model)
end

function range2region(range)
    lb = [range.xmin;range.rmin;range.dmin;range.umin;range.pmin]
    ub = [range.xmax;range.rmax;range.dmax;range.umax;range.pmax]
    return (A = zeros(length(ub), 0), b=zeros(0), lb=lb,ub=ub)
end

function zoh(A,B,Ts)
  nx,nu= size(B);
  M = exp([A*Ts  B*Ts; zeros(nu, nx + nu)])
  return M[1:nx, 1:nx], M[1:nx, nx+1:nx+nu]
end

matrixify(x::Number, n::Int) = diagm(fill(float(x),n))
matrixify(v::AbstractVector, n::Int=length(v)) = diagm(float(v))
matrixify(M::AbstractMatrix, n::Int=size(M,1)) = float(M)

function prettify_parameter_label(label::Symbol)
    str = string(label)
    if endswith(str, "p")
        return string(Symbol(str[1:end-1]))*"^-"
    elseif endswith(str, "r")
        return string(Symbol(str[1:end-1]))*"^r"
    elseif occursin(r"r_\d+$", str)
        base, step = split(str, "_"; limit=2)
        return string(Symbol(base[1:end-1]))*"^r_"*step
    else
        return str
    end
end

# get paramer id of a label
function label2id(mpc, label::Symbol)
    names = get_parameter_names(mpc)
    id = findfirst(x -> x == label, names)
    return isnothing(id) ? (nothing, string(label)) : (id, prettify_parameter_label(label))
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
            if(rm_offset <= length(rm_ids) && i==rm_ids[rm_offset])
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
