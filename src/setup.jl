"""
    setup!(mpc)

Sets up the `mpc` given its current parameters and settings  
Internally, this means generating an mpQP, and setting up a DAQP workspace.
"""
function setup!(mpc::MPC)
    mpc.mpQP = mpc2mpqp(mpc)
    bu,bl = mpc.mpQP.bu[:],mpc.mpQP.bl[:]
    setup_flag,_ = DAQP.setup(mpc.opt_model, mpc.mpQP.H,mpc.mpQP.f[:],mpc.mpQP.A,bu,bl,mpc.mpQP.senses;break_points=mpc.mpQP.break_points)
    if(setup_flag < 0)
        @warn " Cannot setup optimization problem " setup_flag
    else
        # Set up soft weight
        DAQP.settings(mpc.opt_model,Dict(:rho_soft=>1/mpc.settings.soft_weight))
        mpc.mpqp_issetup = true
    end
end

"""
    set_input_bounds!(mpc;umin,umax)

Sets the input bounds umin ≤ u ≤ umax 
"""
function set_input_bounds!(mpc::MPC; umin=zeros(0), umax=zeros(0))
    nmin,nmax = length(umin),length(umax)
    nb = max(nmin,nmax)
    nb == 0 && return
    nb != mpc.model.nu  && @error("# of controls are $(mpc.model.nu), got bounds of dimension $nb")

    mpc.umin = [umin;-1e30*ones(nb-nmin)]
    mpc.umax = [umax;+1e30*ones(nb-nmax)]
    mpc.mpqp_issetup = false
end

"""
    add_constraint!(mpc::MPC;
        Ax, Au, Ar, Aw, Ad, Aup,
        ub, lb, ks, soft, binary, prio)
    add_constraint!(mpc;Ax,Au,ub,lb,
                    ks, soft, binary,prio)

Adds the constraints lb ≤ Ax xₖ + Au uₖ ≤ ub for the time steps k ∈ ks
(additional terms Ar rₖ, Aw wₖ, Ad dₖ, Aup u⁻ₖ are possible)

* `soft` marks if the constraint should be softened (default false)
* `binary` marks if either the upper or lower bounds should be enforced with equality (default false)
* `prio` marks the relative priority of the constraint (default 0)
"""
function add_constraint!(mpc::MPC;
        Ax = nothing, Au= nothing, Ar = zeros(0,0), Aw = zeros(0,0), Ad = zeros(0,0), Aup = zeros(0,0),
        ub = zeros(0), lb = zeros(0),
        ks = 2:mpc.Np, soft=false, binary=false, prio = 0)
    if isnothing(Ax) && isnothing(Au)
        return
    end

    # Get length of constraint
    nlb,nub = length(lb),length(ub) 
    m = max(nlb,nub)
    m == 0 && return

    ub = nub == m ? ub : [ub;1e30*ones(m-nub)]
    lb = nlb == m ? lb : [lb;-1e30*ones(m-nlb)]

    Ax = isnothing(Ax) ? zeros(m,mpc.model.nx) : Ax
    Au = isnothing(Au) ? zeros(m,mpc.model.nu) : Au

    push!(mpc.constraints,Constraint(Au,Ax,Ar,Aw,Ad,Aup,ub,lb,ks,soft,binary,prio))
    mpc.mpqp_issetup = false
end

"""
    set_output_bounds!(mpc;ymin,ymax,
                    ks, soft, binary,prio)

Adds the constraints lb ≤ C x  ≤ ub for the time steps k ∈ ks 

* `soft` marks if the constraint should be softened (default false)
* `binary` marks if either the upper or lower bounds should be enforced with equality (default false)
* `prio` marks the relative priority of the constraint (default 0)
"""
function set_output_bounds!(mpc::MPC; ymin=zeros(0), ymax=zeros(0), ks = 2:mpc.Np, soft = true, binary=false, prio = 0)
    add_constraint!(mpc, Ax = mpc.model.C, Ad = mpc.model.Dd, lb = ymin, ub = ymax;ks,soft,binary,prio)
end

"""
    set_bounds!(mpc;umin,umax,ymin,umax)

Sets the bounds umin ≤ u ≤ umax and ymin ≤ y ≤ umax
"""
function set_bounds!(mpc::MPC; umin=zeros(0), umax=zeros(0), ymin = zeros(0), ymax = zeros(0))
    (!isempty(umin) ||  !isempty(umax)) && set_input_bounds!(mpc;umin,umax)
    (!isempty(ymin) ||  !isempty(ymax)) && set_output_bounds!(mpc;ymin,ymax)
end

"""
    set_objective!(mpc;Q,R,Rr,S,Qf)

Set the weights in the objective function `xN' C' Qf C xN^T + ∑ (C xₖ - rₖ)' Q (C xₖ - rₖ)  + uₖ' R uₖ + Δuₖ' Rr Δuₖ + xₖ' S uₖ

A vector is interpreted as a diagonal matrix.
"""
function set_objective!(mpc::MPC;Q = zeros(0,0), R=zeros(0,0), Rr=zeros(0,0), S= zeros(0,0), Qf=zeros(0,0), Qfx=zeros(0,0))
    isempty(Q)  || (mpc.weights.Q .= matrixify(Q,mpc.model.ny))
    isempty(R)  || (mpc.weights.R .= matrixify(R,mpc.model.nu))
    isempty(Rr) || (mpc.weights.Rr .= matrixify(Rr,mpc.model.nu))
    isempty(S)  || (mpc.weights.S .= float(S))
    isempty(Qf) || (mpc.weights.Qf .= matrixify(Qf,mpc.model.ny))
    isempty(Qfx) || (mpc.weights.Qfx .= matrixify(Qfx,mpc.model.nx))
    mpc.mpqp_issetup = false
end
set_weights! = set_objective! # backwards compatibility

# Terminal ingredients
using MatrixEquations 

"""
    set_terminal_cost!(mpc)

Sets the terminal cost `Qf` to the inifinite horizon LQR cost 
"""
function set_terminal_cost!(mpc)
    if mpc.settings.reference_tracking
        @warn "LQR cost not valid for reference tracking problems. Instead, use set_objective! to set Qf"
        return false
    end
    Qfx, _, _ = ared(mpc.model.F, mpc.model.G, mpc.weights.R, mpc.model.C'*mpc.weights.Q*mpc.model.C) # solve Riccati
    mpc.weights.Qfx .= Qfx
    mpc.mpqp_issetup = false
end

"""
    set_prestabilizing_feedback!(mpc,K)

Sets the prestabilizing feedback `K`
"""
function set_prestabilizing_feedback!(mpc,K::AbstractMatrix)
    mpc.K = K
    mpc.mpqp_issetup = false
end

"""
    set_prestabilizing_feedback!(mpc)

Sets the prestabilizing feedback `K` to the infinte horizon LQR gain`
"""
function set_prestabilizing_feedback!(mpc)
    _, _,mpc.K,_ = ared(mpc.model.F, mpc.model.G, mpc.weights.R+mpc.weights.Rr, mpc.model.C'*mpc.weights.Q*mpc.model.C) # solve Ricatti
    mpc.mpqp_issetup = false
end

"""
    move_block!(mpc,block)

Reduce the number of controls by keeping it constant in blocks.
For example, `block`=[2,1,3] keeps the control constant for 2 time-steps, 1 time step, and 3 time steps.
* if sum(block) ≠ mpc.Np, the resulting block will be padded or clipped
* if `block` is an Int, a vector with constant block size is created
"""
function move_block!(mpc,block)
    block = Int.(copy(block))
    if block isa Number
        block = block <= 0  ? Int[] : fill(block,mpc.Np ÷ block +1)
    end
    if isempty(block)
        mpc.move_blocks = Int[]
        mpc.Nc = mpc.Np
        mpc.mpqp_issetup = false
        return
    end
    Nnew = sum(block)
    if Nnew == mpc.Np
        mpc.move_blocks = block
    elseif(Nnew < mpc.Np) # pad
        mpc.move_blocks = block
        mpc.move_blocks[end] += mpc.Np-Nnew
    elseif Nnew > mpc.Np # clip
        tot,i = 0,1
        while((tot+=block[i]) < mpc.Np) i += 1 end
        mpc.move_blocks = block[1:i]
        mpc.move_blocks[end] += mpc.Np-tot;
    end

    mpc.Nc = sum(mpc.move_blocks[1:end-1])+1;
    mpc.mpqp_issetup = false
end

"""
    set_labels!(mpc;x,u,y,d)
Sets the name of the states `x`, controls `u`, output `u`, disturbance `d` 
"""
function set_labels!(mpc;x=nothing,u=nothing,y=nothing,d=nothing)
    isnothing(x) || (mpc.model.labels.x[:] = x)
    isnothing(u) || (mpc.model.labels.u[:] = u)
    isnothing(y) || (mpc.model.labels.y[:] = y)
    isnothing(d) || (mpc.model.labels.d[:] = d)
end

"""
    set_horizon!(mpc,Np)
Sets the prediction horizon `Np`
"""
function set_horizon!(mpc,Np, Nc = Np)
    mpc.Np = Np
    mpc.Nc = Nc
    mpc.mpqp_issetup = false
end
"""
    set_binary_controls!(mpc,bin_ids)

Makes the controls in bin_ids to binary controls 
"""
function set_binary_controls!(mpc,bin_ids)
    mpc.binary_controls = Int.(copy(bin_ids))
    mpc.mpqp_issetup = false
end
"""
    set_disturbance!(mpc,wmin,wmax)
"""
function set_disturbance!(mpc,wmin,wmax)
    mpc.model.wmin .= wmin
    mpc.model.wmax .= wmax
    mpc.mpqp_issetup = false
end
"""
    settings!(mpc,key1=value1, key2=value2,...)
"""
function settings!(mpc::MPC;kwargs...)
    settings!(mpc,kwargs)
    for (key,val) in kwargs
        key = Symbol(key)
        if hasproperty(mpc.settings,key)
            setproperty!(mpc.settings,key,val)
        else
            @warn("The setting \"$key\" does not exist")
        end
    end
end
function settings!(mpc::MPC, dict)
    for (key,val) in dict
        key = Symbol(key)
        if hasproperty(mpc.settings,key)
            setproperty!(mpc.settings,key,val)
        else
            @warn("The setting \"$key\" does not exist")
        end
    end
end

"""
    set_state_observer!(mpc;F,G,C,Q,R,x0)
Creates a steady-state Kalman filter for estimating the sate.
If `F`,`G`, and `C` are not provided, the model used in `mpc` is used in the filter
"""
function set_state_observer!(mpc::Union{MPC,ExplicitMPC};
        F=nothing,G=nothing,C=nothing,Q=nothing,R=nothing,x0=nothing)
    F = isnothing(F) ? mpc.model.F : F
    G = isnothing(G) ? mpc.model.G : G
    C = isnothing(C) ? mpc.model.C : C
    mpc.state_observer = KalmanFilter(F,G,C;Q,R,x0)
end
