"""
    setup!(mpc)

Sets up the `mpc` given its current parameters and settings  
Internally, this means generating an mpQP, and setting up a DAQP workspace.
"""
function setup!(mpc::MPC)
    mpc.mpQP = mpc2mpqp(mpc)
    mpc.opt_model = DAQP.Model()
    if(mpc.settings.QP_double_sided)
        bu,bl = mpc.mpQP.bu[:],mpc.mpQP.bl[:]
    else
        bu,bl = mpc.mpQP.b[:], -1e30*ones(length(mpc.mpQP.b))
    end
    DAQP.setup(mpc.opt_model, mpc.mpQP.H,mpc.mpQP.f[:],mpc.mpQP.A,bu,bl,mpc.mpQP.senses)
end

"""
    set_bounds!(mpc;umin,umax)

Sets the bounds umin ≤ u ≤ umax 
"""
function set_bounds!(mpc::MPC; umin=zeros(0), umax=zeros(0))
    nmin,nmax = length(umin),length(umax)
    nb = max(nmin,nmax)
    nb == 0 && return
    nb != mpc.model.nu  && @error("# of controls are $(mpc.model.nu), got bounds of dimension $nb")

    mpc.umin = [umin;-1e30*ones(nb-nmin)]
    mpc.umax = [umax;+1e30*ones(nb-nmax)]
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
    set_weights!(mpc;Q,R,Rr,S,rho,Qf)

The weights in the objective function `xN' Qf xN^T + ∑ (C xₖ - rₖ)' Q (C xₖ - rₖ)  + uₖ' R uₖ + Δuₖ' Rr Δuₖ + xₖ' S uₖ 

A vector is interpreted as a diagonal matrix.
* `rho` Is an additional weight for the soft constraints (default value: 1e6)
"""
function set_weights!(mpc::MPC;Q = nothing,R=nothing,Rr=nothing, S= nothing, rho=nothing, Qf=nothing)
    if !isnothing(Q)
        mpc.weights.Q = matrixify(Q,mpc.model.ny)
    end
    if !isnothing(R)
        mpc.weights.R = matrixify(R,mpc.model.nu)
    end
    if !isnothing(Rr)
        mpc.weights.Rr = matrixify(Rr,mpc.model.nu)
    end
    if !isnothing(S)
        mpc.weights.S = float(S)
    end
    if !isnothing(rho)
        mpc.weights.rho = rho
    end
    if !isnothing(Qf)
        mpc.weights.Qf = matrixify(Qf,mpc.model.ny)
    end
end

# Terminal ingredients
using MatrixEquations 

"""
    set_terminal_cost!(mpc)

Sets the terminal cost `Qf` to the inifinite horizon LQR cost 
"""
function set_terminal_cost!(mpc)
    Qf, _, _ = ared(mpc.model.F, mpc.model.G, mpc.weights.R, mpc.C'*mpc.weights.Q*mpc.C) # solve Riccati
    mpc.weights.Qf = Qf
end

"""
    set_prestabilizing_feedback!(mpc,K)

Sets the prestabilizing feedback `K`
"""
function set_prestabilizing_feedback!(mpc,K::AbstractMatrix)
    mpc.K = K
end

"""
    set_prestabilizing_feedback!(mpc)

Sets the prestabilizing feedback `K` to the infinte horizon LQR gain`
"""
function set_prestabilizing_feedback!(mpc)
    _, _,mpc.K,_ = ared(mpc.model.F, mpc.model.G, mpc.weights.R+mpc.weights.Rr, mpc.C'*mpc.weights.Q*mpc.C) # solve Ricatti
end

"""
    move_block!(mpc,block)

Reduce the number of controls by keeping it constant in blocks.
For example, `block`=[2,1,3] keeps the control constant for 2 time-steps, 1 time step, and 3 time steps.
* if sum(block) ≠ mpc.Nc, the resulting block will be padded or clipped
* if `block` is an Int, a vector with constant block size is created
"""
function move_block!(mpc,block::Vector{Int})
    Nnew = sum(block)
    if Nnew == mpc.Nc
        mpc.move_blocks[:] = block
    elseif(Nnew < mpc.Nc) # pad
        mpc.move_blocks[:] = block
        mpc.move_blocks[end] += mpc.Nc-Nnew
    elseif Nnew > mpc.Nc # clip
        tot,i = 0,1
        while((tot+=block[i]) < mpc.Nc) i += 1 end
        mpc.move_blocks = block[1:i]
        mpc.move_blocks[end] += mpc.Nc-tot;
    end
end
function move_block!(mpc,block::Int)
    nb,res  = mpc.Nc ÷ block, mpc.Nc % block
    mpc.move_blocks = fill(block,nb+1)
    mpc.move_blocks[end] = res
end

"""
    set_labels!(mpc;x,u,y,d)
Sets the name of the states `x`, controls `u`, output `u`, disturbance `d` 
"""
function set_labels!(mpc;x=nothing,u=nothing,y=nothing,d=nothing)
    isnothing(x) || (mpc.labels.x[:] = x)
    isnothing(u) || (mpc.labels.u[:] = u)
    isnothing(y) || (mpc.labels.y[:] = y)
    isnothing(d) || (mpc.labels.d[:] = d)
end
