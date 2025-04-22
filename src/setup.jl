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
    nb != mpc.nu  && @error("# of controls are $(mpc.nu), got bounds of dimension $nb")

    mpc.umin = [umin;-1e30*ones(nb-nmin)]
    mpc.umax = [umax;+1e30*ones(nb-nmax)]
end

"""
    add_constraint!(mpc;Ax,Au,ub,lb,
                    ks, soft, binary,prio)

Adds the constraints lb ≤ Au xₖ + Au uₖ ≤ ub for the time steps k ∈ ks 

* `soft` marks if the constraint should be softened (default false)
* `binary` marks if either the upper or lower bounds should be enforced with equality (default false)
* `prio` marks the relative priority of the constraint (default 0)
"""
function add_constraint!(mpc::MPC; Ax = nothing, Au= nothing, ub = zeros(0), lb = zeros(0), 
        ks = 2:mpc.Np, soft=false, binary=false, prio = 0)
    if isnothing(Ax) && isnothing(Bx)
        return
    end

    # Get length of constraint
    nlb,nub = length(lb),length(ub) 
    m = max(nlb,nub)
    m == 0 && return

    ub = nub == m ? ub : [ub;1e30*ones(m-nub)]
    lb = nlb == m ? lb : [lb;-1e30*ones(m-nlb)]

    Ax = isnothing(Ax) ? zeros(m,mpc.nx) : Ax
    Au = isnothing(Au) ? zeros(m,mpc.nu) : Au

    push!(mpc.constraints,Constraint(Au,Ax,ub,lb,ks,soft,binary,prio))

end

"""
    set_output_bounds!(mpc;ymin,ymax,
                    ks, soft, binary,prio)

Adds the constraints lb ≤ C x  ≤ ub for the time steps k ∈ ks 

* `soft` marks if the constraint should be softened (default false)
* `binary` marks if either the upper or lower bounds should be enforced with equality (default false)
* `prio` marks the relative priority of the constraint (default 0)
"""
function set_output_bounds!(mpc::MPC; ymin=nothing, ymax=nothing, ks = nothing, soft = true, binary=false, prio = 0)
    add_constraint!(mpc, Ax = mpc.C, lb = ymin, ub = ymax;ks,soft,binary,prio)
end

"""
    set_weights!(mpc;Q,R,Rr,S,rho,Qf)

The weights in the objective function `xN' Qf xN^T + ∑ (C xₖ - rₖ)' Q (C xₖ - rₖ)  + uₖ' R uₖ + Δuₖ' Rr Δuₖ + xₖ' S uₖ 

A vector is interpreted as a diagonal matrix.
* `rho` Is an additional weight for the soft constraints (default value: 1e6)
"""
function set_weights!(mpc::MPC;Q = nothing,R=nothing,Rr=nothing, S= nothing, rho=nothing, Qf=nothing)
    if !isnothing(Q)
        mpc.weights.Q = Q isa AbstractMatrix ? Q : diagm(Q) 
    end
    if !isnothing(R)
        mpc.weights.R = R isa AbstractMatrix ? R : diagm(R) 
    end
    if !isnothing(Rr)
        mpc.weights.Rr = Rr isa AbstractMatrix ? Rr : diagm(Rr) 
    end
    if !isnothing(S)
        mpc.weights.S = S
    end
    if !isnothing(rho)
        mpc.weights.rho = rho
    end
    if !isnothing(Qf)
        mpc.weights.Qf = Qf isa AbstractMatrix ? Qf : diagm(Qf) 
    end
end

# Terminal ingredients
using MatrixEquations 

"""
    set_terminal_cost!(mpc)

Sets the terminal cost `Qf` to the inifinite horizon LQR cost 
"""
function set_terminal_cost!(mpc)
    Qf, _, _ = ared(mpc.F, mpc.G, mpc.weights.R, mpc.C'*mpc.weights.Q*mpc.C) # solve Riccati
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
    _, _,mpc.K,_ = ared(mpc.F, mpc.G, mpc.weights.R+mpc.weights.Rr, mpc.C'*mpc.weights.Q*mpc.C) # solve Ricatti
end
