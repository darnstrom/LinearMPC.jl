"""
    setup!(mpc)

Sets up the `mpc` given its current parameters and settings  
Internally, this means generating an mpQP, and setting up a DAQP workspace.
"""
function setup!(mpc::MPC)
    mpc.mpqp_issetup = false  # Reset so get_parameter_dims computes from settings
    mpc.mpQP = mpc2mpqp(mpc)
    bu,bl = mpc.mpQP.bu[:],mpc.mpQP.bl[:]
    setup_flag,_ = DAQP.setup(mpc.opt_model, mpc.mpQP.H,mpc.mpQP.f[:],mpc.mpQP.A,bu,bl,
                              mpc.mpQP.senses;break_points=mpc.mpQP.break_points, 
                              is_avi=!mpc.mpQP.is_symmetric)
    if(setup_flag < 0)
        if setup_flag == -1
            @warn " Cannot setup optimization problem - Problem is infeasible"
        elseif setup_flag == -6
            @warn " Cannot setup optimization problem - Equality constraints overdetermined"
        elseif setup_flag == -5
            @warn " Cannot setup optimization problem - Convonvex objective"
        else
            @warn " Cannot setup optimization problem " setup_flag
        end
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
        Ax, Au, Ar, Aw, Ad, Aup, Ap,
        ub, lb, ks, soft, binary, prio)
    add_constraint!(mpc;Ax,Au,ub,lb,
                    ks, soft, binary,prio)

Adds the constraints lb ≤ Ax xₖ + Au uₖ ≤ ub for the time steps k ∈ ks
(additional terms Ar rₖ, Aw wₖ, Ad dₖ, Aup u⁻ₖ, Ap pₖ are possible)

* `soft` marks if the constraint should be softened (default false)
* `binary` marks if either the upper or lower bounds should be enforced with equality (default false)
* `prio` marks the relative priority of the constraint (default 0)
"""
function add_constraint!(mpc::MPC;
        Ax = nothing, Au= nothing, Ar = zeros(0,0), Aw = zeros(0,0), Ad = zeros(0,0), Aup = zeros(0,0), Ap = zeros(0,0),
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

    push!(mpc.constraints,Constraint(Au,Ax,Ar,Aw,Ad,Aup,Ap,ub,lb,ks,soft,binary,prio))
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
    lb = !isempty(ymin) ? ymin-mpc.model.h_offset : zeros(0)
    ub = !isempty(ymax) ? ymax-mpc.model.h_offset : zeros(0)
    add_constraint!(mpc, Ax = mpc.model.C, Ad = mpc.model.Dd, lb = lb, ub = ub; ks,soft,binary,prio)
end

"""
    set_bounds!(mpc;umin,umax,ymin,ymax)

Sets the bounds umin ≤ u ≤ umax and ymin ≤ y ≤ ymax
"""
function set_bounds!(mpc::MPC; umin=zeros(0), umax=zeros(0), ymin = zeros(0), ymax = zeros(0))
    (!isempty(umin) ||  !isempty(umax)) && set_input_bounds!(mpc;umin,umax)
    (!isempty(ymin) ||  !isempty(ymax)) && set_output_bounds!(mpc;ymin,ymax)
end

"""
    set_objective!(mpc;Q,R,Rr,S,Qf,Ex,ex,Eu,eu,Sd)

Set the weights in the objective function `xN' C' Qf C xN^T + ∑ (C xₖ - rₖ)' Q (C xₖ - rₖ)  + uₖ' R uₖ + Δuₖ' Rr Δuₖ + xₖ' S uₖ + dₖ' Sd uₖ + (Ex pₖ + ex)'xₖ + (Eu pₖ + eu)'uₖ

A vector is interpreted as a diagonal matrix.

`Sd` (nd × nu) is a cross term between the measurable disturbance `d` and the control. With `Sd = R` the
control penalty becomes `(uₖ + dₖ)' R (uₖ + dₖ)` up to a constant, which penalizes the control relative to
the input that cancels a disturbance entering through `Bd = B` instead of relative to zero, and thereby
removes the steady-state error a direct penalty on `u` otherwise causes under a persistent disturbance.
"""
function set_objective!(mpc::MPC;Q = zeros(0,0), R=zeros(0,0), Rr=zeros(0,0), S= zeros(0,0),Qf=zeros(0,0), Qfx=zeros(0,0),
        Ex = zeros(0,0), ex = zeros(0), Eu = zeros(0,0), eu = zeros(0), Sd = zeros(0,0))
    Qw = isempty(Q) ? copy(mpc.weights.Q) : matrixify(Q,mpc.model.ny)
    Rw = isempty(R) ? copy(mpc.weights.R) : matrixify(R,mpc.model.nu)
    Rrw = isempty(Rr) ? copy(mpc.weights.Rr) : matrixify(Rr,mpc.model.nu)
    Sw = isempty(S) ? copy(mpc.weights.S) : float(S)
    Qfw = isempty(Qf) ? copy(mpc.weights.Qf) : matrixify(Qf,mpc.model.ny)
    Qfxw = isempty(Qfx) ? copy(mpc.weights.Qfx) : matrixify(Qfx,mpc.model.nx)
    Exw = isempty(Ex) ? copy(mpc.weights.Ex) : float(Ex)
    exw = isempty(ex) ? copy(mpc.weights.ex) : float(ex)
    Euw = isempty(Eu) ? copy(mpc.weights.Eu) : float(Eu)
    euw = isempty(eu) ? copy(mpc.weights.eu) : float(eu)
    Sdw = isempty(Sd) ? copy(mpc.weights.Sd) : matrixify(Sd)
    mpc.weights = MPCWeights(Qw,Rw,Rrw,Sw,Qfw,Qfxw,Exw,exw,Euw,euw,Sdw)
    mpc.mpqp_issetup = false
end


function set_objective!(mpc::MPC, uids::Vector{Int};Q = zeros(0,0), R=zeros(0,0), 
        Rr=zeros(0,0), S= zeros(0,0), Qf=zeros(0,0), Qfx=zeros(0,0),
        Ex = zeros(0,0), ex = zeros(0), Eu = zeros(0,0), eu = zeros(0), Sd = zeros(0,0))
    nu,ny,nx = length(uids), mpc.model.ny, mpc.model.nx
    Q   = isempty(Q)   ? zeros(mpc.model.ny,mpc.model.ny) : matrixify(Q,ny)
    R   = isempty(R)   ? zeros(nu,nu) : matrixify(R,nu)
    Rr  = isempty(Rr)  ? zeros(nu,nu) : matrixify(Rr,nu)
    S   = isempty(S)   ? zeros(nx,nu) : float(S)
    Qf  = isempty(Qf)  ? copy(Q) : matrixify(Qf,ny)
    Qfx = isempty(Qfx) ? zeros(nx,nx) :  matrixify(Qfx,nx)
    Ex  = isempty(Ex)  ? zeros(nx,0) : float(Ex)
    ex  = isempty(ex)  ? zeros(nx) : float(ex)
    Eu  = isempty(Eu)  ? zeros(nu,0) : float(Eu)
    eu  = isempty(eu)  ? zeros(nu) : float(eu)
    Sd  = isempty(Sd)  ? zeros(0,nu) : matrixify(Sd)

    mpc.weights.Rr[uids,uids] .= Rr # To be able to keep track of nuprev
    push!(mpc.objectives, (MPCWeights(Q,R,Rr,S,Qf,Qfx,Ex,ex,Eu,eu,Sd),uids))
    mpc.mpqp_issetup = false
end


set_weights! = set_objective! # backwards compatibility
add_objective! = set_objective!

function empty_objectives!(mpc::MPC)
    empty!(mpc.objectives)
    mpc.mpqp_issetup = false
end

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
function move_block!(mpc,block::Nothing)
    mpc.move_blocks = Vector{Int}[]
    mpc.Nc = mpc.Np
    mpc.mpqp_issetup=false
end

function move_block!(mpc,block::Number)
    block = block <= 0 ? Int[] : fill(Int(block),mpc.Np ÷ block +1)
    move_block!(mpc,block)
end

function move_block!(mpc,block::AbstractVector{<:Number})
    isempty(block) && return move_block!(mpc,nothing)
    return move_block!(mpc,[block for _ in 1:mpc.model.nu])
end

function move_block!(mpc,blocks::Vector{<:AbstractVector{<:Number}})
    length(blocks) == mpc.model.nu || ArgumentError("Need to have blocks for every control input")
    blocks_formated = [format_move_block(mb,mpc.Np) for mb in blocks]
    any(isempty(mb) for mb in blocks_formated) && ArgumentError("One block is empty")

    mpc.move_blocks = blocks_formated 
    mpc.Nc = maximum(sum(mb[1:end-1]) for mb in mpc.move_blocks)+1
    mpc.mpqp_issetup = false
end

"""
    add_logic_constraint!(mpc; delta_ids, Adelta, ub, lb, ks, prio)

Adds inequalities involving only controls that have already been marked as
binary with `set_binary_controls!`.

This is a convenience wrapper for relations of the form
`lb ≤ Adelta*δₖ ≤ ub`. `Aδ` can be used as a shorthand alias for `Adelta`.
The same constraints could also be added directly with `add_constraint!`.
"""
function add_logic_constraint!(mpc::MPC; delta_ids, Adelta=zeros(0,0), Aδ=zeros(0,0), ub=zeros(0), lb=zeros(0), ks = 1:mpc.Np, prio = 0)
    Adelta = isempty(Adelta) ? Aδ : Adelta
    delta_ids = _validate_input_ids(delta_ids, mpc.model.nu, "delta_ids")
    ub = isempty(ub) ? zeros(size(Adelta,1)) : ub
    lb = isempty(lb) ? fill(-1e30, size(Adelta,1)) : lb
    nrows = max(length(lb), length(ub), size(Adelta, 1))
    nrows == 0 && return
    Au = _expand_input_block(Adelta, delta_ids, mpc.model.nu, nrows, "Adelta")
    ub = length(ub) == nrows ? ub : [ub; fill(1e30, nrows - length(ub))]
    lb = length(lb) == nrows ? lb : [lb; fill(-1e30, nrows - length(lb))]
    add_constraint!(mpc; Au, ub, lb, ks, prio)
end

"""
    add_indicator_constraint!(mpc, delta_id;
        Ax, Au, c, m, M, ϵ, sense, ks, prio)

Adds a big-M indicator relation between a control `u[delta_id]` and an affine
expression in the state and controller inputs. The selected input should
already be marked as binary with `set_binary_controls!`.

With `sense = :le`, the constraint encodes
`u[delta_id] = 1 ↔ Ax*xₖ + Au*uₖ + c ≤ 0`.
With `sense = :ge`, it encodes
`u[delta_id] = 1 ↔ Ax*xₖ + Au*uₖ + c ≥ 0`.

The scalars `m` and `M` are the lower and upper big-M bounds for the affine
expression, and `ϵ` is the strictness margin used in the reverse implication.
"""
function add_indicator_constraint!(mpc::MPC, delta_id::Integer;
        Ax = zeros(1, mpc.model.nx), Au = zeros(1, mpc.model.nu), c = zeros(size(Ax,1)),
        m, M, ϵ = sqrt(eps(Float64)), sense::Symbol = :le, ks = 1:mpc.Np, prio = 0)
    _validate_input_id(delta_id, mpc.model.nu, "delta_id")
    Mv = fill(float(M), size(Ax,1))
    mv = fill(float(m), size(Ax,1))
    A1 = copy(Au)
    A2 = copy(Au)
    if sense == :le
        A1[:, delta_id] .+= Mv
        A2 = -copy(Au)
        A2[:, delta_id] .+= mv .- ϵ
        add_constraint!(mpc; Ax, Au = A1, ub = Mv .- c, lb = fill(-1e30, size(Ax,1)), ks, prio)
        add_constraint!(mpc; Ax = -Ax, Au = A2, ub = fill(-ϵ, size(Ax,1)) .- c, lb = fill(-1e30, size(Ax,1)), ks, prio)
    elseif sense == :ge
        A1[:, delta_id] .-= Mv
        A2 = -copy(Au)
        A2[:, delta_id] .+= mv .- ϵ
        add_constraint!(mpc; Ax, Au = A1, ub = -c, lb = fill(-1e30, size(Ax,1)), ks, prio)
        add_constraint!(mpc; Ax = -Ax, Au = A2, ub = c .- mv, lb = fill(-1e30, size(Ax,1)), ks, prio)
    end
end

"""
    add_product_constraint!(mpc, z_id, delta_id;
        Ax, Au, c, m, M, ks, prio)

Adds a mixed-integer reformulation of the product
`u[z_id] = u[delta_id] * (Ax*xₖ + Au*uₖ + c)`, where `u[delta_id]` is expected
to be one of the controls marked as binary with `set_binary_controls!`.

The scalars `m` and `M` must bound the affine factor over the relevant domain.
"""
function add_product_constraint!(mpc::MPC, z_id::Integer, delta_id::Integer;
        Ax = zeros(1, mpc.model.nx), Au = zeros(1, mpc.model.nu), c = zeros(size(Ax,1)),
        m, M, ks = 1:mpc.Np, prio = 0)
    _validate_input_id(delta_id, mpc.model.nu, "delta_id")
    _validate_input_id(z_id, mpc.model.nu, "z_id")
    Mv = fill(float(M), size(Ax,1))
    mv = fill(float(m), size(Ax,1))

    A1 = zeros(size(Ax,1), mpc.model.nu)
    A1[:, delta_id] .-= Mv
    A1[:, z_id] .+= 1.0
    add_constraint!(mpc; Au = A1, ub = zeros(size(Ax,1)), lb = fill(-1e30, size(Ax,1)), ks, prio)

    A2 = zeros(size(Ax,1), mpc.model.nu)
    A2[:, delta_id] .+= mv
    A2[:, z_id] .-= 1.0
    add_constraint!(mpc; Au = A2, ub = zeros(size(Ax,1)), lb = fill(-1e30, size(Ax,1)), ks, prio)

    A3 = -copy(Au)
    A3[:, delta_id] .-= mv
    A3[:, z_id] .+= 1.0
    add_constraint!(mpc; Ax = -Ax, Au = A3, ub = c .- mv, lb = fill(-1e30, size(Ax,1)), ks, prio)

    A4 = copy(Au)
    A4[:, delta_id] .+= Mv
    A4[:, z_id] .-= 1.0
    add_constraint!(mpc; Ax, Au = A4, ub = Mv .- c, lb = fill(-1e30, size(Ax,1)), ks, prio)
end

function add_ifthenelse_relation!(mpc::MPC, Au_out, delta_id::Integer;
        Ax_then, Au_then, c_then, Ax_else, Au_else, c_else,
        m_then, M_then, m_else, M_else, ks = 1:mpc.Np, prio = 0)
    _validate_input_id(delta_id, mpc.model.nu, "delta_id")
    nrows = size(Ax_then, 1)

    M1 = fill(float(M_then), nrows)
    m1 = fill(float(m_then), nrows)
    M2 = fill(float(M_else), nrows)
    m2 = fill(float(m_else), nrows)

    A1 = Au_out - Au_else
    A1[:, delta_id] .+= m2 .- M1
    add_constraint!(mpc; Ax = -Ax_else, Au = A1, ub = c_else, lb = fill(-1e30, nrows), ks, prio)

    A2 = Au_else - Au_out
    A2[:, delta_id] .+= m1 .- M2
    add_constraint!(mpc; Ax = Ax_else, Au = A2, ub = -c_else, lb = fill(-1e30, nrows), ks, prio)

    A3 = Au_out - Au_then
    A3[:, delta_id] .+= M2 .- m1
    add_constraint!(mpc; Ax = -Ax_then, Au = A3, ub = c_then .+ (M2 .- m1), lb = fill(-1e30, nrows), ks, prio)

    A4 = Au_then - Au_out
    A4[:, delta_id] .+= M1 .- m2
    add_constraint!(mpc; Ax = Ax_then, Au = A4, ub = -c_then .+ (M1 .- m2), lb = fill(-1e30, nrows), ks, prio)
end

"""
    add_ifthenelse_constraint!(mpc, z_id, delta_id;
        Ax_then, Au_then, c_then, Ax_else, Au_else, c_else,
        m_then, M_then, m_else, M_else, ks, prio)

Adds a mixed-integer `if/then/else` relation for an auxiliary input `u[z_id]`.

If `u[delta_id] = 1`, then
`u[z_id] = Ax_then*xₖ + Au_then*uₖ + c_then`.
Otherwise,
`u[z_id] = Ax_else*xₖ + Au_else*uₖ + c_else`.

The `m_*` and `M_*` arguments bound the corresponding affine branch values.
"""
function add_ifthenelse_constraint!(mpc::MPC, z_id::Integer, delta_id::Integer;
        Ax_then = zeros(1, mpc.model.nx), Au_then = zeros(1, mpc.model.nu), c_then = zeros(size(Ax_then,1)),
        Ax_else = zeros(size(Ax_then,1), mpc.model.nx), Au_else = zeros(size(Ax_then,1), mpc.model.nu), c_else = zeros(size(Ax_then,1)),
        m_then, M_then, m_else, M_else, ks = 1:mpc.Np, prio = 0)
    _validate_input_id(z_id, mpc.model.nu, "z_id")
    Au_out = zeros(size(Ax_then,1), mpc.model.nu)
    Au_out[:, z_id] .= 1.0
    add_ifthenelse_relation!(mpc, Au_out, delta_id;
                             Ax_then, Au_then, c_then, Ax_else, Au_else, c_else,
                             m_then, M_then, m_else, M_else, ks, prio)
end

"""
    add_ifthenelse_input_constraint!(mpc, u_id, delta_id;
        Ax_then, Au_then, c_then, Ax_else, Au_else, c_else,
        m_then, M_then, m_else, M_else, ks, prio)

Adds a mixed-integer `if/then/else` relation for an input `u[u_id]`.

If `u[delta_id] = 1`, then
`u[u_id] = Ax_then*xₖ + Au_then*uₖ + c_then`.
Otherwise,
`u[u_id] = Ax_else*xₖ + Au_else*uₖ + c_else`.

The `m_*` and `M_*` arguments bound the corresponding affine branch values.
"""
function add_ifthenelse_input_constraint!(mpc::MPC, u_id::Integer, delta_id::Integer;
        Ax_then = zeros(1, mpc.model.nx), Au_then = zeros(1, mpc.model.nu), c_then = zeros(size(Ax_then,1)),
        Ax_else = zeros(size(Ax_then,1), mpc.model.nx), Au_else = zeros(size(Ax_then,1), mpc.model.nu), c_else = zeros(size(Ax_then,1)),
        m_then, M_then, m_else, M_else, ks = 1:mpc.Np, prio = 0)
    _validate_input_id(u_id, mpc.model.nu, "u_id")
    Au_out = zeros(size(Ax_then,1), mpc.model.nu)
    Au_out[:, u_id] .= 1.0
    add_ifthenelse_relation!(mpc, Au_out, delta_id;
                             Ax_then, Au_then, c_then, Ax_else, Au_else, c_else,
                             m_then, M_then, m_else, M_else, ks, prio)
end

function format_move_block(block::AbstractVector{<:Number},Np::Int)
    block = Int.(copy(block))
    isempty(block) && return Int[]
    Nnew = sum(block)
    if(Nnew < Np) # pad
        block[end] += Np-Nnew
    elseif Nnew > Np # clip
        tot,i = 0,1
        while((tot+=block[i]) < Np) i += 1 end
        block = block[1:i]
        block[end] += Np-tot;
    end
    return block
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
function set_horizon!(mpc,Np, Nc = Np, Nc_binary = mpc.Nc_binary)
    mpc.Np = Np
    mpc.Nc = Nc
    mpc.Nc_binary= Nc_binary
    mpc.mpqp_issetup = false
end
"""
    set_binary_controls!(mpc,bin_ids, Nc_binary=nothing)

Makes the controls in bin_ids to binary controls.
Nc_binary is the "binary control horizon" (default = control horizon) 
"""
function set_binary_controls!(mpc,bin_ids,Nc_binary=-1)
    mpc.binary_controls = Int.(copy(bin_ids))
    mpc.Nc_binary = Nc_binary
    mpc.mpqp_issetup = false
end

function _validate_input_ids(ids, nu::Int, name::AbstractString)
    ids = Int.(collect(ids))
    isempty(ids) && return ids
    return ids
end

function _expand_input_block(Ablock, ids, nu::Int, nrows::Int, name::AbstractString)
    ids = _validate_input_ids(ids, nu, name)
    isempty(Ablock) && return zeros(nrows, nu)
    Ablock = float(Ablock)
    if isempty(ids)
        return Ablock
    end
    A = zeros(nrows, nu)
    A[:, ids] .= Ablock
    return A
end

function _validate_input_id(id::Integer, nu::Int, name::AbstractString)
    1 <= id <= nu || throw(ArgumentError("$name must be between 1 and $nu"))
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
    set_x0_uncertainty!(mpc,wmin,wmax)
"""
function set_x0_uncertainty!(mpc,x0_uncertainty)
    mpc.Δx0 .= x0_uncertainty 
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
            mpc.mpqp_issetup = false
        else
            @warn("The setting \"$key\" does not exist")
        end
    end
end

"""
    set_state_observer!(mpc;F,G,Gd,C,Dd,Q,R,x0)
Creates a steady-state Kalman filter for estimating the sate.
If `F`,`G`, and `C` are not provided, the model used in `mpc` is used in the filter
"""
function set_state_observer!(mpc::Union{MPC,ExplicitMPC};
        F=nothing,G=nothing,Gd=nothing,C=nothing,Dd=nothing,
        f_offset=nothing, h_offset=nothing,
        Q=nothing,R=nothing,x0=nothing)
    F = isnothing(F) ? mpc.model.F : F
    G = isnothing(G) ? mpc.model.G : G
    Gd = isnothing(Gd) ? mpc.model.Gd : Gd
    C = isnothing(C) ? mpc.model.C : C
    Dd = isnothing(Dd) ? mpc.model.Dd : Dd
    f_offset = isnothing(f_offset) ? mpc.model.f_offset : f_offset 
    h_offset = isnothing(h_offset) ? mpc.model.h_offset : h_offset 
    mpc.state_observer = KalmanFilter(F,G,C;Gd,Dd,f_offset,h_offset,Q,R,x0)
end

function normalize_offset_free_method(method::Symbol)
    aliases = Dict(
        :state => :state_disturbance,
        :state_disturbance => :state_disturbance,
        :velocity => :velocity,
        :output => :output_disturbance,
        :output_disturbance => :output_disturbance,
        :general => :general,
    )
    haskey(aliases, method) || throw(ArgumentError("Unknown offset-free method $method"))
    return aliases[method]
end

function rebuild_model(model::Model, Gd, Dd, disturbance_labels)
    labels = Labels(model.labels.x, model.labels.u, model.labels.y, Symbol.(disturbance_labels))
    nd = size(Gd, 2)
    return Model(model.F, model.G, float(Gd), model.f_offset, model.xo, model.uo,
                 model.wmin, model.wmax, model.C, float(Dd), model.h_offset,
                 model.true_dynamics, model.true_h,
                 model.nx, model.nu, model.ny, nd, model.Ts, labels)
end

function strip_offset_free_model(model::Model, nd_measured::Int)
    return rebuild_model(model, model.Gd[:,1:nd_measured], model.Dd[:,1:nd_measured], model.labels.d[1:nd_measured])
end

function append_offset_free_model(model::Model, Bd, Cd, disturbance_labels)
    return rebuild_model(model, [model.Gd Bd], [model.Dd Cd], [model.labels.d; disturbance_labels])
end

function default_offset_free_labels(method::Symbol, nd::Int)
    prefix = method == :output_disturbance ? "yoff" : "dof"
    return Symbol.(prefix .* string.(1:nd))
end

function nominal_observer_gain(F, C; Q=nothing, R=nothing)
    nx, ny = size(F, 1), size(C, 1)
    return KalmanFilter(F, zeros(nx, ny), C; Q, R).K
end

function validate_offset_free_model(F, C, Bd, Cd)
    nx = size(F, 1)
    nd = size(Bd, 2)
    ny = size(C, 1)
    size(Bd, 1) == nx || throw(ArgumentError("Bd must have $nx rows"))
    size(Cd) == (ny, nd) || throw(ArgumentError("Cd must have size ($ny, $nd)"))
    rank([F - Matrix{Float64}(I, nx, nx) Bd; C Cd]) == nx + nd ||
        throw(ArgumentError("Offset-free disturbance model violates rank([F-I Bd; C Cd]) = nx + nd"))
end

function build_offset_free_observer(model::Model, nd_measured::Int, method::Symbol;
        Q=nothing, R=nothing, K=nothing, Bd=nothing, Cd=nothing,
        Kx=nothing, Kd=nothing, x0=nothing, d0=nothing)

    F, G, C = model.F, model.G, model.C
    method = normalize_offset_free_method(method)
    nx, ny = model.nx, model.ny

    if method == :state_disturbance || method == :velocity
        K = isnothing(K) ? nominal_observer_gain(F, C; Q, R) : float(K)
        size(K) == (nx, ny) || throw(ArgumentError("K must have size ($nx, $ny)"))
        Bd = K
        Cd = Matrix{Float64}(I, ny, ny) - C*K
        Kx = K
        Kd = Matrix{Float64}(I, ny, ny)
    elseif method == :output_disturbance
        Bd = zeros(nx, ny)
        Cd = Matrix{Float64}(I, ny, ny)
    else
        isnothing(Bd) && throw(ArgumentError("Method :general requires Bd"))
        isnothing(Cd) && throw(ArgumentError("Method :general requires Cd"))
        Bd = float(Bd)
        Cd = float(Cd)
    end

    Bd = float(Bd)
    Cd = float(Cd)
    validate_offset_free_model(F, C, Bd, Cd)
    ndo = size(Bd, 2)

    x0 = isnothing(x0) ? zeros(nx) : float(x0)
    d0 = isnothing(d0) ? zeros(ndo) : float(d0)
    length(x0) == nx || throw(ArgumentError("x0 must have length $nx"))
    length(d0) == ndo || throw(ArgumentError("d0 must have length $ndo"))

    Faug = [F Bd; zeros(ndo, nx) Matrix{Float64}(I, ndo, ndo)]
    Gaug = [G; zeros(ndo, model.nu)]
    Gdaug = [model.Gd[:,1:nd_measured]; zeros(ndo, nd_measured)]
    Caug = [C Cd]
    xaug0 = [x0; d0]
    faug = [model.f_offset; zeros(ndo)]

    estimator = if !isnothing(Kx) || !isnothing(Kd) || method == :state_disturbance || method == :velocity
        Kx = isnothing(Kx) ? zeros(nx, ny) : float(Kx)
        Kd = isnothing(Kd) ? zeros(ndo, ny) : float(Kd)
        size(Kx) == (nx, ny) || throw(ArgumentError("Kx must have size ($nx, $ny)"))
        size(Kd) == (ndo, ny) || throw(ArgumentError("Kd must have size ($ndo, $ny)"))
        KalmanFilter(Faug, Gaug, Gdaug, faug, Caug, model.Dd[:,1:nd_measured], model.h_offset,
                     [Kx; Kd], xaug0)
    else
        KalmanFilter(Faug, Gaug, Caug; Gd=Gdaug, Dd=model.Dd[:,1:nd_measured],
                     f_offset=faug, h_offset=model.h_offset, Q, R, x0=xaug0)
    end

    return OffsetFreeObserver(estimator, model.C, model.Dd[:,1:nd_measured], model.h_offset,
                              nx, nd_measured, ndo, method), Bd, Cd
end

"""
    set_offset_free_observer!(mpc; method=:state_disturbance, Q, R, K, Bd, Cd, Kx, Kd, x0, d0, disturbance_labels)

Create an offset-free observer/controller pair following the formulations reviewed by
Pannocchia (2015). The controller model is augmented with constant disturbance channels,
and the observer estimates both the nominal state and the disturbance.

Supported methods are:
- `:state_disturbance` (default), using the equivalent disturbance-model realization from Theorem 9
- `:velocity`, using the equivalent disturbance-model realization from Theorem 15 with `Ke = I`
- `:output_disturbance`, using a pure output-bias disturbance model
- `:general`, using user-provided `Bd` and `Cd`

For `:state_disturbance` and `:velocity`, the nominal observer gain `K` can be provided
directly or obtained from the existing steady-state Kalman filter tuning through `Q` and `R`.
"""
function set_offset_free_observer!(mpc::MPC;
        method::Symbol=:state_disturbance,
        Q=nothing, R=nothing, K=nothing,
        Bd=nothing, Cd=nothing,
        Kx=nothing, Kd=nothing,
        x0=nothing, d0=nothing,
        disturbance_labels=nothing)

    nd_measured = mpc.state_observer isa OffsetFreeObserver ? mpc.state_observer.nd_measured : mpc.model.nd
    mpc.model = strip_offset_free_model(mpc.model, nd_measured)

    observer, Bd, Cd = build_offset_free_observer(mpc.model, nd_measured, method;
                                                  Q, R, K, Bd, Cd, Kx, Kd, x0, d0)

    labels = isnothing(disturbance_labels) ? default_offset_free_labels(observer.formulation, size(Bd, 2)) : disturbance_labels
    length(labels) == size(Bd, 2) || throw(ArgumentError("Need $(size(Bd, 2)) disturbance labels"))

    mpc.model = append_offset_free_model(mpc.model, Bd, Cd, labels)
    mpc.state_observer = observer
    mpc.mpqp_issetup = false
    return observer
end

"""
    set_operating_point!(mpc;xo,uo)
Sets the operating point to the state xo and control uo and linearize
"""
function set_operating_point!(mpc;xo=nothing,uo=nothing,relinearize=true)
    !isnothing(xo) && (mpc.model.xo[:] = xo)
    !isnothing(uo) && (mpc.model.uo[:] = uo)

    if !isnothing(xo) || !isnothing(uo)
        mpc.model = LinearMPC.Model(mpc.model.true_dynamics,mpc.model.true_h,
                                    mpc.model.xo,mpc.model.uo)
        mpc.mpqp_issetup = false
    end
end


"""
    set_offset!(mpc;xo,uo,doff,fo,ho)
Set bias terms in dynamics and measurements.

Concretely we have that
`f_offet = fo - F * xo - G * uo - Gd * doff` and
`h_offet = ho - C * xo - Dd * doff`,
which adds a constant term to the dynamics and measurement function, respectively.
Note that if the system is linearized, these offsets are set automatically.
If some of the offset are not entered, they are interpreted as zero.
"""
function set_offset!(mpc;xo=zeros(0),uo=zeros(0),doff=zeros(0),fo=zeros(0),ho=zeros(0))
    isempty(xo) && (xo = zeros(mpc.model.nx))
    isempty(uo) && (uo= zeros(mpc.model.nu))
    isempty(fo) && (fo = zeros(mpc.model.nx))
    isempty(ho) && (ho = zeros(mpc.model.ny))
    isempty(doff) && (doff = zeros(mpc.model.nd))

    mpc.model.xo .= xo
    mpc.model.uo .= uo
    mpc.uprev .= uo

    mpc.model.f_offset .= fo - mpc.model.F*xo - mpc.model.G*uo - mpc.model.Gd*doff
    mpc.model.h_offset .= ho - mpc.model.C*xo - mpc.model.Dd*doff

    mpc.mpqp_issetup=false
end
