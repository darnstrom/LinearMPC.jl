function set_bounds!(mpc::MPC; umin=zeros(0), umax=zeros(0))
    nmin,nmax = length(umin),length(umax)
    nb = max(nmin,nmax)
    nb == 0 && return
    nb != mpc.nu  && @error("# of controls are $(mpc.nu), got bounds of dimension $nb")

    mpc.umin = [umin;-1e30*ones(nb-nmin)]
    mpc.umax = [umax;+1e30*ones(nb-nmax)]
end


function add_constraint!(mpc::MPC; Ax = nothing, Au= nothing, ub = zeros(0), lb = zeros(0), 
        ks = 2:mpc.Np, soft=false, binary=false)
    if isnothing(Ax) && isnothing(Bx)
        return
    end

    # Get length of constraint
    nlb,nub  =length(lb),length(ub) 
    m = max(nlb,nub)
    m == 0 && return

    ub = nub == m ? ub : [ub;1e30*ones(m-nub)]
    lb = nlb == m ? lb : [lb;-1e30*ones(m-nub)]

    Ax = isnothing(Ax) ? zeros(m,mpc.nx) : Ax
    Au = isnothing(Au) ? zeros(m,mpc.nu) : Au

    push!(mpc.constraints,Constraint(Au,Ax,ub,lb,ks,soft,binary))

end

function set_output_bounds!(mpc::MPC; ymin=nothing, ymax=nothing, ks = nothing, soft = true, binary=false)
    add_constraint!(mpc, Ax = mpc.C, lb = ymin, ub = ymax;ks,soft,binary)
end

function set_weights!(mpc::MPC;Q = nothing,R=nothing,Rr=nothing,rho=nothing, Qf=nothing)
    if !isnothing(Q)
        mpc.weights.Q = Q isa AbstractMatrix ? Q : diagm(Q) 
    end
    if !isnothing(R)
        mpc.weights.R = R isa AbstractMatrix ? R : diagm(R) 
    end
    if !isnothing(Rr)
        mpc.weights.Rr = Rr isa AbstractMatrix ? Rr : diagm(Rr) 
    end
    if !isnothing(rho)
        mpc.weights.rho = rho
    end
    if !isnothing(Qf)
        mpc.weights.Qf = Qf isa AbstractMatrix ? Qf : diagm(Qf) 
    end
end

function update_dynamics!(mpc::MPC,F,G)
    mpc.F,mpc.G = F,G
    mpc.mpQP = mpc2mpqp(mpc)
    return
end

