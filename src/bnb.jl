# Warm start of the branch and bound of DAQP for problems with binary controls (setting `bnb_warm_start`).
#
# The binary decision variables are those whose simple bounds have the binary sense. A binary decision variable
# is fixed at a value v for one solve by setting both of its bounds to v (DAQP treats equal bounds as an equality
# constraint until they differ again).

# Sets up the data of the warm start from mpc.mpQP: the binary decision variables and, for each of them, the
# decision variable that holds its control one step later. The stored solution is discarded.
function setup_bnb!(mpc::MPC)
    mpQP = mpc.mpQP
    n, nu = length(mpQP.f), mpc.model.nu
    ms = length(mpQP.bu)-size(mpQP.A,1) # Number of simple bounds
    # (With a prestabilizing feedback, the control bounds are general constraints and there are none)
    binary_ids = [j for j in 1:min(ms,n) if (mpQP.senses[j] & DAQP.BINARY) != 0]
    shift_ids = Int[]
    if !isempty(binary_ids)
        # The controls over the control horizon are T times the decision variables
        T = isempty(mpc.move_blocks) ? Matrix{Float64}(I,nu*mpc.Nc,n) : first(move_block_matrix(mpc))
        if size(T) == (nu*mpc.Nc,n)
            var_of_row = [findfirst(!iszero,view(T,r,:)) for r in 1:size(T,1)]
            for j in binary_ids
                r = findfirst(!iszero,view(T,:,j)) # Row of the first step of the decision variable
                ctrl,step = mod1(r,nu),cld(r,nu)
                push!(shift_ids,var_of_row[(min(step+1,mpc.Nc)-1)*nu+ctrl])
            end
        else
            empty!(binary_ids)
        end
    end
    mpc.bnb = BnBData(binary_ids,shift_ids,Float64[])
end

"""
    reset_bnb_warm_start!(mpc)

Discards the solution that is stored for the warm start of the branch and bound (setting `bnb_warm_start`),
so that the next solve does not form a candidate from it. This is appropriate when the problem changes
between two calls such that the previous solution is no longer representative, for example after a jump of
the state or the reference. [`setup!`](@ref) and the start of a `Simulation` discard the stored solution as well.
"""
function reset_bnb_warm_start!(mpc::MPC)
    empty!(mpc.bnb.xprev)
    return nothing
end

use_bnb_strategy(mpc::MPC) = mpc.settings.bnb_warm_start && !isempty(mpc.bnb.binary_ids)

# Data of one call of solve_bnb
mutable struct BnBSolve
    bu::Vector{Float64}             # Bounds for the parameter θ of the call
    bl::Vector{Float64}
    settings::DAQP.DAQPSettings     # Settings of the DAQP model at the start of the call
    time_limit::Float64             # Time limit of the whole call [s] (0 if there is none)
    t0::UInt64                      # Start of the call [ns]
    offset::Float64                 # Internal objective of DAQP minus the objective J, measured at the last solve
    iterations::Int
    nodes::Int
    qp_count::Int
end

function BnBSolve(mpc::MPC,θ)
    mpQP = mpc.mpQP
    mul!(mpQP._bth, mpQP.W, θ)
    mul!(mpQP._f, mpQP.f_theta, θ)
    mpQP._f .+= mpQP.f
    settings = DAQP.settings(mpc.opt_model)
    return BnBSolve(mpQP.bu .+ mpQP._bth, mpQP.bl .+ mpQP._bth, settings, settings.time_limit,
                    time_ns(), NaN, 0, 0, 0)
end

bnb_elapsed(s::BnBSolve) = (time_ns()-s.t0)/1e9

# One solve of DAQP for the parameter of `s`, with the decision variables `fix_ids` fixed at `fix_vals`.
# With `cutoff`, only solutions with an objective below `cutoff` are accepted (the objective of a solution
# that has been found for the same parameter, from which the offset of the internal objective of DAQP is
# known). A search (`search = true`) is limited to the time that remains of the time limit of the call; the
# other solves, in which all binary decision variables are fixed, are solved without a time limit, since
# they provide the integer-feasible fallback of the call.
function bnb_qp!(mpc::MPC, s::BnBSolve; fix_ids=Int[], fix_vals=Float64[], cutoff=nothing, search=false)
    mpQP,model = mpc.mpQP,mpc.opt_model
    changes = Dict{Symbol,Any}()
    if search && s.time_limit > 0
        remaining = s.time_limit-bnb_elapsed(s)
        remaining > 0 || return (x=fill(NaN,length(mpQP.f)), λ=zeros(length(mpQP.bu)), fval=NaN, flag=DAQP.TIMELIMIT)
        changes[:time_limit] = remaining
    elseif !search && s.time_limit > 0
        changes[:time_limit] = 0.0
    end
    if !isnothing(cutoff)
        # fval_bound is compared with the internal objective of DAQP, which exceeds J by an offset that depends
        # on the parameter (half of f'H⁻¹f without equality reduction)
        isnan(s.offset) && throw(ArgumentError("A cutoff requires a previous solve for the same parameter"))
        changes[:fval_bound] = cutoff+s.offset
    end
    mpQP._bu .= s.bu
    mpQP._bl .= s.bl
    for (j,v) in zip(fix_ids,fix_vals)
        mpQP._bu[j] = v
        mpQP._bl[j] = v
    end
    isempty(changes) || DAQP.settings(model,changes)
    DAQP.update(model,nothing,mpQP._f,nothing,mpQP._bu,mpQP._bl,nothing)
    x,fval,flag,info = DAQP.solve(model)
    isempty(changes) || DAQP.settings(model,s.settings)
    flag >= 1 && (s.offset = 0.5*unsafe_load(model.work).fval-fval)
    s.iterations += info.iterations
    s.nodes += info.nodes
    s.qp_count += 1
    return (x=x, λ=copy(info.λ), fval=fval, flag=Int(flag))
end

# The bound (l or u) that is nearest to v
nearest_bound(v,l,u) = abs(v-l) <= abs(u-v) ? l : u

# Candidate from the solution of the previous call: each binary decision variable takes the value of its
# control one step later in the previous solution (the last step is repeated), rounded to the nearest bound,
# and the continuous decision variables are optimized for these values. Returns nothing if there is no
# previous solution or if the candidate is infeasible.
function bnb_candidate(mpc::MPC, s::BnBSolve)
    bnb = mpc.bnb
    length(bnb.xprev) == length(mpc.mpQP.f) || return nothing
    vals = [nearest_bound(bnb.xprev[k],s.bl[j],s.bu[j]) for (j,k) in zip(bnb.binary_ids,bnb.shift_ids)]
    r = bnb_qp!(mpc,s; fix_ids=bnb.binary_ids, fix_vals=vals)
    return r.flag >= 1 ? r : nothing
end

# Solve of a problem with binary controls with the warm start of the branch and bound (see solve)
function solve_bnb(mpc::MPC,θ)
    s = BnBSolve(mpc,θ)
    cand = bnb_candidate(mpc,s)
    # The search only accepts solutions that are better than the candidate
    res = bnb_qp!(mpc,s; cutoff=isnothing(cand) ? nothing : cand.fval, search=true)
    source = :search
    if res.flag < 1 && !isnothing(cand)
        # No better solution than the candidate, or the time limit has been reached
        res,source = cand,:candidate
    end
    return bnb_result(mpc,s,res,source,cand)
end

function bnb_result(mpc::MPC, s::BnBSolve, res, source, cand)
    mpc.bnb.xprev = res.flag >= 1 ? copy(res.x) : Float64[]
    info = (x=res.x, λ=res.λ, fval=res.fval, exitflag=res.flag,
            status=get(DAQP.flag2status,res.flag,:Unknown),
            solve_time=bnb_elapsed(s), setup_time=0.0,
            iterations=s.iterations, nodes=s.nodes,
            source=source, candidate_fval=isnothing(cand) ? NaN : cand.fval, qp_count=s.qp_count)
    return copy(res.x),res.fval,Cint(res.flag),info
end
