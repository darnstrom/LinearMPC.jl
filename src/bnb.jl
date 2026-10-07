# Solution strategies of the branch and bound of DAQP for problems with binary controls: a warm start from the
# solution of the previous call (setting `bnb_warm_start`), and groups of binary controls that are resolved after
# the others (defer_binary_controls!).
#
# The binary decision variables are those whose simple bounds have the binary sense. A binary decision variable
# is fixed at a value v for one solve by setting both of its bounds to v (DAQP treats equal bounds as an equality
# constraint until they differ again). The deferred binary decision variables are relaxed by an update of the
# senses of the DAQP model in which their binary sense is cleared, and restored by an update with the original
# senses. This requires that an update of the senses updates the binary constraints of the branch and bound of
# DAQP (darnstrom/daqp#208), see daqp_updates_binaries.

# Sets up the data of the solution strategies from mpc.mpQP: the binary decision variables, for each of them
# the decision variable that holds its control one step later, and the binary decision variables of the
# groups of deferred binary controls. The stored solution is discarded.
function setup_bnb!(mpc::MPC)
    mpQP = mpc.mpQP
    n, nu = length(mpQP.f), mpc.model.nu
    ms = length(mpQP.bu)-size(mpQP.A,1) # Number of simple bounds
    # (With a prestabilizing feedback, the control bounds are general constraints and there are none)
    binary_ids = [j for j in 1:min(ms,n) if (mpQP.senses[j] & DAQP.BINARY) != 0]
    shift_ids,ctrls,steps,durations = Int[],Int[],Int[],Float64[]
    if !isempty(binary_ids)
        # The controls over the control horizon are T times the decision variables
        T = isempty(mpc.move_blocks) ? Matrix{Float64}(I,nu*mpc.Nc,n) : first(move_block_matrix(mpc))
        if size(T) == (nu*mpc.Nc,n)
            var_of_row = [findfirst(!iszero,view(T,r,:)) for r in 1:size(T,1)]
            for j in binary_ids
                r = findfirst(!iszero,view(T,:,j)) # Row of the first step of the decision variable
                ctrl,step = mod1(r,nu),cld(r,nu)
                push!(shift_ids,var_of_row[(min(step+1,mpc.Nc)-1)*nu+ctrl])
                push!(ctrls,ctrl)
                push!(steps,step)
                push!(durations,count(!iszero,view(T,:,j)))
            end
        else
            empty!(binary_ids)
        end
    end
    deferred = mpc.bnb.deferred
    groups = BnBGroup[]
    if !isempty(binary_ids)
        for d in deferred
            g = setup_deferred_group(d,ctrls,steps,durations)
            isnothing(g) || push!(groups,g)
        end
    end
    is_deferred = fill(false,length(binary_ids))
    for g in groups
        is_deferred[g.pos] .= true
    end
    groups = [g.kind == :rounding ? BnBGroup(g.kind,g.pos,g.weights,g.durations,
                                             bnb_logic_rows(mpQP,binary_ids,is_deferred,g.pos)) : g
              for g in groups]
    relaxed_senses = copy(mpQP.senses)
    for j in binary_ids[is_deferred]
        relaxed_senses[j] &= ~Cint(DAQP.BINARY)
    end
    mpc.bnb = BnBData(binary_ids,shift_ids,Float64[],deferred,groups,is_deferred,relaxed_senses)
end

# The binary decision variables of the group d, given the control (ctrls), the first step (steps) and the
# number of time steps (durations) of each binary decision variable. Returns nothing if d has none.
function setup_deferred_group(d::DeferredBinaryGroup, ctrls, steps, durations)
    poss = [sort(findall(==(c),ctrls); by=p->steps[p]) for c in d.ids]
    missing_ids = d.ids[isempty.(poss)]
    isempty(missing_ids) || @warn "The deferred controls $missing_ids are not binary"
    all(isempty,poss) && return nothing
    if d.resolution == :enumerate
        return BnBGroup(:enumeration,reduce(vcat,poss),Float64[],Float64[],BnBLogicRow[])
    elseif !isempty(d.weights)
        if all(p -> steps[p] == steps[first(poss)], poss)
            # Step by step, each step in the order of the controls
            return BnBGroup(:encoding,vec(permutedims(reduce(hcat,poss))),d.weights,Float64[],BnBLogicRow[])
        end
        @warn "The deferred controls $(d.ids) do not have the same binary steps, all combinations of their binary decision variables are evaluated instead of the integers that they encode"
        return BnBGroup(:enumeration,reduce(vcat,poss),Float64[],Float64[],BnBLogicRow[])
    end
    pos = first(poss) # A single control (see defer_binary_controls!)
    return BnBGroup(:rounding,pos,Float64[],Float64.(durations[pos]),BnBLogicRow[])
end

# The logic constraints of the binary decision variables at the positions `pos`: the general constraints that
# are neither soft nor binary, whose bounds do not depend on the parameter, and that involve only these and the
# binary decision variables that are not deferred, at least one of the former. With all of these binary
# decision variables fixed, they determine whether a candidate is infeasible without a solve.
function bnb_logic_rows(mpQP, binary_ids, is_deferred, pos)
    ms = length(mpQP.bu)-size(mpQP.A,1)
    pos_of = Dict(j => p for (p,j) in enumerate(binary_ids) if !is_deferred[p] || p in pos)
    rows = BnBLogicRow[]
    for r in 1:size(mpQP.A,1)
        i = ms+r
        (mpQP.senses[i] & (DAQP.SOFT|DAQP.BINARY)) == 0 || continue
        iszero(view(mpQP.W,i,:)) || continue
        nz = findall(!iszero,view(mpQP.A,r,:))
        all(j -> haskey(pos_of,j), nz) && any(j -> pos_of[j] in pos, nz) || continue
        push!(rows,BnBLogicRow([pos_of[j] for j in nz],mpQP.A[r,nz],mpQP.bl[i],mpQP.bu[i]))
    end
    return rows
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

"""
    defer_binary_controls!(mpc, ids; weights=nothing, resolution=:auto)

Adds a group of deferred binary controls `ids`. The branch and bound first runs with the binary controls of
all groups relaxed, and these are resolved afterwards with the other binary controls fixed at their values in
the relaxed solution (see [`solve`](@ref LinearMPC.solve)). This reduces the number of nodes when the deferred
binary controls are numerous but have a small effect on the decisions of the other binary controls, for
example bits that encode an integer setpoint.

A group is resolved at the binary steps of its controls (the steps of the binary horizon, or the move blocks
within it), from its normalized relaxed values `b ∈ [0,1]` (0 at the lower and 1 at the upper bound):
* With `weights`, the controls `ids` encode the integer `∑ weights[i]*b[i]` at each binary step, for example
  `weights = [1,2,4]` for a setpoint `0,…,7` in three bits. The candidates are the encodable integers just
  below and above the relaxed value at each binary step, combined over the steps. An integer that several
  combinations of the bits encode is represented by the first of them in the order of the binary numbers
  `0:2^length(ids)-1`. If the controls affect the problem only through the encoded integers (for example
  a setpoint), the optimal objective for fixed other binary controls is a jointly convex function of the
  encoded integers of all steps and groups (the minimum of a convex QP over the continuous decision
  variables, with the encoded integers as further variables). In one dimension, its minimum over the
  integers is at one of the two integers next to the minimizer; in several dimensions, the best
  combination of these neighbours is not necessarily the optimal one.
* A single control without `weights` is resolved by sum-up rounding (Sager, Bock and Diehl, Math. Program.
  133:1–23, 2012) of its relaxed values `a_k` over its binary steps: `b_k = 1` if
  `∑_{i≤k} a_i Δ_i - ∑_{i<k} b_i Δ_i ≥ Δ_k/2`, where `Δ_k` is the number of time steps of the move block of
  step `k` (1 without move blocks). The candidates are this sequence, the sequence in which the first step
  that is off is on, and the sequence in which the last step that is on is off. A step is only changed if
  the logic constraints of the control allow it: the general constraints that involve only its binary steps
  and binary controls that are not deferred, with bounds that do not depend on the parameter (for example
  `δ_g ≤ δ` from [`add_logic_constraint!`](@ref), with δ not deferred). Other infeasible candidates are
  discarded after their QP. At a single binary step, the candidates are the values 0 and 1.
* `resolution = :enumerate` evaluates all combinations of the binary decision variables of the group.

Each combination of the candidates of the groups is a QP with all binary controls fixed; infeasible ones are
discarded. If there are more than the setting `deferred_max_combinations`, the groups are resolved in turn,
in the order of their declaration: the candidates of a group are evaluated with the groups before it at
their best and the groups after it at their first candidate. A group with more candidates than
`deferred_max_combinations` has only its first candidate, the rounded relaxed values.

The best combination is accepted if its objective exceeds that of the relaxed search by at most the setting
`deferred_tol`. Otherwise, the branch and bound runs without relaxation, with the best combination as cutoff.
Since this check compares with a lower bound of the optimal objective, the returned solution is within
`deferred_tol` of the optimum (in addition to the tolerances `abs_subopt` and `rel_subopt` of DAQP), also where
the candidates above do not contain the optimal combination, for example where the rounding of the groups
separately is not optimal in several dimensions, or where the other binary controls of the relaxed solution
are not optimal.

The regularization that is added to the objective for each binary control at each binary step,
`(u-umin)(u-umax)/2` up to a constant, is equal at both bounds and lower between them. The relaxed objective
therefore remains a lower bound, but it is lower than without the regularization, by at most
`(umax-umin)^2/8` for each binary step of a relaxed control, and the comparison with `deferred_tol` includes
this difference. Its quadratic term is part of the Hessian, which only an update of the Hessian (a new
factorization) could change for the relaxed search. A change of its linear term alone cannot reduce the
difference: with the same quadratic term, any other linear term exceeds the objective at one of the bounds,
and the relaxed objective would no longer be a lower bound.

The relaxation requires a DAQP library in which an update of the senses updates the binary constraints of the
branch and bound (darnstrom/daqp#208); [`solve`](@ref LinearMPC.solve) throws an error otherwise. The C code
that [`codegen`](@ref LinearMPC.codegen) generates does not depend on it.

The controls `ids` must be binary controls (see [`set_binary_controls!`](@ref)) and belong to no other group.
Several controls require `weights` or `resolution = :enumerate`. [`clear_deferred_binary_controls!`](@ref)
removes all groups.
"""
function defer_binary_controls!(mpc::MPC, ids; weights=nothing, resolution=:auto)
    ids = Int.(collect(ids))
    isempty(ids) && throw(ArgumentError("ids must contain at least one control"))
    for id in ids
        _validate_input_id(id, mpc.model.nu, "ids")
    end
    allunique(ids) || throw(ArgumentError("ids must not contain a control more than once"))
    for d in mpc.bnb.deferred
        common = intersect(d.ids,ids)
        isempty(common) || throw(ArgumentError("The controls $common belong to a group already"))
    end
    resolution in (:auto,:enumerate) || throw(ArgumentError("resolution must be :auto or :enumerate"))
    weights = isnothing(weights) ? Float64[] : Float64.(collect(weights))
    if resolution == :enumerate
        isempty(weights) || throw(ArgumentError("weights are not used with resolution = :enumerate"))
    elseif isempty(weights)
        length(ids) == 1 || throw(ArgumentError("Several controls require weights or resolution = :enumerate"))
    else
        length(weights) == length(ids) ||
            throw(ArgumentError("weights must have one entry per control in ids, got $(length(weights)) for $(length(ids)) controls"))
        length(weights) <= 20 || throw(ArgumentError("At most 20 controls can encode an integer"))
    end
    push!(mpc.bnb.deferred, DeferredBinaryGroup(ids,weights,resolution))
    mpc.mpqp_issetup = false
    return nothing
end

"""
    clear_deferred_binary_controls!(mpc)

Removes all groups of deferred binary controls (see [`defer_binary_controls!`](@ref)).
"""
function clear_deferred_binary_controls!(mpc::MPC)
    empty!(mpc.bnb.deferred)
    mpc.mpqp_issetup = false
    return nothing
end

use_bnb_strategy(mpc::MPC) = !isempty(mpc.bnb.binary_ids) &&
    (mpc.settings.bnb_warm_start || !isempty(mpc.bnb.groups))

# Whether an update of the senses of a DAQP model updates the binary constraints of its branch and bound
# (darnstrom/daqp#208), which the relaxation of deferred binary controls requires. Without it, DAQP keeps the
# binary constraints of its setup (unless the equality reduction is active), so that a relaxed search would be
# the exact one. This is determined once, from the problem min x'x/2 - x₁ - x₂ s.t. x₁ + x₂ ≤ 0.8 with
# x₁ ∈ {0,1}: its solution has x₁ = 0, whereas x₁ = 0.4 with x₁ relaxed.
const DAQP_UPDATES_BINARIES = Ref{Union{Nothing,Bool}}(nothing)
function daqp_updates_binaries()
    if isnothing(DAQP_UPDATES_BINARIES[])
        model = DAQP.Model()
        senses = Cint[DAQP.BINARY,0,0]
        flag,_ = DAQP.setup(model,Matrix{Float64}(I,2,2),[-1.0,-1.0],[1.0 1.0],[1.0,10.0,0.8],[0.0,-10.0,-1e30],senses)
        supported = flag >= 0
        for (s,x1) in ((Cint[0,0,0],0.4),(senses,0.0)) # Relaxed and restored
            supported || break
            DAQP.update(model,nothing,nothing,nothing,nothing,nothing,s)
            x,_,flag,_ = DAQP.solve(model)
            supported = flag >= 1 && abs(x[1]-x1) < 1e-6
        end
        DAQP_UPDATES_BINARIES[] = supported
    end
    return DAQP_UPDATES_BINARIES[]
end

# Data of one call of solve_bnb
mutable struct BnBSolve
    bu::Vector{Float64}             # Bounds for the parameter θ of the call
    bl::Vector{Float64}
    settings::DAQP.DAQPSettings     # Settings of the DAQP model at the start of the call
    time_limit::Float64             # Time limit of the whole call [s] (0 if there is none)
    t0::UInt64                      # Start of the call [ns]
    offset::Float64                 # Internal objective of DAQP minus the objective J, measured at the last solve
    relaxed::Bool                   # Whether the deferred binary decision variables are relaxed in the DAQP model
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
                    time_ns(), NaN, false, 0, 0, 0)
end

bnb_elapsed(s::BnBSolve) = (time_ns()-s.t0)/1e9
bnb_at_time_limit(s::BnBSolve) = s.time_limit > 0 && bnb_elapsed(s) >= s.time_limit

# One solve of DAQP for the parameter of `s`, with the decision variables `fix_ids` fixed at `fix_vals`, and
# with the deferred binary decision variables relaxed if `relaxed = true` (`relaxed = nothing` keeps the senses
# of the previous solve, for solves in which all binary decision variables are fixed). With `cutoff`, only
# solutions with an objective below `cutoff` are accepted (the objective of a solution that has been found for
# the same parameter, from which the offset of the internal objective of DAQP is known). A search
# (`search = true`) is limited to the time that remains of the time limit of the call; the other solves, in
# which all binary decision variables are fixed or relaxed, are solved without a time limit, since they provide
# the integer-feasible fallback of the call.
function bnb_qp!(mpc::MPC, s::BnBSolve; fix_ids=Int[], fix_vals=Float64[], relaxed=false, cutoff=nothing, search=false)
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
    senses = nothing
    if !isnothing(relaxed) && relaxed != s.relaxed
        senses = relaxed ? mpc.bnb.relaxed_senses : mpQP.senses
        s.relaxed = relaxed
    end
    isempty(changes) || DAQP.settings(model,changes)
    DAQP.update(model,nothing,mpQP._f,nothing,mpQP._bu,mpQP._bl,senses)
    x,fval,flag,info = DAQP.solve(model)
    isempty(changes) || DAQP.settings(model,s.settings)
    # The internal objective is info.fval_ldp if DAQPBase provides it, and otherwise read from the workspace
    flag >= 1 && (s.offset = (haskey(info,:fval_ldp) ? info.fval_ldp : 0.5*unsafe_load(model.work).fval)-fval)
    s.iterations += info.iterations
    s.nodes += info.nodes
    s.qp_count += 1
    return (x=x, λ=copy(info.λ), fval=fval, flag=Int(flag))
end

# Restores the senses of the DAQP model if the deferred binary decision variables are relaxed
function bnb_restore!(mpc::MPC, s::BnBSolve)
    s.relaxed || return nothing
    DAQP.update(mpc.opt_model,nothing,nothing,nothing,nothing,nothing,mpc.mpQP.senses)
    s.relaxed = false
    return nothing
end

# The bound (l or u) that is nearest to v
nearest_bound(v,l,u) = abs(v-l) <= abs(u-v) ? l : u

# Tolerance of the rounding of the normalized relaxed values of deferred binary decision variables: a value
# within it of 1/2 is rounded down, and the sum-up rounding switches on within it of its threshold. This makes
# the rounding of values at 1/2, toward which the regularization of the binary controls draws relaxed binary
# decision variables, independent of rounding errors (also in the generated C code).
const BNB_ROUND_TOL = 1e-6

# Candidate from the solution of the previous call: each binary decision variable takes the value of its
# control one step later in the previous solution (the last step is repeated), rounded to the nearest bound,
# and the continuous decision variables are optimized for these values. The deferred binary decision
# variables are instead relaxed and then resolved for the other ones (see bnb_resolve_deferred). Returns
# nothing if the warm start is not enabled, if there is no previous solution or if the candidate is infeasible.
function bnb_candidate(mpc::MPC, s::BnBSolve)
    bnb = mpc.bnb
    mpc.settings.bnb_warm_start || return nothing
    length(bnb.xprev) == length(mpc.mpQP.f) || return nothing
    vals = [nearest_bound(bnb.xprev[k],s.bl[j],s.bu[j]) for (j,k) in zip(bnb.binary_ids,bnb.shift_ids)]
    if isempty(bnb.groups)
        r = bnb_qp!(mpc,s; fix_ids=bnb.binary_ids, fix_vals=vals)
        return r.flag >= 1 ? r : nothing
    end
    O = .!bnb.is_deferred
    r = bnb_qp!(mpc,s; fix_ids=bnb.binary_ids[O], fix_vals=vals[O], relaxed=true)
    return r.flag >= 1 ? bnb_resolve_deferred(mpc,s,r) : nothing
end

# The value Σ w[i]*b[i] that the bits b of the integer c encode
encoded_value(w,c) = sum(w[i] for i in eachindex(w) if (c >> (i-1)) & 1 == 1; init=0.0)

# Candidates of a group (vectors of the values of its binary decision variables, true for the upper bound),
# given the normalized relaxed values z and the values vals of all binary decision variables, in which those
# that are not deferred are fixed. The first candidate rounds the relaxed values.
function bnb_group_candidates(g::BnBGroup, z, vals, s::BnBSolve, l, u, max_combinations)
    if g.kind == :encoding
        # At each step, the codes of the encodable integers just below and above the relaxed value, the
        # nearest first
        w = g.weights
        m = length(w)
        lo,hi = sum(min.(w,0.0)),sum(max.(w,0.0))
        tol = 1e-6*(1+maximum(abs,w))
        near,alt = Int[],Int[]
        for k in 1:length(g.pos)÷m
            p = g.pos[(k-1)*m+1:k*m]
            v = clamp(sum(w[i]*z[p[i]] for i in 1:m),lo,hi)
            below,above = (-Inf,-1),(Inf,-1)
            for c in 0:2^m-1
                val = encoded_value(w,c)
                val <= v+tol && val > below[1] && (below = (val,c))
                val >= v-tol && val < above[1] && (above = (val,c))
            end
            nearest,other = v-below[1] <= above[1]-v+tol ? (below,above) : (above,below)
            push!(near,nearest[2])
            push!(alt,below[1] == above[1] ? -1 : other[2])
        end
        two = findall(>=(0),alt) # The steps with two candidates
        ncand = length(two) < 30 && 2^length(two) <= max_combinations ? 2^length(two) : 1
        return map(0:ncand-1) do c
            codes = copy(near)
            for (t,k) in enumerate(two)
                (c >> (t-1)) & 1 == 1 && (codes[k] = alt[k])
            end
            [(codes[(q-1)÷m+1] >> ((q-1)%m)) & 1 == 1 for q in eachindex(g.pos)]
        end
    elseif g.kind == :rounding
        # Sum-up rounding, and the sequences with the first step that is off on and the last step that is on
        # off, if the logic constraints allow these changes
        a,Δ = z[g.pos],g.durations
        b,d = fill(false,length(a)),0.0
        for k in eachindex(a)
            d += a[k]*Δ[k]
            b[k] = d >= (0.5-BNB_ROUND_TOL)*Δ[k]
            b[k] && (d -= Δ[k])
        end
        cands = [b]
        v = copy(vals)
        v[g.pos] .= b
        for (state,ks) in ((false,eachindex(b)),(true,reverse(eachindex(b))))
            for k in ks
                b[k] == state || continue
                v[g.pos[k]] = !state
                feasible = bnb_logic_feasible(g,v,l,u,s.settings.primal_tol)
                v[g.pos[k]] = state
                if feasible
                    push!(cands,setindex!(copy(b),!state,k))
                    break
                end
            end
        end
        return length(cands) <= max_combinations ? cands : cands[1:1]
    else # :enumeration
        n = length(g.pos)
        near = z[g.pos] .> 0.5+BNB_ROUND_TOL
        ncand = n < 30 && 2^n <= max_combinations ? 2^n : 1
        return [near .⊻ digits(Bool,c;base=2,pad=n) for c in 0:ncand-1]
    end
end

# Whether the values `vals` of the binary decision variables (true for the upper bound u, otherwise the lower
# bound l) satisfy the logic constraints of the group g
function bnb_logic_feasible(g::BnBGroup, vals, l, u, tol)
    for row in g.rows
        act = sum(c*(vals[p] ? u[p] : l[p]) for (p,c) in zip(row.pos,row.coef))
        (act > row.upper+tol || act < row.lower-tol) && return false
    end
    return true
end

# Resolves the deferred binary decision variables of the solution `r`, in which the other binary decision
# variables are integer feasible: these are fixed at their values in `r`, and the combinations of the
# candidates of the groups are evaluated (see defer_binary_controls!). Returns the best of the resulting
# solutions (`r` if its deferred binary decision variables are integer feasible), or nothing if all are
# infeasible.
function bnb_resolve_deferred(mpc::MPC, s::BnBSolve, r)
    bnb = mpc.bnb
    B,D = bnb.binary_ids,bnb.is_deferred
    l,u = s.bl[B],s.bu[B]
    xb = r.x[B]
    all(min.(xb[D].-l[D],u[D].-xb[D]) .<= s.settings.primal_tol) && return r
    max_combinations = mpc.settings.deferred_max_combinations
    max_combinations >= 1 || throw(ArgumentError("The setting deferred_max_combinations must be positive"))
    z = clamp.((xb.-l)./(u.-l),0,1)
    vals = z .> 0.5+BNB_ROUND_TOL
    cands = [bnb_group_candidates(g,z,vals,s,l,u,max_combinations) for g in bnb.groups]
    best = nothing
    function evaluate(choice)
        for (k,c) in enumerate(choice)
            vals[bnb.groups[k].pos] .= cands[k][c]
        end
        rk = bnb_qp!(mpc,s; fix_ids=B, fix_vals=ifelse.(vals,u,l), relaxed=nothing)
        rk.flag >= 1 || return Inf
        (isnothing(best) || rk.fval < best.fval) && (best = rk)
        return rk.fval
    end
    nc = length.(cands)
    total = 1
    for n in nc
        total = min(total*n,max_combinations+1)
    end
    if total <= max_combinations
        # All combinations, with the first group changing fastest
        for choice in Iterators.product((1:n for n in nc)...)
            evaluate(choice)
        end
    else
        # The groups in turn: the candidates of a group with the groups before it at their best and the
        # groups after it at their first candidate
        choice = ones(Int,length(nc))
        current = evaluate(choice)
        for g in eachindex(nc)
            best_c = 1
            for c in 2:nc[g]
                choice[g] = c
                f = evaluate(choice)
                f < current && ((current,best_c) = (f,c))
            end
            choice[g] = best_c
        end
    end
    return best
end

# Solve of a problem with binary controls with the solution strategies of the branch and bound (see solve)
function solve_bnb(mpc::MPC,θ)
    if !isempty(mpc.bnb.groups) && !daqp_updates_binaries()
        error("Deferred binary controls require a DAQP library in which an update of the senses updates the binary constraints of the branch and bound (darnstrom/daqp#208). With the loaded library, the deferred binary controls would not be relaxed.")
    end
    s = BnBSolve(mpc,θ)
    try
        cand = bnb_candidate(mpc,s)
        isempty(mpc.bnb.groups) || return solve_deferred(mpc,s,cand)
        # The search only accepts solutions that are better than the candidate
        res = bnb_qp!(mpc,s; cutoff=isnothing(cand) ? nothing : cand.fval, search=true)
        source = :search
        if res.flag < 1 && !isnothing(cand)
            # No better solution than the candidate, or the time limit has been reached
            res,source = cand,:candidate
        end
        return bnb_result(mpc,s,res,source,cand)
    finally
        bnb_restore!(mpc,s)
    end
end

# Solve with deferred binary controls (and the candidate `cand` of the warm start, if any)
function solve_deferred(mpc::MPC, s::BnBSolve, cand)
    tol = mpc.settings.deferred_tol
    tol >= 0 || throw(ArgumentError("The setting deferred_tol must be nonnegative"))
    # 1. Search with the deferred binary decision variables relaxed. Its objective is a lower bound (within the
    #    suboptimality tolerances of DAQP) if it completes; with the candidate as cutoff, it only accepts
    #    solutions that are better than the candidate.
    a = bnb_qp!(mpc,s; relaxed=true, cutoff=isnothing(cand) ? nothing : cand.fval, search=true)
    if a.flag < 1
        # No solution of the relaxed problem is better than the candidate, or the time limit has been reached
        isnothing(cand) || return bnb_result(mpc,s,cand,:candidate,cand)
        a.flag == DAQP.TIMELIMIT && return bnb_result(mpc,s,a,:search,cand)
        full = bnb_qp!(mpc,s; search=true)
        return bnb_result(mpc,s,full,:search,cand)
    end
    at_time_limit = bnb_at_time_limit(s)
    # 2. The deferred binary decision variables resolved for the other binary decision variables of the
    #    relaxed solution
    resolved = bnb_resolve_deferred(mpc,s,a)
    best,source = resolved,:deferred
    if !isnothing(cand) && (isnothing(best) || cand.fval <= best.fval)
        best,source = cand,:candidate
    end
    if at_time_limit
        # The relaxed search has reached the time limit, so its objective is not a lower bound
        isnothing(best) && return bnb_result(mpc,s,merge(a,(flag=DAQP.TIMELIMIT,)),:search,cand,a.fval)
        return bnb_result(mpc,s,best,source,cand,a.fval)
    end
    # 3. The best solution is accepted if its objective exceeds that of the relaxed search by at most tol,
    #    otherwise the full search runs with it as cutoff
    !isnothing(best) && best.fval-a.fval <= tol && return bnb_result(mpc,s,best,source,cand,a.fval)
    full = bnb_qp!(mpc,s; cutoff=isnothing(best) ? nothing : best.fval, search=true)
    (full.flag >= 1 || isnothing(best)) && return bnb_result(mpc,s,full,:search,cand,a.fval)
    return bnb_result(mpc,s,best,source,cand,a.fval)
end

function bnb_result(mpc::MPC, s::BnBSolve, res, source, cand, relaxed_fval=NaN)
    mpc.bnb.xprev = res.flag >= 1 ? copy(res.x) : Float64[]
    info = (x=res.x, λ=res.λ, fval=res.fval, exitflag=res.flag,
            status=get(DAQP.flag2status,res.flag,:Unknown),
            solve_time=bnb_elapsed(s), setup_time=0.0,
            iterations=s.iterations, nodes=s.nodes,
            source=source, candidate_fval=isnothing(cand) ? NaN : cand.fval,
            relaxed_fval=relaxed_fval, qp_count=s.qp_count)
    return copy(res.x),res.fval,Cint(res.flag),info
end
