struct Scenario
    x0::Vector{Float64}
    T::Float64
    N::Int
    r::VecOrMat{Float64}
    d::VecOrMat{Float64}
    l::VecOrMat{Float64}
    callback::Function
    dynamics::Union{Function,Nothing}
    get_measurement::Union{Function,Nothing}
end

struct Simulation
    ts::Vector{Float64}
    ys::Matrix{Float64}
    us::Matrix{Float64}
    xs::Matrix{Float64}
    rs::Matrix{Float64}
    ds::Matrix{Float64}
    xhats::Matrix{Float64}
    yms::Matrix{Float64}

    solve_times::Vector{Float64}

    mpc::Union{MPC,ExplicitMPC}
    scenario::Scenario
end

function Scenario(x0;T=-1.0,N=1000,r=nothing,d=nothing,l=nothing,
        callback = (x,u,d,k)->nothing, get_measurement=nothing, dynamics = nothing)
    r = isnothing(r) ? zeros(0,0) : Array{Float64}(r)
    d = isnothing(d) ? zeros(0,0) : Array{Float64}(d)
    l = isnothing(l) ? zeros(0,0) : Array{Float64}(l)
    return Scenario(Array{Float64}(x0),Float64(T),Int(N),r,d,l,callback,dynamics,get_measurement)
end

function Simulation(mpc::Union{MPC,ExplicitMPC}, scenario::Scenario)
    # Get number of time steps (prioritize T if it is entered)
    N = scenario.T < 0 ? scenario.N : Int(abs(ceil(scenario.T/mpc.model.Ts)))
    dynamics = isnothing(scenario.dynamics) ? mpc.model.true_dynamics : scenario.dynamics
    # Check if MPC has observer
    has_observer = !isnothing(mpc.state_observer)
    ny = has_observer ? size(mpc.state_observer.C,1) : mpc.model.ny
    if isnothing(scenario.get_measurement)
        get_measurement = if has_observer 
            (x,d) -> mpc.state_observer.C*x+mpc.state_observer.Dd*d+mpc.state_observer.h_offset
        else
            (x,d) -> mpc.model.C*x+mpc.model.Dd*d+mpc.model.h_offset
        end
    else
        get_measurement = scenario.get_measurement
    end

    x,u = scenario.x0,zeros(mpc.model.nu)
    
    xs = zeros(mpc.model.nx,N);
    ys = zeros(mpc.model.ny,N);
    rs = repeat(mpc.model.C*mpc.model.xo,1,N)
    ds = zeros(mpc.model.nd,N);
    ls = zeros(mpc.model.nu,N);
    us = zeros(mpc.model.nu,N)
    xhats = zeros(mpc.model.nx,N);
    # ys = yms if no observer
    yms = has_observer ? zeros(size(mpc.state_observer.C,1),N) : zeros(mpc.model.ny,N)

    solve_times = zeros(N)
    # Setup reference 
    if(!isempty(scenario.r))
        Nr = min(N,size(scenario.r,2))
        rs[:,1:Nr].= scenario.r[:,1:Nr]
        rs[:,size(scenario.r,2)+1:end] .= scenario.r[:,end] # hold last value
    end
    r_preview = mpc.settings.reference_preview && !isempty(scenario.r)

    # Setup disturbance
    if(!isempty(scenario.d))
        ds[:,1:size(scenario.d,2)].= scenario.d
        ds[:,size(scenario.d,2)+1:end] .= scenario.d[:,end] # hold last
    end

    # Setup linear cost trajectory
    if(!isempty(scenario.l))
        nl = size(scenario.l,2)
        Nl = min(N,nl)
        ls[:,1:Nl].= scenario.l[:,1:Nl]
        ls[:,nl+1:end] .= scenario.l[:,end] # hold last
    end
    l_preview = mpc.settings.linear_cost && !isempty(scenario.l)

    # Start the simulation
    has_observer && set_state!(mpc,scenario.x0)
    for k = 1:N
        xs[:,k], yms[:,k] = x, get_measurement(x,ds[:,k])
        ys[:,k] = has_observer ? mpc.model.C*x + mpc.model.Dd*ds[:,k] : yms[:,k]

        # Get state estimate
        xhat = has_observer ? correct_state!(mpc,yms[:,k],ds[:,k]) : x
        xhats[:,k] = xhat

        # Get linear cost and reference preview
        rk = r_preview ? get_preview(rs, k, mpc.Np) : rs[:,k]
        lk = l_preview ? get_preview(ls, k, mpc.Nc) : nothing

        solve_times[k] = @elapsed u = compute_control(mpc,xhat;r=rk, d=ds[:,k],l=lk)

        has_observer && predict_state!(mpc,u,ds[:,k])

        x = dynamics(x,u,ds[:,k])
        scenario.callback(x,u,ds[:,k],k)
        us[:,k] = u
    end
    Ts = mpc.model.Ts < 0.0 ? 1.0 : mpc.model.Ts
    return Simulation(collect(Ts*(0:1:N-1)),ys,us,xs,rs,ds,xhats,yms,solve_times,mpc,scenario)
end

function Simulation(dynamics, mpc::Union{MPC,ExplicitMPC};x0=zeros(mpc.model.nx),T=-1.0, N=1000, r=nothing,d=nothing, l=nothing, callback=(x,u,d,k)->nothing, get_measurement= nothing)
    r = isnothing(r) ? zeros(0,0) : Array{Float64}(r)
    d = isnothing(d) ? zeros(0,0) : Array{Float64}(d)
    l = isnothing(l) ? zeros(0,0) : Array{Float64}(l)
    return Simulation(mpc,Scenario(Array{Float64}(x0),Float64(T),Int(N),r,d,l,
                                   callback,dynamics,get_measurement))
end

Simulation(mpc::Union{MPC,ExplicitMPC}; kwargs...) = Simulation(mpc.model.true_dynamics, mpc; kwargs...)

function get_preview(rs,k,Nc)
    preview = zeros(size(rs,1), Nc)
    @views for i in 1:Nc
        preview[:, i] .= rs[:, min(k + i, end)]
    end
    return preview
end
"""
    get_reference_preview(rs, k, Np)

Extract reference preview from reference trajectory starting at time step k.
"""
get_reference_preview(rs, k, Np) = get_preview(rs, k , Np)

"""
    get_linear_cost_preview(ls, k, Nc)

Extract linear cost preview from linear cost trajectory starting at time step k.
"""
get_linear_cost_preview(ls, k, Nc) = get_prewview(ls,k-1,Nc)


using RecipesBase

@recipe function f(sim::Simulation; yids=1:sim.mpc.model.ny,uids=1:sim.mpc.model.nu,xids=[])


    layout = Tuple{Int,Int}[]
    length(yids) == 0 || push!(layout,(length(yids), 1))
    length(uids) == 0 || push!(layout,(length(uids), 1))
    length(xids) == 0 || push!(layout,(length(xids), 1))
    layout := reshape(layout,1,length(layout))

    if sim.mpc isa MPC
        umin,umax = sim.mpc.umin, sim.mpc.umax
    else
        umin,umax = sim.mpc.bst.clipping[:,1],sim.mpc.bst.clipping[:,2]
    end
    # Plot y
    id = 1 
    for i in yids 
        @series begin
            linecolor := :black 
            subplot   --> id
            linestyle := :dash
            linewidth --> 0.4
            label     --> "reference"
            primary   --> false
            sim.ts, sim.rs[i,:]
        end
        @series begin
            yguide  --> latexify(make_subscript(string(sim.mpc.model.labels.y[i])))
            color   --> 1
            subplot --> id 
            #label--> latexify(make_subscript(string(sim.mpc.model.labels.y[i])))
            legend --> false
            if i == length(yids) 
                xguide --> (sim.mpc.model.Ts < 0 ? "Time step" : "Time [s]")
            end
            sim.ts, sim.ys[i, :]
        end

        id+=1
    end

    # Plot u 
    for i in uids 
        # lower bound
        if(length(umin) >= i && umin[i] > -1e12)
            @series begin
                color     := :black 
                subplot   --> id
                linestyle --> :dash
                linewidth --> 0.8
                primary   --> false
                sim.ts[[1,end]], fill(umin[i],2) 
            end
        end

        # upper bound
        if(length(umax) >= i && umax[i] < 1e12)
            @series begin
                color     := :black 
                subplot   --> id
                linestyle --> :dash
                linewidth --> 0.8 
                primary   --> false
                sim.ts[[1,end]], fill(umax[i],2) 
            end
        end
        @series begin
            yguide--> latexify(make_subscript(string(sim.mpc.model.labels.u[i])))
            color      --> 1
            subplot    --> id
            seriestype --> :steppost
            legend --> false
            #label--> latexify(make_subscript(string(sim.mpc.model.labels.u[i])))
            if i == length(uids) 
                xguide --> (sim.mpc.model.Ts < 0 ? "Time step" : "Time [s]")
            end
            sim.ts, sim.us[i, :]
        end
        id+=1
    end

    ## Plot x 
    for i in xids 
        @series begin
            yguide  --> latexify(make_subscript(string(sim.mpc.model.labels.x[i])))
            color  --> 1
            subplot--> id
            legend --> false
            if i == length(xids)
                xguide --> (sim.mpc.model.Ts < 0 ? "Time step" : "Time [s]")
            end
            sim.ts, sim.xs[i, :]
        end
        id+=1
    end
end

"""
    evaluate_cost(mpc,sim;Q,Rr,S)
Compute the cost 0.5 ∑ x'*Q x + u' R u + Δu' Rr Δu + x' S u
"""
function evaluate_cost(mpc::MPC,sim::Simulation;
        Q=mpc.weights.Q, R = mpc.weights.R, Rr = mpc.weights.Rr, S = mpc.weights.S)
    return evaluate_cost(mpc,sim.xs,sim.us,sim.rs;Q,R,Rr,S)
end
