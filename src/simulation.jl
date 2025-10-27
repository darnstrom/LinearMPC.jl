struct Simulation
    ts::Vector{Float64}
    ys::Matrix{Float64}
    us::Matrix{Float64}
    xs::Matrix{Float64}
    rs::Matrix{Float64}
    ds::Matrix{Float64}
    xhats::Matrix{Float64}

    solve_times ::Vector{Float64}

    mpc::Union{MPC,ExplicitMPC}
end

function Simulation(dynamics, mpc::Union{MPC,ExplicitMPC}; x0=zeros(mpc.model.nx),N=1000, r=nothing,d=nothing, callback=(x,u,d,k)->nothing, get_measurement= nothing)

    # Check if MPC has observer
    has_observer = !isnothing(mpc.state_observer)
    ny = has_observer ? size(mpc.state_observer.C,1) : mpc.model.ny
    if isnothing(get_measurement)
        if has_observer 
            get_measurement = (x,d) -> mpc.state_observer.C*x+mpc.model.Dd*d
        else
            get_measurement = (x,d) -> mpc.model.C*x+mpc.model.Dd*d
        end
    end

    x,u = x0,zeros(mpc.model.nu)
    
    xs = zeros(mpc.model.nx,N);
    ys = zeros(ny,N);
    rs = zeros(mpc.model.ny,N);
    ds = zeros(mpc.model.nd,N);
    us = zeros(mpc.model.nu,N)
    xhats = zeros(mpc.model.nx,N);

    solve_times = zeros(N)
    # Setup reference 
    if(!isnothing(r))
        Nr = min(N,size(r,2))
        rs[:,1:Nr].= r[:,1:Nr]
        rs[:,size(r,2)+1:end] .= r[:,end] # hold last value of reference
    end

    # Setup disturbance
    if(!isnothing(d))
        ds[:,1:size(d,2)].= d
        ds[:,size(d,2)+1:end] .= d[:,end] # hold last
    end

    # Start the simulation
    has_observer && set_state!(mpc.state_observer,x0)
    for k = 1:N
        xs[:,k], ys[:,k] = x, get_measurement(x,ds[:,k])

        # Get state estimate
        xhat = has_observer ? correct!(mpc.state_observer,ys[:,k]) : x
        xhats[:,k] = xhat
        # Get reference
        if mpc.settings.reference_preview && !isnothing(r)
            # Reference preview mode: provide future references
            r = get_reference_preview(rs, k, mpc.Np)
        else # Standard mode: provide current reference
            r = rs[:,k]
        end

        solve_times[k] = @elapsed u = compute_control(mpc,xhat;r,d=ds[:,k])

        has_observer && predict!(mpc.state_observer,u)
        
        x = dynamics(x,u,ds[:,k])
        callback(x,u,ds[:,k],k)
        us[:,k] = u
    end
    Ts = mpc.model.Ts < 0.0 ? 1.0 : mpc.model.Ts
    return Simulation(collect(Ts*(0:1:N-1)),ys,us,xs,rs,ds,xhats,solve_times,mpc)
end

"""
    get_reference_preview(rs, k, Np)

Extract reference preview from reference trajectory starting at time step k.
"""
function get_reference_preview(rs, k, Np)
    ny, N = size(rs)
    r_preview = zeros(ny, Np)
    
    @views for i in 1:Np
        if k + i <= N
            r_preview[:, i] .= rs[:, k + i]
        else
            # Use last available reference
            r_preview[:, i] .= rs[:, end]
        end
    end
    
    return r_preview
end

function Simulation(mpc::Union{MPC,ExplicitMPC}; kwargs...)
    return Simulation((xk,uk,dk)-> mpc.model.F*xk +mpc.model.G*uk+mpc.model.Gd*dk, mpc; kwargs...)
end

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
