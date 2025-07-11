struct Simulation
    ts::Vector{Float64}
    ys::Matrix{Float64}
    us::Matrix{Float64}
    xs::Matrix{Float64}
    rs::Matrix{Float64}
    ds::Matrix{Float64}

    mpc::Union{MPC,ExplicitMPC}
end

function Simulation(dynamics, mpc::Union{MPC,ExplicitMPC}; x0=zeros(mpc.model.nx),N=1000, r=nothing,d=nothing, callback=(x,u,d,k)->nothing)
    x,u = x0,zeros(mpc.model.nu)
    
    xs = zeros(mpc.model.nx,N); xs[:,1] = x0
    ys = zeros(mpc.model.ny,N);
    rs = zeros(mpc.model.ny,N);
    ds = zeros(mpc.model.nd,N);
    us = zeros(mpc.model.nu,N)

    # Setup reference 
    if(!isnothing(r))
        rs[:,1:size(r,2)].= r
        rs[:,size(r,2)+1:end] .= r[:,end] # hold last value of reference
    end

    # Setup disturbance
    if(!isnothing(d))
        ds[:,1:size(d,2)].= d
        ds[:,size(d,2)+1:end] .= d[:,end] # hold last
    end

    # Start the simulation
    for k = 1:N
        xs[:,k], ys[:,k] = x, mpc.model.C*x+mpc.model.Dd*ds[:,k]
        
        # Prepare reference for controller
        if !isnothing(r)
            # Reference preview mode: provide future references
            r_preview = get_reference_preview(rs, k, mpc.Np)
            u = compute_control(mpc,x;r=r_preview,d=ds[:,k])
        else
            # Standard mode: provide current reference
            u = compute_control(mpc,x;r=rs[:,k],d=ds[:,k])
        end
        
        x = dynamics(x,u,ds[:,k])
        callback(x,u,ds[:,k],k)
        us[:,k] = u
    end
    Ts = mpc.model.Ts < 0.0 ? 1.0 : mpc.model.Ts
    return Simulation(collect(Ts*(0:1:N-1)),ys,us,xs,rs,ds,mpc)
end

"""
    get_reference_preview(rs, k, Np)

Extract reference preview from reference trajectory starting at time step k.
"""
function get_reference_preview(rs, k, Np)
    ny, N = size(rs)
    r_preview = zeros(ny, Np)
    
    @views for i in 1:Np
        if k + i - 1 <= N
            r_preview[:, i] .= rs[:, k + i - 1]
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

@recipe function f(sim::Simulation; yids=1:sim.mpc.model.ny,uids=sim.mpc.model.nu,xids=[])


    layout = Tuple{Int,Int}[]
    length(yids) == 0 || push!(layout,(length(yids), 1))
    length(uids) == 0 || push!(layout,(length(uids), 1))
    length(xids) == 0 || push!(layout,(length(xids), 1))
    layout := reshape(layout,1,length(layout))

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
            legend  --> true
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
        if(length(sim.mpc.umin) > i && sim.mpc.umin[i] > -1e12)
            @series begin
                color     := :black 
                subplot   --> id
                linestyle --> :dash
                linewidth --> 1.0 
                primary   --> false
                sim.ts[[1,end]], fill(sim.mpc.umin[i],2) 
            end
        end

        # upper bound
        if(length(sim.mpc.umax) > i &&sim.mpc.umax[i] < 1e12)
            @series begin
                color     := :black 
                subplot   --> id
                linestyle --> :dash
                linewidth --> 1 
                primary   --> false
                sim.ts[[1,end]], fill(sim.mpc.umax[i],2) 
            end
        end
        @series begin
            yguide--> latexify(make_subscript(string(sim.mpc.model.labels.u[i])))
            color      --> 1
            subplot    --> id
            seriestype --> :steppost
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
            legend --> true
            if i == length(xids)
                xguide --> (sim.mpc.model.Ts < 0 ? "Time step" : "Time [s]")
            end
            sim.ts, sim.xs[i, :]
        end
        id+=1
    end
end
