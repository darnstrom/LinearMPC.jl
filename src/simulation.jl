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
        u = compute_control(mpc,x;r=rs[:,k],d=ds[:,k],uprev=u)
        x = dynamics(x,u,ds[:,k])
        callback(x,u,ds[:,k],k)
        us[:,k] = u
    end
    return Simulation(collect(mpc.model.Ts*(0:1:N-1)),ys,us,xs,rs,ds,mpc)
end

function Simulation(mpc::Union{MPC,ExplicitMPC}; kwargs...)
    return Simulation((xk,uk,dk)-> mpc.model.F*xk +mpc.model.G*uk+mpc.model.Gd*dk, mpc; kwargs...)
end

using RecipesBase

@recipe function f(sim::Simulation; show_y=true,show_u=true,show_x=false)

    ny = show_y ?  sim.mpc.model.ny : 0 
    nu = show_u ?  sim.mpc.model.nu : 0
    nx = show_x ?  sim.mpc.model.nx : 0

    layout = Tuple{Int,Int}[]
    ny == 0 || push!(layout,(ny, 1))
    nu == 0 || push!(layout,(nu, 1))
    nx == 0 || push!(layout,(nx, 1))
    layout := reshape(layout,1,length(layout))

    # Plot y
    id = 1 
    for i in 1:ny
        @series begin
            yguide  --> sim.mpc.model.labels.y[i]
            color   --> 1
            subplot --> id 
            label--> latexify(make_subscript(string(sim.mpc.model.labels.y[i])))
            legend  --> true
            sim.ts, sim.ys[i, :]
        end
        @series begin
            color     --> 2
            subplot   --> id
            linestyle --> :dash
            linewidth --> 0.5
            label     --> "reference"
            primary   --> false
            i == ny && (xguide --> "Time [s]")
            sim.ts, sim.rs[i, :]
        end

        id+=1
    end

    # Plot u 
    for i in 1:nu
        @series begin
            yguide     --> sim.mpc.model.labels.u[i]
            color      --> 1
            subplot    --> id
            seriestype --> :steppost
            label--> latexify(make_subscript(string(sim.mpc.model.labels.u[i])))
            legend     --> true
            sim.ts, sim.us[i, :]
        end
        # lower bound
        if(sim.mpc.umin[i] > -1e12)
            @series begin
                color     --> 3 
                subplot   --> id
                linestyle --> :dash
                linewidth --> 1.0 
                primary   --> false
                sim.ts, fill(sim.mpc.umin[i],length(sim.ts)) 
            end
        end

        # upper bound
        if(sim.mpc.umax[i] < 1e12)
            @series begin
                color     --> 3 
                subplot   --> id
                linestyle --> :dash
                linewidth --> 1 
                primary   --> false
                i == nu && (xguide --> "Time [s]")
                sim.ts, fill(sim.mpc.umax[i],length(sim.ts)) 
            end
        end
        id+=1
    end

    ## Plot x 
    for i in 1:nx
        @series begin
            label  --> latexify(make_subscript(string(sim.mpc.model.labels.x[i])))
            color  --> 1
            subplot--> id
            legend --> true
            i == nx && (xguide --> "Time [s]")
            sim.ts, sim.xs[i, :]
        end
        id+=1
    end
end
