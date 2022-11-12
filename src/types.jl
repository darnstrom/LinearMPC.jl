# Constraints for the OCP
mutable struct MPCConstraints
    # lb <= u <=ub
    lb::Vector{Float64}
    ub::Vector{Float64}
    Ncc::Int64
    binary_controls::Vector{Int64}
    #sense_cc::Vector{Int64}

    # lby <= Cy x_i <= uby for i ∈ 1,2, …, Ncy
    Cy::Vector{Matrix{Float64}}
    lby::Vector{Vector{Float64}}
    uby::Vector{Vector{Float64}}
    Ncy::Vector{AbstractVector{Int64}} # TODO change name...
    #sense_cy::Vector{Int64}


    # Au u_i + Ax x_i <= bg for i ∈ 0,…,Ncy
    Au::Matrix{Float64}
    Ax::Matrix{Float64}
    bg::Vector{Float64}
    Ncg::Int64
    #sense_cg::Vector{Int64}
end

function MPCConstraints(nx,nu)
    return MPCConstraints(zeros(0),zeros(0),0,zeros(0),
                          Matrix{Float64}[],Vector{Float64}[],Vector{Float64}[],[0:0],
                          zeros(0,nu),zeros(0,nx),zeros(0),0)
end

# Weights used to define the objective function of the OCP
mutable struct MPCWeights
    Q::Matrix{Float64}
    R::Matrix{Float64}
    Rr::Matrix{Float64}
    rho::Float64
    Qf
end

function MPCWeights(nu,nr)
    return MPCWeights(Matrix{Float64}(I,nr,nr),Matrix{Float64}(I,nu,nu),Matrix{Float64}(I,nu,nu),1e6,nothing)
end

Base.@kwdef mutable struct MPCSettings
    QP_double_sided::Bool = true 
    reference_tracking::Bool= true
    soft_constraints::Bool= true
    explicit_soft::Bool= false
    solver_opts::Dict{Symbol,Any} = Dict()
end


# MPC controller
mutable struct MPC 
    # Plant
    F::Matrix{Float64}
    G::Matrix{Float64}
    Ts::Float64

    # Dims
    nx::Int64
    nu::Int64

    # OCP data
    Np::Int64
    Nc::Int64
    C::Matrix{Float64}
    weights::MPCWeights
    constraints::MPCConstraints

    # Settings
    settings::MPCSettings

    #Optimization problem
    mpQP
end

function MPC(F,G,C,N)
    nx,nu = size(G);
    nr = size(C,1);
    MPC(F,G,1,nx,nu,N,N,C,MPCWeights(nu,nr),MPCConstraints(nx,nu),MPCSettings(),nothing);
end


mutable struct MPQP
    H::Matrix{Float64}
    f::Matrix{Float64}
    f_theta::Matrix{Float64}
    H_theta::Matrix{Float64}

    A::Matrix{Float64}
    b::Matrix{Float64}
    W::Matrix{Float64}

    bounds_table::Vector{Int64}
    senses::Vector{Cint}
    MPQP()=new()
    MPQP(H,f,f_theta,H_theta,A,b,W,bounds_table,senses) = new(H,f,f_theta,H_theta,A,b,W,bounds_table,senses) 
end
