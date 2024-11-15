# lb <= Au uk + Ax xk <= ub for k ∈ ks 
mutable struct Constraint
    Au::Matrix{Float64}
    Ax::Matrix{Float64}
    ub::Vector{Float64}
    lb::Vector{Float64}
    ks::AbstractVector{Int64}
    soft::Bool
    binary::Bool
end

# Constraints for the OCP
mutable struct MPCConstraints

end

function MPCConstraints(nx,nu)
    return MPCConstraints(zeros(0),zeros(0),0,zeros(0),Constraint[])
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
    return MPCWeights(Matrix{Float64}(I,nr,nr),Matrix{Float64}(I,nu,nu),zeros(nu,nu),1e6,nothing)
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
    F::Union{Matrix{Float64},Vector{Matrix{Float64}}}
    G::Union{Matrix{Float64},Vector{Matrix{Float64}}}
    Ts::Float64

    # Dims
    nx::Int
    nu::Int

    # Horizons 
    Np::Int # Prediction
    Nc::Int # Control
    Nb::Int # Bound

    C::Matrix{Float64}

    ## 
    weights::MPCWeights

    # lb <= u <=ub
    umin::Vector{Float64}
    umax::Vector{Float64}
    binary_controls::Vector{Int64}

    # General constraints 
    constraints::Vector{Constraint}

    # Settings
    settings::MPCSettings

    #Optimization problem
    mpQP
end

function MPC(F,G,C,Np; Nc = Np, Nb = Nc)
    nx,nu = size(G);
    nr = size(C,1);
    MPC(F,G,1,nx,nu,
        Np,Nc,Nc,C,MPCWeights(nu,nr),
        zeros(0),zeros(0),zeros(0),
        Constraint[],MPCSettings(),nothing);
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
