# lb <= Au uk + Ax xk <= ub for k ∈ ks
# (additional terms Ar rₖ, Aw wₖ, Ad dₖ, Aup u⁻ₖ)

mutable struct Constraint
    Au::Matrix{Float64}
    Ax::Matrix{Float64}
    Ar::Matrix{Float64}
    Aw::Matrix{Float64}
    Ad::Matrix{Float64}
    Aup::Matrix{Float64}
    ub::Vector{Float64}
    lb::Vector{Float64}
    ks::AbstractVector{Int64}
    soft::Bool
    binary::Bool
    prio::Int
end

# Weights used to define the objective function of the OCP
mutable struct MPCWeights
    Q::Matrix{Float64}
    R::Matrix{Float64}
    Rr::Matrix{Float64}
    S::Matrix{Float64}
    rho::Float64
    Qf::Matrix{Float64}
end

function MPCWeights(nu,nx,nr)
    return MPCWeights(Matrix{Float64}(I,nr,nr),Matrix{Float64}(I,nu,nu),zeros(nu,nu),zeros(nx,nu),1e6,zeros(0,0))
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

    # Disturbance model
    Gd::Matrix{Float64}
    Dd::Matrix{Float64}

    # Dims
    nx::Int
    nu::Int
    ny::Int

    nr::Int
    nd::Int
    nuprev::Int

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

    # DAQP optimization model
    opt_model 

    # Prestabilizing feedback
    K::Matrix{Float64}

    # Move blocks
    move_blocks::Vector{Int}
end

function MPC(F,G;
        C=nothing,Np=10, Nc = Np, Nb = Nc, Ts = 1.0, Gd = nothing, Dd = nothing)
    G = reshape(G,size(G,1),:) 
    nx,nu = size(G)
    C = isnothing(C) ? Matrix{Float64}(I,nx,nx) : float(C)
    ny = size(C,1);
    (size(C,2)==size(F,1)==nx) || throw(ArgumentError("Dimensions of ss-model incompatible"))
    # disturbance
    Gd = isnothing(Gd) ? zeros(nx,0) : Gd
    Dd = isnothing(Dd) ? zeros(ny,0) : Dd 
    nd = max(size(Gd,2),size(Dd,2))
    Gd = [Gd zeros(nx,nd-size(Gd,2))]
    Dd = [Dd zeros(ny,nd-size(Dd,2))]

    MPC(float(F),float(G),Ts,Gd,Dd,
        nx,nu,ny,0,nd,0,
        Np,Nc,Nc,C,MPCWeights(nu,nx,ny),
        zeros(0),zeros(0),zeros(0),
        Constraint[],MPCSettings(),nothing,nothing,zeros(nu,nx),Int[])
end

function MPC(A,B, Ts::Float64;
        C=nothing,Np=10, Nc = Np, Nb = Nc, Bd = nothing)
    (size(A,1)==size(B,1)) || throw(ArgumentError("Dimensions of ss-model incompatible"))
    F,G,Gd=zoh(A,B,Ts;Bd)
    MPC(F,G;C,Np,Nc,Nb,Ts,Gd)
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

struct ParameterRange
    xmin::Vector{Float64}
    xmax::Vector{Float64}

    rmin::Vector{Float64}
    rmax::Vector{Float64}

    dmin::Vector{Float64}
    dmax::Vector{Float64}

    umin::Vector{Float64}
    umax::Vector{Float64}
end

function ParameterRange(mpc)

    nx,nr,nd,nuprev = get_parameter_dims(mpc);

    # states
    xmin,xmax = -100*ones(nx), 100*ones(nx)

    # references 
    rmin,rmax = -100*ones(nr),100*ones(nr)

    # constant disturbance
    dmin,dmax = -100*ones(nd),100*ones(nd)

    # previous control (for Δu)
    if(nuprev > 0)
        nmin,nmax = length(mpc.umin),length(mpc.umax)
        nb = max(nmin,nmax)
        umin = [mpc.umin;-100*ones(nb-nmin)]
        umax = [mpc.umax;+100*ones(nb-nmax)]
    else
        umin,umax = zeros(0),zeros(0)
    end

    return ParameterRange(xmin,xmax,
                          rmin,rmax,
                          dmin,dmax,
                          umin,umax)
end
