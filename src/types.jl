# lb <= Au uk + Ax xk <= ub for k ∈ ks
# (additional terms Ar rₖ, Aw wₖ, Ad dₖ, Aup u⁻ₖ)

struct Constraint
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
struct MPCWeights
    Q::Matrix{Float64}
    R::Matrix{Float64}
    Rr::Matrix{Float64}
    S::Matrix{Float64}
    Qf::Matrix{Float64}
    Qfx::Matrix{Float64}
end

function MPCWeights(nu,nx,nr)
    return MPCWeights(Matrix{Float64}(I,nr,nr),Matrix{Float64}(I,nu,nu),zeros(nu,nu),
                      zeros(nx,nu),zeros(nr,nr),zeros(nx,nx))
end

"""
MPC controller settings.

# Fields
- `QP_double_sided::Bool = true`: Use double-sided QP formulation
- `reference_condensation::Bool = false`: Collapse reference trajectory to setpoint 
- `reference_tracking::Bool = true`: Enable reference tracking
- `reference_preview::Bool = false`: Enable time-varying reference preview
- `soft_constraints::Bool = true`: Allow soft constraint violations
- `explicit_soft::Bool = false`: Use explicit slack variables for soft constraints
- `soft_weight::Float64 = 1e6`: Penalty weight for soft constraint violations
- `solver_opts::Dict{Symbol,Any}`: Additional solver options
"""
Base.@kwdef mutable struct MPCSettings
    QP_double_sided::Bool = true 
    reference_condensation::Bool= false
    reference_tracking::Bool= true
    reference_preview::Bool = false
    soft_constraints::Bool= true
    explicit_soft::Bool= false
    soft_weight::Float64= 1e6
    solver_opts::Dict{Symbol,Any} = Dict()
end

# MPC controller
mutable struct MPC 

    model::Model

    # parameters 
    nr::Int
    nuprev::Int

    # Horizons 
    Np::Int # Prediction
    Nc::Int # Control

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
    opt_model::DAQPBase.Model

    # Prestabilizing feedback
    K::Matrix{Float64}

    # Move blocks
    move_blocks::Vector{Int}

    mpqp_issetup::Bool

    uprev::Vector{Float64}

    traj2setpoint::Matrix{Float64}
end

function MPC(model::Model;Np=10,Nc=Np)
    MPC(model,0,0,Np,Nc,
        MPCWeights(model.nu,model.nx,model.ny),
        zeros(0),zeros(0),zeros(0),
        Constraint[],MPCSettings(),nothing,
        DAQP.Model(),zeros(model.nu,model.nx),Int[],false, zeros(model.nu),zeros(0,0))
end

function MPC(F,G;Gd=zeros(0,0), C=zeros(0,0), Dd= zeros(0,0), Ts= -1.0, Np=10, Nc = Np)
    MPC(Model(F,G;Gd,C,Dd,Ts);Np,Nc);
end

function MPC(A,B,Ts::Float64; Bd = zeros(0,0), C = zeros(0,0), Dd = zeros(0,0), Np=10, Nc=Np)
    MPC(Model(A,B,Ts;Bd,C,Dd);Np,Nc)
end

function MPC(sys; Ts=1.0, Np=10, Nc=Np)
    MPC(Model(sys;Ts);Np,Nc)
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


function ParameterRange(mpc::MPC)

    nx,nr,nd,nuprev = get_parameter_dims(mpc);

    xmin,xmax = -100*ones(nx),100*ones(nx)
    rmin,rmax = -100*ones(nr),100*ones(nr)
    dmin,dmax = -100*ones(nd),100*ones(nd)
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
