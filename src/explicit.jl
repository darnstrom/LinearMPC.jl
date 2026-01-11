mutable struct ExplicitMPC
    model::Model

    nr::Int
    nuprev::Int
    nl::Int
    Np::Int
    Nc::Int
    move_blocks::Vector{Vector{Int}}

    solution::ParametricDAQP.Solution
    mpQP
    TH::NamedTuple{(:A,:b,:lb,:ub),Tuple{Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Float64}}}
    bst::Union{Nothing,ParametricDAQP.BinarySearchTree}
    settings::MPCSettings
    K::Matrix{Float64} # Prestabilizing feedback
    uprev::Vector{Float64}
    traj2setpoint::Matrix{Float64}
    state_observer
end

function ExplicitMPC(mpc::MPC; range=nothing, build_tree=false, opts=ParametricDAQP.Settings(), single_soft=true)
    mpQP = make_singlesided(mpc2mpqp(mpc);single_soft,soft_weight=mpc.settings.soft_weight)
    if mpQP.has_binaries
        @warn("Explicit controllers currently not supported for hybrid systems")
        return nothing
    end
    if(range==nothing)
        @warn("No parameter range defined. Using default limits [-100 and 100]. If you want a bigger/smaller region, create a ParameterRange")
        range = ParameterRange(mpc)
    end

    TH = range2region(range)

    mpQP = merge(mpQP,(out_inds=1:mpc.model.nu,)) # Only compute control at first time step
    # Compute mpQP solution
    opts.daqp_settings = merge(Dict(:sing_tol => 1e-11),opts.daqp_settings)
    sol,info = ParametricDAQP.mpsolve(mpQP, TH;opts)
    nx,nr,nd,nuprev,nl = get_parameter_dims(mpc)
    empc = ExplicitMPC(mpc.model, nr, nuprev, nl, mpc.Np, mpc.Nc, mpc.move_blocks, sol, mpQP, TH,
                       nothing, mpc.settings, mpc.K, zeros(mpc.model.nu), mpc.traj2setpoint, mpc.state_observer)

    # Build binary search tree
    build_tree && build_tree!(empc)

    return empc
end

function get_parameter_dims(mpc::ExplicitMPC)
    return mpc.model.nx, mpc.nr, mpc.model.nd, mpc.nuprev, mpc.nl
end

function form_parameter(mpc::Union{MPC,ExplicitMPC},x,r,d,uprev,l=nothing)
    # Setup parameter vector θ
    nx,nr,nd,nuprev,nl = get_parameter_dims(mpc)
    r = format_reference(mpc, r)
    d = isnothing(d) ? zeros(nd) : d
    uprev = isnothing(uprev) ? mpc.uprev[1:nuprev] : uprev[1:nuprev]
    l_vec = format_linear_cost(mpc, l)
    return [x;r;d;uprev;l_vec]
end

function build_tree!(mpc::ExplicitMPC; daqp_settings=nothing, clipping=true)
    mpc.bst  = ParametricDAQP.build_tree(mpc.solution;daqp_settings, clipping)
    for i in 1:length(mpc.bst.feedbacks) # Correct for prestabilizing feedback
        mpc.bst.feedbacks[i][1:mpc.model.nx,:] -= mpc.K'
    end
    return mpc.bst
end
