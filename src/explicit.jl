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

## Parameter utils 
using Latexify
function get_parameter_plot(mpc,l1,l2)
    id1,l1 = label2id(mpc,l1)
    isnothing(id1) && throw(ArgumentError("Unknown label $l1"))  
    id2,l2 = label2id(mpc,l2)
    isnothing(id2) && throw(ArgumentError("Unknown label $l2"))
    # lower id is always plotted on x-axis
    if id1 < id2
        lx,ly,idx,idy = l1,l2,id1,id2
    else
        lx,ly,idx,idy = l2,l1,id2,id1
    end
    fix_ids = setdiff(1:sum(get_parameter_dims(mpc)),[id1,id2])
    return [idx,idy],fix_ids,make_subscript(lx),make_subscript(ly)
end

function fixed_parameters_string(mpc,fix_ids,fix_vals;show_zero = false)
    xlabels = [string(l) for l in mpc.model.labels.x]
    rlabels = [string(l)*"^r " for l in mpc.model.labels.y[1:mpc.nr]]
    dlabels = [string(l) for l in mpc.model.labels.d]
    ulabels = [string(l)*"^- " for l in mpc.model.labels.u[1:mpc.nuprev]]
    ls = [xlabels;rlabels;dlabels;ulabels]
    str = [latexify(make_subscript(ls[fix_ids[i]])*"="*string(fix_vals[i]))*", "
           for i in 1:length(fix_ids) if show_zero || fix_vals[i] != 0 ]
    return isempty(str) ? "" : reduce(*,str)[1:end-2]
end

## Plots.jl
using RecipesBase
@recipe function f(mpc::ExplicitMPC; parameters=[],control=nothing, show_fixed=true, show_zero=false, 
        x=nothing,r=nothing,d=nothing,uprev=nothing)

    # Make sure these are removed from plots directory Plots.jl arguments 
    r = pop!(plotattributes, :rotation,r)
    x = pop!(plotattributes, :x,x)
    plotattributes[:xrotation] = 0
    plotattributes[:yrotation] = 0
    plotattributes[:zrotation] = 0

    x= isnothing(x) ? zeros(mpc.model.nx) : x
    r= isnothing(r) ? zeros(mpc.nr) : r
    d= isnothing(d) ? zeros(mpc.model.nd) : d
    uprev= isnothing(uprev) ? zeros(mpc.nuprev) : uprev

    parameters,control = Symbol.(parameters),Symbol(control)
    length(parameters) != 2 && throw(ArgumentError("parameters has to contain two elements"))
    u_id = something(findfirst(x->x==control,mpc.model.labels.u),0)

    free_ids,fix_ids,lx,ly= get_parameter_plot(mpc,parameters[1],parameters[2])
    fix_vals = [x;r;d;uprev][fix_ids]

    xlabel --> latexify(lx)
    ylabel --> latexify(ly)
    if u_id != 0
        zlabel --> latexify(make_subscript(string(mpc.model.labels.u[u_id])))
    end

    if show_fixed
        titlefontsize --> 10
        title --> fixed_parameters_string(mpc,fix_ids,fix_vals;show_zero)
    end
    plotattributes[:CR_attr] = (u_id ,free_ids,fix_ids,fix_vals)
    return ParametricDAQP.get_critical_regions(mpc.solution)
end
