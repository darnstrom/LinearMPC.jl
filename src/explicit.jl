mutable struct ExplicitMPC 
    model::Model

    nr::Int
    nuprev::Int
    Np::Int

    solution::ParametricDAQP.Solution
    mpQP
    TH::NamedTuple{(:A,:b,:lb,:ub),Tuple{Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Float64}}}
    bst::Union{Nothing,ParametricDAQP.BinarySearchTree}
    settings::MPCSettings
    K::Matrix{Float64} # Prestabilizing feedback
    uprev::Vector{Float64}
end

function ExplicitMPC(mpc::MPC; range=nothing, build_tree=false, opts=ParametricDAQP.Settings(), single_soft=true)
    mpQP =  mpc2mpqp(mpc;singlesided=true,single_soft)
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
    nx,nr,nd,nuprev = get_parameter_dims(mpc)
    empc = ExplicitMPC(mpc.model, nr,nuprev, mpc.Np, sol, mpQP,TH,
                       nothing, mpc.settings, mpc.K, zeros(mpc.model.nu))

    # Build binary search tree
    build_tree && build_tree!(empc)

    return empc
end

function get_parameter_dims(mpc::ExplicitMPC)
    return mpc.model.nx, mpc.nr, mpc.model.nd, mpc.nuprev
end

function build_tree!(mpc::ExplicitMPC; daqp_settings=nothing, clipping=true)
    mpc.bst  = ParametricDAQP.build_tree(mpc.solution;daqp_settings, clipping)
    for i in 1:length(mpc.bst.feedbacks) # Correct for prestabilizing feedback
        mpc.bst.feedbacks[i][1:mpc.model.nx,:] -= mpc.K'
    end
    return mpc.bst
end

function plot_regions(mpc::ExplicitMPC;fix_ids=nothing,fix_vals=nothing,opts=Dict{Symbol,Any}())
    ParametricDAQP.plot_regions(mpc.solution;fix_ids,fix_vals,opts)
end

function plot_feedback(mpc::ExplicitMPC;u_id=1,fix_ids=nothing,fix_vals=nothing,opts=Dict{Symbol,Any}())
    push!(opts,:zlabel=>"\\large\$u_{"*string(u_id)*"}\$")
    ParametricDAQP.plot_solution(mpc.solution;z_id=u_id,fix_ids,fix_vals,opts)
end

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

function plot_regions(mpc::ExplicitMPC,lth1,lth2; show_fixed=true, show_zero = false,
        x=nothing, r=nothing, d=nothing, uprev=nothing)
    lth1,lth2 = Symbol(lth1),Symbol(lth2)

    x= isnothing(x) ? zeros(mpc.model.nx) : x
    r= isnothing(r) ? zeros(mpc.nr) : r
    d= isnothing(d) ? zeros(mpc.model.nd) : d
    uprev= isnothing(uprev) ? zeros(mpc.nuprev) : uprev

    free_ids,fix_ids,lx,ly= get_parameter_plot(mpc,lth1,lth2)
    fix_vals = [x;r;d;uprev][fix_ids]

    opts = Dict{Symbol,Any}()
    push!(opts,:xlabel=>latexify(lx))
    push!(opts,:ylabel=>latexify(ly))
    show_fixed && push!(opts,:title=>fixed_parameters_string(mpc,fix_ids,fix_vals;show_zero))

    ParametricDAQP.plot_regions(mpc.solution;fix_ids,fix_vals,opts)
end

function plot_regions(mpc::ExplicitMPC,ths; show_fixed=true, show_zero = false,
        x=nothing, r=nothing, d=nothing, uprev=nothing)
    plot_region(mpc,first(ths),last(ths):x,r,d,uprev,show_fixed,show_zero)
end

function plot_feedback(mpc::ExplicitMPC,lu1,lth1,lth2; show_fixed=true, show_zero = false,
        x=nothing, r=nothing, d=nothing, uprev=nothing)
    lu1,lth1,lth2 = Symbol(lu1),Symbol(lth1),Symbol(lth2)

    x= isnothing(x) ? zeros(mpc.model.nx) : x
    r= isnothing(r) ? zeros(mpc.nr) : r
    d= isnothing(d) ? zeros(mpc.model.nd) : d
    uprev= isnothing(uprev) ? zeros(mpc.nuprev) : uprev

    u_id = findfirst(x->x==lu1,mpc.model.labels.u)
    isnothing(u_id) && throw(ArgumentError("Unknown control $lu1"))
    lz = make_subscript(string(mpc.model.labels.u[u_id]))

    free_ids,fix_ids,lx,ly= get_parameter_plot(mpc,lth1,lth2)
    fix_vals = [x;r;d;uprev][fix_ids]

    fixed_parameters_string(mpc,fix_ids,fix_vals)

    opts = Dict{Symbol,Any}()
    push!(opts,:xlabel=>latexify(lx))
    push!(opts,:ylabel=>latexify(ly))
    push!(opts,:zlabel=>latexify(lz))

    show_fixed && push!(opts,:title=>fixed_parameters_string(mpc,fix_ids,fix_vals;show_zero))

    ParametricDAQP.plot_solution(mpc.solution;z_id=u_id,fix_ids,fix_vals,opts)
end

function plot_feedback(mpc::ExplicitMPC,lu1,ths; show_fixed=true, show_zero = false,
        x=nothing, r=nothing, d=nothing, uprev=nothing)
    plot_feedback(mpc,lu1,first(ths),last(ths);x,r,d,uprev,show_fixed, show_zero)
end

function fixed_parameters_string(mpc,fix_ids,fix_vals;show_zero = false)
    xlabels = [string(l) for l in mpc.model.labels.x]
    rlabels = [string(l)*"^r " for l in mpc.model.labels.y[1:mpc.nr]]
    dlabels = [string(l) for l in mpc.model.labels.d]
    ulabels = [string(l)*"^- " for l in mpc.model.labels.u[1:mpc.nuprev]]
    ls = [xlabels;rlabels;dlabels;ulabels]
    str = [latexify(make_subscript(ls[fix_ids[i]])*"="*string(fix_vals[i]))*", "
           for i in 1:length(fix_ids) if show_zero || fix_vals[i] != 0 ]
    return isempty(str) ? "" : "\\tiny "*reduce(*,str)[1:end-2]
end
