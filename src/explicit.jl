mutable struct ExplicitMPC 
    model::Model

    nr::Int
    nuprev::Int

    solution::ParametricDAQP.Solution
    mpQP
    TH
    bst
    settings::MPCSettings
    K::Matrix{Float64} # Prestabilizing feedback
end

function ExplicitMPC(mpc::MPC; range=nothing, build_tree=false)
    mpQP = isnothing(mpc.mpQP) ? mpc2mpqp(mpc) : mpc.mpQP
    if(range==nothing)
        @warn("No parameter range defined. Using default limits [-100 and 100]. If you want a bigger/smaller region, create a ParameterRange")
        range = ParameterRange(mpc)
    end

    TH = range2region(range)

    # Transform into single-sided QP 
    if(mpc.settings.QP_double_sided) 
        ncstr = length(mpQP.bu);
        n_bounds = ncstr-size(mpQP.A,1);
        bounds_table=[collect(ncstr+1:2*ncstr);collect(1:ncstr)]
        A = [I(n_bounds) zeros(n_bounds,size(mpQP.A,2)-n_bounds);mpQP.A]
        A = [A;-A]
        if(mpc.settings.explicit_soft && any(senses.==DAQP.SOFT))# Correct sign for slack
            A[:,end].= -abs.(A[:,end])
        end
        b = [mpQP.bu;-mpQP.bl]
        W = [mpQP.W;-mpQP.W]
        senses = [mpQP.senses;mpQP.senses]
        mpQP = (H=mpQP.H,f=mpQP.f, H_theta = mpQP.H_theta, f_theta=mpQP.f_theta,
                A=Matrix{Float64}(A), b=b, W=W, senses=senses,
                bounds_table=bounds_table)
    end
    mpQP = merge(mpQP,(out_inds=1:mpc.model.nu,)) # Only compute control at first time step
    # Compute mpQP solution
    sol,info = ParametricDAQP.mpsolve(mpQP, TH)
    nx,nr,nd,nuprev = get_parameter_dims(mpc)
    empc = ExplicitMPC(mpc.model,nr,nuprev, sol,mpQP, TH, nothing,mpc.settings,mpc.K,)

    # Build binary search tree
    build_tree && build_tree!(empc)

    return empc
end

function get_parameter_dims(mpc::ExplicitMPC)
    return mpc.model.nx, mpc.nr, mpc.model.nd, mpc.nuprev
end

function build_tree!(mpc::ExplicitMPC)
    mpc.bst  = ParametricDAQP.build_tree(mpc.solution)
    for i in 1:length(mpc.bst.feedbacks) # Correct for prestabilizing feedback
        mpc.bst.feedbacks[i][1:mpc.model.nx,:] -= mpc.K'
    end
end

function plot_regions(mpc::ExplicitMPC;fix_ids=nothing,fix_vals=nothing,opts=Dict{Symbol,Any}())
    ParametricDAQP.plot_regions(mpc.solution;fix_ids,fix_vals,opts)
end

function plot_feedback(mpc::ExplicitMPC;u_id=1,fix_ids=nothing,fix_vals=nothing,opts=Dict{Symbol,Any}())
    push!(opts,:zlabel=>"\\large\$u_{"*string(u_id)*"}\$")
    ParametricDAQP.plot_solution(mpc.solution;z_id=u_id,fix_ids,fix_vals,opts)
end

using Latexify

function get_parameter_plot(mpc::ExplicitMPC,l1,l2)
    id1,l1 = label2id(mpc,l1)
    isnothing(id1) && throw(ArgumentError("Unknown label $l1"))  
    id2,l2 = label2id(mpc,l2)
    isnothing(id2) && throw(ArgumentError("Unknown label $l2"))
    # lower id is always plotted on x-axis
    if id1 < id2
        lx,ly = l1,l2
    else
        lx,ly = l2,l1
    end
    fix_ids = setdiff(1:sum(get_parameter_dims(mpc)),[id1,id2])
    return [id1,id2],fix_ids,make_subscript(lx),make_subscript(ly)
end

function plot_regions(mpc::ExplicitMPC,lth1,lth2; show_fixed=true, show_zero = false,
        x=zeros(mpc.model.nx), r=zeros(mpc.nr), d=zeros(mpc.model.nd), uprev=zeros(mpc.nuprev))
    lth1,lth2 = Symbol(lth1),Symbol(lth2)

    free_ids,fix_ids,lx,ly= get_parameter_plot(mpc,lth1,lth2)
    fix_vals = [x;r;d;uprev][fix_ids]

    opts = Dict{Symbol,Any}()
    push!(opts,:xlabel=>latexify(lx))
    push!(opts,:ylabel=>latexify(ly))
    show_fixed && push!(opts,:title=>fixed_parameters_string(mpc,fix_ids,fix_vals;show_zero))

    ParametricDAQP.plot_regions(mpc.solution;fix_ids,fix_vals,opts)
end

function plot_regions(mpc::ExplicitMPC,ths; show_fixed=true, show_zero = false,
        x=zeros(mpc.model.nx), r=zeros(mpc.nr), d=zeros(mpc.model.nd), uprev=zeros(mpc.nuprev))
    plot_region(mpc,first(ths),last(ths):x,r,d,uprev,show_fixed,show_zero)
end

function plot_feedback(mpc::ExplicitMPC,lu1,lth1,lth2; show_fixed=true, show_zero = false,
        x=zeros(mpc.model.nx), r=zeros(mpc.nr), d=zeros(mpc.model.nd), uprev=zeros(mpc.nuprev))
    lu1,lth1,lth2 = Symbol(lu1),Symbol(lth1),Symbol(lth2)

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
        x=zeros(mpc.model.nx), r=zeros(mpc.nr), d=zeros(mpc.model.nd), uprev=zeros(mpc.nuprev))
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
