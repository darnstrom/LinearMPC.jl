mutable struct ExplicitMPC 
    nx::Int
    nu::Int
    ny::Int

    nr::Int
    nw::Int
    nd::Int
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
    mpQP = merge(mpQP,(out_inds=1:mpc.nu,)) # Only compute control at first time step
    # Compute mpQP solution
    sol,info = ParametricDAQP.mpsolve(mpQP, TH)
    nx,nr,nw,nd,nuprev = get_parameter_dims(mpc)
    empc = ExplicitMPC(mpc.nx,mpc.nu,mpc.ny,nr,nw,nd,nuprev,
                       sol,mpQP, TH, nothing,mpc.settings,mpc.K)

    # Build binary search tree
    build_tree && build_tree!(empc)

    return empc
end

function get_parameter_dims(mpc::ExplicitMPC)
    return mpc.nx, mpc.nr, mpc.nw, mpc.nd, mpc.nuprev
end

function build_tree!(mpc::ExplicitMPC)
    mpc.bst  = ParametricDAQP.build_tree(mpc.solution)
    for i in 1:length(mpc.bst.feedbacks) # Correct for prestabilizing feedback
        mpc.bst.feedbacks[i][1:mpc.nx,:] -= mpc.K'
    end
end

function plot_regions(mpc::ExplicitMPC;fix_ids=nothing,fix_vals=nothing,opts=Dict{Symbol,Any}())
    ParametricDAQP.plot_regions(mpc.solution;fix_ids,fix_vals,opts)
end

function plot_feedback(mpc::ExplicitMPC;u_id=1,fix_ids=nothing,fix_vals=nothing,opts=Dict{Symbol,Any}())
    push!(opts,:zlabel=>"\\large\$u_{"*string(u_id)*"}\$")
    ParametricDAQP.plot_solution(mpc.solution;z_id=u_id,fix_ids,fix_vals,opts)
end
