mutable struct ExplicitMPC 
    nx::Int
    nu::Int
    ny::Int
    nth::Int 
    solution::ParametricDAQP.Solution
    mpQP
    TH
    bst
    settings::MPCSettings
end

function ExplicitMPC(mpc::MPC; range=nothing, build_tree=false)
    mpQP = mpc.mpQP
    if isnothing(mpQP)
        mpQP = mpc2mpqp(mpc)
    end
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
        if(mpc.settings.explicit_soft && any(sense.==DAQP.SOFT))# Correct sign for slack
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
    # Build binary search tree
    bst = build_tree ? ParametricDAQP.build_tree(sol) : nothing 
    return ExplicitMPC(mpc.nx,mpc.nu,mpc.ny,length(TH.ub),
                       sol,mpQP, TH, bst,mpc.settings)
end


function plot_regions(mpc::ExplicitMPC;fix_ids=nothing,fix_vals=nothing,opts=Dict{Symbol,Any}())
    ParametricDAQP.plot_regions(mpc.solution;fix_ids,fix_vals,opts)
end

function plot_feedback(mpc::ExplicitMPC;u_id=0,fix_ids=nothing,fix_vals=nothing,opts=Dict{Symbol,Any}())
    ParametricDAQP.plot_solution(mpc.solution;z_id=u_id,fix_ids,fix_vals,opts)
end
