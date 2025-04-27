struct CertificationResult
    mpc::MPC
    max_iterations::Int
    partition::Vector{ASCertain.Region}
    region::NamedTuple{(:A,:b,:lb,:ub),Tuple{Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Float64}}}
end

Base.:display(result::CertificationResult ) = display(result.partition)

"""
    certify(mpc; range, AS0, opts)
    
Provide certificates on the iteration complexity of DAQP for solving the resulting optimization problems. 
* `range` is the parameter range over which the certification should be done
* `AS0` is the starting working set in DAQP (defaults to empty)
* `settings` the settings used in the certification (see ASCertain.CertSettings()) 
"""
function certify(mpc::MPC; range=nothing, AS0 = Int[], settings = nothing)
    isnothing(mpc.mpQP) && setup!(mpc) # ensure mpQP is setup 
    settings = isnothing(settings) ? ASCertain.CertSettings() : settings

    if mpc.settings.QP_double_sided
        mpQP = make_singlesided(mpc.mpQP; explicit_soft = mpc.settings.explicit_soft)
    else
        mpQP = mpc.mpQP
    end

    if isnothing(range)
        @warn("No parameter range defined. Using default limits [-100 and 100]."* 
              "If you want a bigger/smaller region, create a ParameterRange")
        range = ParameterRange(mpc)
    end

    region = range2region(range)
    part, iter_max = ASCertain.certify(mpQP,region,convert(Vector{Int},AS0);opts=settings)
    return CertificationResult(mpc,iter_max,part,region)
end

function plot(cert::CertificationResult,lth1,lth2; show_fixed=true, show_zero = false, 
        x=nothing, r=nothing, d=nothing, uprev=nothing)
    lth1,lth2 = Symbol(lth1),Symbol(lth2)

    x= isnothing(x) ? zeros(cert.mpc.model.nx) : x
    r= isnothing(r) ? zeros(cert.mpc.nr) : r
    d= isnothing(d) ? zeros(cert.mpc.model.nd) : d
    uprev= isnothing(uprev) ? zeros(cert.mpc.nuprev) : uprev

    free_ids,fix_ids,lx,ly= get_parameter_plot(cert.mpc,lth1,lth2)
    fix_vals = [x;r;d;uprev][fix_ids]

    opts = Dict{Symbol,Any}()
    push!(opts,:xlabel=>latexify(lx))
    push!(opts,:ylabel=>latexify(ly))
    show_fixed && push!(opts,:title=>fixed_parameters_string(cert.mpc,fix_ids,fix_vals;show_zero))

    push!(opts, :xticklabels=>"{$(cert.region.lb[free_ids[1]]),$(cert.region.ub[free_ids[1]])}")
    push!(opts, :yticklabels=>"{$(cert.region.lb[free_ids[2]]),$(cert.region.ub[free_ids[2]])}")

    ASCertain.pplot(cert.partition;key=nothing, fix_ids, fix_vals, opts)
end
