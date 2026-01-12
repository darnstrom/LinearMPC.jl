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
function certify(mpc::MPC; range=nothing, AS0 = Int[], settings = ASCertain.CertSettings(), single_soft=true)
    mpQP = make_singlesided(mpc2mpqp(mpc);single_soft,soft_weight=mpc.settings.soft_weight)
    if isnothing(range)
        @warn("No parameter range defined. Using default limits [-100 and 100]."* 
              "If you want a bigger/smaller region, create a ParameterRange")
        range = ParameterRange(mpc)
    end

    region = range2region(range)
    part, iter_max = ASCertain.certify(mpQP,region,convert(Vector{Int},AS0);opts=settings)
    return CertificationResult(mpc,iter_max,part,region)
end
## Plots.jl
@recipe function f(cert::CertificationResult; parameters=[], show_fixed=true, show_zero=false, 
        x=nothing,r=nothing,d=nothing,uprev=nothing)
    # Make sure these are removed from plots directory Plots.jl arguments 
    r = pop!(plotattributes, :rotation,r)
    x = pop!(plotattributes, :x,x)
    plotattributes[:xrotation] = 0
    plotattributes[:yrotation] = 0
    plotattributes[:zrotation] = 0

    x= isnothing(x) ? zeros(cert.mpc.model.nx) : x
    r= isnothing(r) ? zeros(cert.mpc.nr) : r
    d= isnothing(d) ? zeros(cert.mpc.model.nd) : d
    uprev= isnothing(uprev) ? zeros(cert.mpc.nuprev) : uprev

    parameters = Symbol.(parameters)
    length(parameters) != 2 && throw(ArgumentError("parameters has to contain two elements"))

    free_ids,fix_ids,lx,ly= get_parameter_plot(cert.mpc,parameters[1],parameters[2])
    fix_vals = [x;r;d;uprev][fix_ids]

    xlabel --> latexify(lx)
    ylabel --> latexify(ly)

    if show_fixed
        colorbar_title --> fixed_parameters_string(cert.mpc,fix_ids,fix_vals;show_zero)
    end
    plotattributes[:CR_attr] = (0,free_ids,fix_ids,fix_vals)
    return cert.partition
end
