function codegen(mpc::MPC;fname="mpc_workspace", dir="codegen", opt_settings=nothing, src=true)
    length(dir)==0 && (dir="codegen")
    dir[end] != '/' && (dir*="/") ## Make sure it is a correct directory path
    ## Generate mpQP
    mpc.settings.QP_double_sided = true # force double-sided constraint for generated code
    setup!(mpc)
    # Generate QP workspace
    d = mpc.opt_model
    if(!isnothing(opt_settings))
        DAQP.settings(d,opt_settings)
    end
    # TODO: move this to DAQP.jl
    DAQP.codegen(d;fname,dir,src)

    # Append MPC-specific data/functions
    render_mpc_workspace(mpc;fname,dir, fmode="a")

    @info "Generated code for MPC controller" dir fname
end

function codegen(mpc::ExplicitMPC;fname="empc", dir="codegen", opt_settings=nothing, src=true,float_type="float")
    length(dir)==0 && (dir="codegen")
    dir[end] != '/' && (dir*="/") ## Make sure it is a correct directory path
    # Generate code for explicit solution
    ParametricDAQP.codegen(mpc.solution;dir,fname,float_type)

    # Generate code for MPC
    nth = mpc.nx+mpc.nr+mpc.nw+mpc.nuprev

    # HEADER
    fh = open(joinpath(dir,"mpc_compute_control.h"), "w")
    hguard = "MPC_COMPUTE_CONTROL_H"
    @printf(fh, "#ifndef %s\n",   hguard);
    @printf(fh, "#define %s\n\n", hguard);

    write(fh, "typedef $float_type c_float;\n")
    @printf(fh, "#define N_STATE %d\n",mpc.nx);
    @printf(fh, "#define N_REFERENCE %d\n",mpc.nr);
    @printf(fh, "#define N_IN_DISTURBANCE %d\n",mpc.nw);
    @printf(fh, "#define N_OUT_DISTURBANCE %d\n",mpc.nd);
    @printf(fh, "#define N_CONTROL_PREV %d\n",mpc.nuprev);

    @printf(fh, "extern c_float mpc_parameter[%d];\n", nth);

    @printf(fh, "int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* in_dist, c_float* out_dist); \n");
    @printf(fh, "#endif // ifndef %s\n", hguard);
    close(fh)

    # SOURCE
    fsrc = open(joinpath(dir,"mpc_compute_control.c"), "w")
    write(fsrc, """
#include "mpc_compute_control.h"
#include "$fname.h"

int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* in_dist, c_float* out_dist){
    int i,j;
    c_float mpc_parameter[$nth];
    // update parameter
    for(i=0,j=0;j<N_STATE;i++, j++) mpc_parameter[i] = state[j];
    for(j=0;j<N_REFERENCE;i++, j++) mpc_parameter[i] = reference[j];
    for(j=0;j<N_IN_DISTURBANCE;i++, j++) mpc_parameter[i] = in_dist[j];
    for(j=0;j<N_OUT_DISTURBANCE;i++, j++) mpc_parameter[i] = out_dist[j];
    for(j=0;j<N_CONTROL_PREV;i++, j++) mpc_parameter[i] = control[j];

    // Get the solution at the parameter
    $(fname)_evaluate(mpc_parameter,control);

    return 1;
}
          """)
    close(fsrc)

    @info "Generated code for EMPC controller" dir fname
end

function render_mpc_workspace(mpc;fname="mpc_workspace",dir="",fmode="w")
    mpLDP = qp2ldp(mpc.mpQP,mpc.nu) 
    mpLDP.Xth[1:mpc.nx,:] -= mpc.K' #Account for prestabilizing feedback
    # Get dimensions
    nth,m = size(mpLDP.Dth)

    # Setup files
    fh = open(dir*fname*".h", fmode)
    fsrc = open(dir*fname*".c", fmode)

    # HEADER 
    hguard = uppercase(fname)*"_MPC_H"
    @printf(fh, "#ifndef %s\n",   hguard);
    @printf(fh, "#define %s\n\n", hguard);

    @printf(fh, "#define N_THETA %d\n",nth);
    @printf(fh, "#define N_STATE %d\n",mpc.nx);
    @printf(fh, "#define N_REFERENCE %d\n",mpc.nr);
    @printf(fh, "#define N_IN_DISTURBANCE %d\n",mpc.nw);
    @printf(fh, "#define N_OUT_DISTURBANCE %d\n",mpc.nd);
    @printf(fh, "#define N_CONTROL_PREV %d\n",mpc.nuprev);

    @printf(fh, "#define N_CONTROL %d\n\n",mpc.nu);

    @printf(fh, "extern c_float mpc_parameter[%d];\n", nth);

    @printf(fh, "extern c_float Dth[%d];\n", nth*m);
    @printf(fh, "extern c_float du[%d];\n", m);
    @printf(fh, "extern c_float dl[%d];\n\n", m);

    @printf(fh, "extern c_float Xth[%d];\n\n", mpc.nu*nth);
    @printf(fh, "extern c_float uscaling[%d];\n\n", mpc.nu);


    # SRC 
    write_float_array(fsrc,zeros(nth),"mpc_parameter");
    write_float_array(fsrc,mpLDP.Dth[:],"Dth");
    write_float_array(fsrc,mpLDP.du[:],"du");
    write_float_array(fsrc,mpLDP.dl[:],"dl");
    write_float_array(fsrc,mpLDP.Xth[:],"Xth");
    write_float_array(fsrc,mpLDP.uscaling[:],"uscaling");

    fmpc_h = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_update_qp.h"), "r");
    write(fh, read(fmpc_h))

    @printf(fsrc, "#include \"%s.h\"\n",fname);
    fmpc_src = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_update_qp.c"), "r");
    write(fsrc, read(fmpc_src))

    @printf(fh, "#endif // ifndef %s\n", hguard);

    close(fh)
    close(fsrc)
    close(fmpc_h)
    close(fmpc_src)
end

function write_float_array(f,a::Vector{<:Real},name::String)
    N = length(a)
    @printf(f, "c_float %s[%d] = {\n", name, N);
    for el in a 
        @printf(f, "(c_float)%.20f,\n", el);
    end
    @printf(f, "};\n");
end


function qp2ldp(mpQP,n_control;normalize=true)
    n = size(mpQP.H,1)
    nb = length(mpQP.bu)-size(mpQP.A,1)
    R = cholesky((mpQP.H+mpQP.H')/2)
    Mext = [Matrix{Float64}(I(n)[1:nb,:]); mpQP.A]/R.U
    Vth = (R.L)\mpQP.f_theta
    v = R.L\mpQP.f#  usually zero since mpQP.f = 0 for MPC
    Dth = mpQP.W + Mext*Vth
    Δd = Mext*v
    du = mpQP.bu[:]+Δd 
    dl = mpQP.bl[:]+Δd

    if(normalize)
        # Normalize
        norm_factors = zeros(size(Mext,1))
        for i in 1:size(Mext,1)
            norm_factor = norm(Mext[i,:],2)
            norm_factors[i] = norm_factor
            Mext[i,:]./=norm_factor
            Dth[i,:]./=norm_factor
            du[i]/= norm_factor
            dl[i]/= norm_factor
        end
        uscaling = norm_factors[1:n_control]
    else
        uscaling = ones(n_control)
    end
    Xth = -mpQP.H\mpQP.f_theta;
    Xth = Xth[1:n_control,:];

    # col major => row major
    Dth = Dth'[:,:]
    Xth = Xth'[:,:]

    return (M=Mext[n+1:end,:], Dth=Dth, du=du, dl=dl, Xth=Xth, uscaling=uscaling)
end
