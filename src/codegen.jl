function control_codegen_signature(mpc, name="mpc_compute_control")
    args = ["c_float* control", "c_float* state", "c_float* reference", "c_float* disturbance"]
    mpc.np > 0 && push!(args, "c_float* affine_parameter")
    return "int $name(" * join(args, ", ") * ");\n"
end

function control_codegen_args(mpc)
    args = ["control", "state", "reference", "disturbance"]
    mpc.np > 0 && push!(args, "affine_parameter")
    return join(args, ", ")
end

function control_codegen_definition(mpc, name="mpc_compute_control")
    args = ["c_float* control", "c_float* state", "c_float* reference", "c_float* disturbance"]
    mpc.np > 0 && push!(args, "c_float* affine_parameter")
    return "int $name(" * join(args, ", ") * ")"
end

"""
    codegen(mpc; fname="mpc_workspace", dir="codegen", opt_settings=nothing, src=true,
            float_type="double", warm_start=false, bnb_warm_start=mpc.settings.bnb_warm_start)

Generates C code for `mpc` in the directory `dir`. The control is computed by the C function
`mpc_compute_control`.

* `opt_settings`: settings of DAQP in the generated code (a `Dict`, see `DAQP.settings`)
* `src`: copy the source files of DAQP into `dir`
* `float_type`: `"double"` or `"float"`
* `warm_start`: start the solve of a QP from the working set of the previous call (problems without
  binary controls)
* `bnb_warm_start`: warm start of the branch and bound of problems with binary controls from the
  solution of the previous call, as in [`solve`](@ref LinearMPC.solve) with the setting
  `bnb_warm_start`. The generated header then defines `DAQP_BNB_WARMSTART`, the variable
  `bnb_candidate_used` tells whether the latest call of `mpc_compute_control` returned the
  candidate, and `mpc_reset_bnb_warm_start()` discards the stored solution. It has no effect
  without binary controls or with a prestabilizing feedback.
"""
function codegen(mpc::MPC;fname="mpc_workspace", dir="codegen", opt_settings=nothing, src=true, float_type="double",warm_start=false,
                 bnb_warm_start=mpc.settings.bnb_warm_start)
    length(dir)==0 && (dir="codegen")
    dir[end] != '/' && (dir*="/") ## Make sure it is a correct directory path
    ## Generate mpQP
    setup!(mpc)
    # Generate QP workspace
    d = mpc.opt_model
    if(!isnothing(opt_settings))
        DAQP.settings(d,opt_settings)
    end
    DAQP.codegen(d;fname,dir,src)

    if( float_type == "float" || float_type == "single")
        # Append #define DAQP_SINGLE_PRECISION at the top of types
        mv(joinpath(dir,"types.h"),joinpath(dir,"types_old.h"))
        fold = open(joinpath(dir,"types_old.h"),"r")
        s = read(fold, String)
        fnew = open(joinpath(dir,"types.h"),"w")
        write(fnew, "#ifndef DAQP_SINGLE_PRECISION\n # define DAQP_SINGLE_PRECISION\n#endif \n"*s);
        close(fold)
        close(fnew)
        rm(joinpath(dir,"types_old.h"))
    end

    # Append preamble at the start of H file
    mv(joinpath(dir,fname*".h"),joinpath(dir,"mpc_old.h"))
    fold = open(joinpath(dir,"mpc_old.h"),"r")
    fpre = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/preamble.c"), "r");
    s = read(fold, String)
    pre = read(fpre, String)
    fnew = open(joinpath(dir,fname*".h"),"w")
    write(fnew, pre*s)
    close(fpre)
    close(fold)
    close(fnew)
    rm(joinpath(dir,"mpc_old.h"))

    # Append MPC-specific data/functions
    render_mpc_workspace(mpc;fname,dir,float_type, fmode="a",warm_start,bnb_warm_start)

    @info "Generated code for MPC controller" dir fname
end

function codegen(mpc::ExplicitMPC;fname="empc", dir="codegen", opt_settings=nothing, src=true,float_type="double")
    length(dir)==0 && (dir="codegen")
    dir[end] != '/' && (dir*="/") ## Make sure it is a correct directory path
    # Generate code for explicit solution
    ParametricDAQP.codegen(mpc.solution;dir,fname,float_type)

    # Generate code for MPC
    nth = sum(get_parameter_dims(mpc)) 

    # HEADER
    fh = open(joinpath(dir,"mpc_compute_control.h"), "w")

    fpre = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/preamble.c"), "r");
    write(fh, read(fpre))
    close(fpre)

    hguard = "MPC_COMPUTE_CONTROL_H"
    @printf(fh, "#ifndef %s\n",   hguard);
    @printf(fh, "#define %s\n\n", hguard);

    write(fh, "typedef $float_type c_float;\n")
    @printf(fh, "#define N_STATE %d\n",mpc.model.nx);
    @printf(fh, "#define N_REFERENCE %d\n",mpc.nr);
    @printf(fh, "#define N_DISTURBANCE_BASE %d\n",mpc.model.nd);
    @printf(fh, "#define N_DISTURBANCE %d\n",mpc.nd);
    mpc.settings.disturbance_preview && @printf(fh, "#define N_DISTURBANCE_PREVIEW_HORIZON %d\n",mpc.Np);
    @printf(fh, "#define N_CONTROL_PREV %d\n",mpc.nuprev);
    @printf(fh, "#define N_AFFINE_PARAMETER %d\n",mpc.np);
    @printf(fh, "#define N_AFFINE_PARAMETER_BASE %d\n", get_affine_parameter_base_dim(mpc));
    mpc.settings.parameter_preview && @printf(fh, "#define N_AFFINE_PARAMETER_HORIZON %d\n", mpc.Np);

    @printf(fh, "extern c_float mpc_parameter[%d];\n", nth);

    write(fh, control_codegen_signature(mpc))

    # SOURCE
    fsrc = open(joinpath(dir,"mpc_compute_control.c"), "w")

    write(fsrc, "#include \"mpc_compute_control.h\"\n")
    write(fsrc, "#include \"$fname.h\"\n")
    if(mpc.settings.reference_condensation)
        @printf(fsrc, "#define N_PREVIEW_HORIZON %d\n",mpc.Np)
        write_float_array(fsrc,mpc.traj2setpoint[:],"traj2setpoint");
    end

    # Update parameter
    fmpc_para = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_update_parameter.c"), "r");
    write(fsrc, read(fmpc_para))
    close(fmpc_para)

    # Compute control
    write(fsrc, """
$(control_codegen_definition(mpc)){
    c_float mpc_parameter[$nth];
    // update parameter
    mpc_update_parameter(mpc_parameter,$(control_codegen_args(mpc)));

    // Get the solution at the parameter
    $(fname)_evaluate(mpc_parameter,control);

    return 1;
}
          """)

    if !isnothing(mpc.state_observer)
        mpc.settings.disturbance_preview && throw(ArgumentError("Codegeneration not supported for disturbance preview with a state observer."))
        @printf(fh, "#define N_CONTROL %d\n",mpc.model.nu);
        codegen(mpc.state_observer,mpc,fh,fsrc)
    end

    @printf(fh, "#endif // ifndef %s\n", hguard);
    close(fh)
    close(fsrc)

    @info "Generated code for EMPC controller" dir fname
end

function render_mpc_workspace(mpc;fname="mpc_workspace",dir="",fmode="w", float_type="double", warm_start=false, bnb_warm_start=false)
    mpLDP = qp2ldp(mpc.mpQP,mpc.model.nu) 
    # Without binary decision variables (also with a prestabilizing feedback, for which the control
    # bounds are general constraints), the code is the same as without the warm start
    bnb_warm_start = bnb_warm_start && !isempty(mpc.bnb.binary_ids)
    mpLDP.Uth_offset[1:mpc.model.nx,:] -= mpc.K' #Account for prestabilizing feedback
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
    @printf(fh, "#define N_STATE %d\n",mpc.model.nx);
    @printf(fh, "#define N_REFERENCE %d\n",mpc.nr);
    @printf(fh, "#define N_DISTURBANCE_BASE %d\n",mpc.model.nd);
    @printf(fh, "#define N_DISTURBANCE %d\n",mpc.nd);
    mpc.settings.disturbance_preview && @printf(fh, "#define N_DISTURBANCE_PREVIEW_HORIZON %d\n",mpc.Np);
    @printf(fh, "#define N_CONTROL_PREV %d\n",mpc.nuprev);
    @printf(fh, "#define N_AFFINE_PARAMETER %d\n",mpc.np);
    @printf(fh, "#define N_AFFINE_PARAMETER_BASE %d\n", get_affine_parameter_base_dim(mpc));
    mpc.settings.parameter_preview && @printf(fh, "#define N_AFFINE_PARAMETER_HORIZON %d\n", mpc.Np);

    @printf(fh, "#define N_CONTROL %d\n\n",mpc.model.nu);

    if warm_start
        @printf(fh, "#define DAQP_WARMSTART %d\n\n", 1)
    end

    if bnb_warm_start
        @printf(fh, "#define DAQP_BNB_WARMSTART %d\n", 1)
        @printf(fh, "#define N_BNB_BINARY %d\n\n", length(mpc.bnb.binary_ids))
    end

    @printf(fh, "extern c_float mpc_parameter[%d];\n", nth);

    @printf(fh, "extern c_float Dth[%d];\n", nth*m);
    @printf(fh, "extern c_float du[%d];\n", m);
    @printf(fh, "extern c_float dl[%d];\n\n", m);

    @printf(fh, "extern c_float Uth_offset[%d];\n\n", mpc.model.nu*nth);
    @printf(fh, "extern c_float u_offset[%d];\n\n", mpc.model.nu);
    @printf(fh, "extern c_float uscaling[%d];\n\n", mpc.model.nu);


    # SRC 
    write_float_array(fsrc,zeros(nth),"mpc_parameter");
    write_float_array(fsrc,mpLDP.Dth[:],"Dth");
    write_float_array(fsrc,mpLDP.du[:],"du");
    write_float_array(fsrc,mpLDP.dl[:],"dl");
    write_float_array(fsrc,mpLDP.Uth_offset[:],"Uth_offset");
    write_float_array(fsrc,mpLDP.u_offset[:],"u_offset");
    write_float_array(fsrc,mpLDP.uscaling[:],"uscaling");

    if(mpc.settings.reference_condensation)
        @printf(fsrc, "#define N_PREVIEW_HORIZON %d\n",mpc.Np)
        @printf(fsrc, "extern c_float traj2setpoint[%d];\n", length(mpc.traj2setpoint));
        write_float_array(fsrc,mpc.traj2setpoint[:],"traj2setpoint");
    end

    bnb_warm_start && render_bnb_warm_start_data(mpc,mpLDP,fh,fsrc)

    fmpc_h = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_update_qp.h"), "r");
    write(fh, read(fmpc_h))
    close(fmpc_h)

    @printf(fsrc, "#include \"%s.h\"\n",fname);
    fmpc_para = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_update_parameter.c"), "r");
    write(fsrc, read(fmpc_para))
    close(fmpc_para)
    fmpc_src = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_update_qp.c"), "r");
    qp_src = read(fmpc_src, String)
    close(fmpc_src)
    if bnb_warm_start
        # mpc_compute_control calls the branch and bound with the warm start (mpc_bnb_warm_start.c)
        # instead of daqp_bnb, so that the code is unchanged without the warm start
        bnb_call = "daqp_bnb(&daqp_work)"
        count(bnb_call, qp_src) == 1 || error("mpc_update_qp.c is expected to call $bnb_call once")
        qp_src = replace(qp_src, bnb_call => "mpc_bnb_warm_start()")
    end
    write(fsrc, qp_src)
    if bnb_warm_start
        fbnb_h = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_bnb_warm_start.h"), "r");
        write(fh, read(fbnb_h))
        close(fbnb_h)
        fbnb_src = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_bnb_warm_start.c"), "r");
        write(fsrc, read(fbnb_src))
        close(fbnb_src)
    end

    if !isnothing(mpc.state_observer)
        mpc.settings.disturbance_preview && throw(ArgumentError("Codegeneration not supported for disturbance preview with a state observer."))
        codegen(mpc.state_observer,mpc,fh,fsrc)
    end

    @printf(fh, "#endif // ifndef %s\n", hguard);

    close(fh)
    close(fsrc)
end

# Data of the warm start of the branch and bound in the generated code (see mpc_bnb_warm_start.c):
# for each binary decision variable, its index and bounds, and the index, the upper bound and the
# normalization of the LDP row of the decision variable that holds its control one step later
function render_bnb_warm_start_data(mpc,mpLDP,fh,fsrc)
    mpQP = mpc.mpQP
    binary_ids,shift_ids = mpc.bnb.binary_ids,mpc.bnb.shift_ids
    # The stored solution and the bounds of the candidate are in QP variables, which requires
    # bounds that do not depend on the parameter (as for the control bounds without a
    # prestabilizing feedback)
    if !iszero(mpQP.W[union(binary_ids,shift_ids),:])
        throw(ArgumentError("The warm start of the branch and bound in the generated code requires bounds of the binary controls that do not depend on the parameter"))
    end
    nb = length(binary_ids)
    for (type,name) in (("int","bnb_binary_ids"),("int","bnb_shift_ids"),
                        ("c_float","bnb_lower"),("c_float","bnb_upper"),
                        ("c_float","bnb_shift_upper"),("c_float","bnb_shift_scaling"))
        @printf(fh, "extern %s %s[%d];\n", type, name, nb);
    end
    @printf(fh, "\n");
    write_int_array(fsrc,binary_ids.-1,"bnb_binary_ids");
    write_int_array(fsrc,shift_ids.-1,"bnb_shift_ids");
    write_float_array(fsrc,mpQP.bl[binary_ids],"bnb_lower");
    write_float_array(fsrc,mpQP.bu[binary_ids],"bnb_upper");
    write_float_array(fsrc,mpQP.bu[shift_ids],"bnb_shift_upper");
    write_float_array(fsrc,mpLDP.scaling[shift_ids],"bnb_shift_scaling");
end

function write_float_array(f,a::Vector{<:Real},name::String)
    N = length(a)
    @printf(f, "c_float %s[%d] = {\n", name, N);
    for el in a
        @printf(f, "(c_float)%.20f,\n", el);
    end
    @printf(f, "};\n");
end

function write_int_array(f,a::Vector{<:Integer},name::String)
    N = length(a)
    @printf(f, "int %s[%d] = {\n", name, N);
    for el in a
        @printf(f, "%d,\n", el);
    end
    @printf(f, "};\n");
end


function qp2ldp(mpQP,n_control;normalize=true)
    n = size(mpQP.H,1)
    nb = length(mpQP.bu)-size(mpQP.A,1)
    R = hessian_factor(mpQP)
    Mext = [Matrix{Float64}(I(n)[1:nb,:]); mpQP.A]/R.U
    Vth = (R.L)\mpQP.f_theta
    v = R.L\mpQP.f#  usually zero since mpQP.f = 0 for MPC
    Dth = mpQP.W + Mext*Vth
    Δd = Mext*v
    du = mpQP.bu[:]+Δd 
    dl = mpQP.bl[:]+Δd

    uscaling = ones(n_control)
    norm_factors = ones(size(Mext,1))
    if(normalize)
        # Normalize
        for i in 1:size(Mext,1)
            norm_factor = norm(Mext[i,:],2)
            norm_factors[i] = norm_factor
            if(norm_factor>0)
                Mext[i,:]./=norm_factor
                Dth[i,:]./=norm_factor
                du[i]/= norm_factor
                dl[i]/= norm_factor
            end
        end
        if(nb > 0) # XXX might be wrong if dim(mpc.umax) < dim(u)
            uscaling[1:n_control]  .= norm_factors[1:n_control]
        end
    end
    Uth_offset = -(R\mpQP.f_theta);
    Uth_offset = Uth_offset[1:n_control,:];

    u_offset = -(R\mpQP.f);
    u_offset = u_offset[1:n_control]
    # col major => row major
    Dth = Dth'[:,:]
    Uth_offset = Uth_offset'[:,:]

    return (M=Mext[n+1:end,:], Dth=Dth, du=du, dl=dl, 
            Uth_offset=Uth_offset, u_offset = u_offset, uscaling=uscaling, scaling=norm_factors)
end
