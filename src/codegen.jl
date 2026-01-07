function codegen(mpc::MPC;fname="mpc_workspace", dir="codegen", opt_settings=nothing, src=true, float_type="double",warm_start=false)
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

    # Append MPC-specific data/functions
    render_mpc_workspace(mpc;fname,dir,float_type, fmode="a",warm_start)

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
    hguard = "MPC_COMPUTE_CONTROL_H"
    @printf(fh, "#ifndef %s\n",   hguard);
    @printf(fh, "#define %s\n\n", hguard);

    write(fh, "typedef $float_type c_float;\n")
    @printf(fh, "#define N_STATE %d\n",mpc.model.nx);
    @printf(fh, "#define N_REFERENCE %d\n",mpc.nr);
    @printf(fh, "#define N_DISTURBANCE %d\n",mpc.model.nd);
    @printf(fh, "#define N_CONTROL_PREV %d\n",mpc.nuprev);
    @printf(fh, "#define N_LINEAR_COST %d\n",mpc.nl);

    # Generate move blocking data for linear cost averaging
    if mpc.nl > 0 && !isempty(mpc.move_blocks)
        if any(x -> x != mpcmove_blocks[1], mpc.move_blocks)
            throw(ArgumentError("Codegeneration not supported for parametric linear cost + varying move blocks."))
        end
        @printf(fh, "#define N_MOVE_BLOCKS %d\n", length(mpc.move_blocks[1]))
        @printf(fh, "#define N_PREDICTION_HORIZON %d\n", mpc.Np)
        @printf(fh, "extern int move_blocks[%d];\n", length(mpc.move_blocks[1]))
    end

    @printf(fh, "extern c_float mpc_parameter[%d];\n", nth);

    if mpc.nl > 0
        @printf(fh, "int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance, c_float* linear_cost);\n");
    else
        @printf(fh, "int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance);\n");
    end

    # SOURCE
    fsrc = open(joinpath(dir,"mpc_compute_control.c"), "w")
    write(fsrc, "#include \"mpc_compute_control.h\"\n")
    write(fsrc, "#include \"$fname.h\"\n")
    if(mpc.settings.reference_condensation)
        @printf(fsrc, "#define N_PREVIEW_HORIZON %d\n",mpc.Np)
        write_float_array(fsrc,mpc.traj2setpoint[:],"traj2setpoint");
    end

    # Generate move_blocks array for linear cost averaging
    if mpc.nl > 0 && !isempty(mpc.move_blocks)
        @printf(fsrc, "#define N_CONTROL %d\n", mpc.model.nu)
        write_int_array(fsrc, mpc.move_blocks[1], "move_blocks")
    end

    # Update parameter
    fmpc_para = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_update_parameter.c"), "r");
    write(fsrc, read(fmpc_para))
    close(fmpc_para)

    # Compute control
    if mpc.nl > 0
        write(fsrc, """
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance, c_float* linear_cost){
    c_float mpc_parameter[$nth];
    // update parameter
    mpc_update_parameter(mpc_parameter,control,state,reference,disturbance,linear_cost);

    // Get the solution at the parameter
    $(fname)_evaluate(mpc_parameter,control);

    return 1;
}
          """)
    else
        write(fsrc, """
int mpc_compute_control(c_float* control, c_float* state, c_float* reference, c_float* disturbance){
    c_float mpc_parameter[$nth];
    // update parameter
    mpc_update_parameter(mpc_parameter,control,state,reference,disturbance);

    // Get the solution at the parameter
    $(fname)_evaluate(mpc_parameter,control);

    return 1;
}
          """)
    end

    if !isnothing(mpc.state_observer)
        @printf(fh, "#define N_CONTROL %d\n",mpc.model.nu);
        codegen(mpc.state_observer,fh,fsrc)
    end

    @printf(fh, "#endif // ifndef %s\n", hguard);
    close(fh)
    close(fsrc)

    @info "Generated code for EMPC controller" dir fname
end

function render_mpc_workspace(mpc;fname="mpc_workspace",dir="",fmode="w", float_type="double", warm_start=false)
    mpLDP = qp2ldp(mpc.mpQP,mpc.model.nu) 
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
    @printf(fh, "#define N_DISTURBANCE %d\n",mpc.model.nd);
    @printf(fh, "#define N_CONTROL_PREV %d\n",mpc.nuprev);
    @printf(fh, "#define N_LINEAR_COST %d\n",mpc.nl);

    @printf(fh, "#define N_CONTROL %d\n\n",mpc.model.nu);

    if warm_start
        @printf(fh, "#define DAQP_WARMSTART %d\n\n")
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

    # Generate move blocking data for linear cost averaging
    if mpc.nl > 0 && !isempty(mpc.move_blocks)
        if any(x -> x != mpc.move_blocks[1], mpc.move_blocks)
            throw(ArgumentError("Codegeneration not supported for parametric linear cost + varying move blocks."))
        end
        @printf(fh, "#define N_MOVE_BLOCKS %d\n", length(mpc.move_blocks[1]))
        @printf(fh, "#define N_PREDICTION_HORIZON %d\n", mpc.Np)
        @printf(fh, "extern int move_blocks[%d];\n", length(mpc.move_blocks[1]))
        write_int_array(fsrc, mpc.move_blocks[1], "move_blocks")
    end

    fmpc_h = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_update_qp.h"), "r");
    write(fh, read(fmpc_h))
    close(fmpc_h)

    @printf(fsrc, "#include \"%s.h\"\n",fname);
    fmpc_para = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_update_parameter.c"), "r");
    write(fsrc, read(fmpc_para))
    close(fmpc_para)
    fmpc_src = open(joinpath(dirname(pathof(LinearMPC)),"../codegen/mpc_update_qp.c"), "r");
    write(fsrc, read(fmpc_src))
    close(fmpc_src)

    if !isnothing(mpc.state_observer)
        codegen(mpc.state_observer,fh,fsrc)
    end

    @printf(fh, "#endif // ifndef %s\n", hguard);

    close(fh)
    close(fsrc)
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
            if(norm_factor>0)
                Mext[i,:]./=norm_factor
                Dth[i,:]./=norm_factor
                du[i]/= norm_factor
                dl[i]/= norm_factor
            end
        end
        uscaling = norm_factors[1:n_control]
    else
        uscaling = ones(n_control)
    end
    Uth_offset = -mpQP.H\mpQP.f_theta;
    Uth_offset = Uth_offset[1:n_control,:];

    u_offset = -mpQP.H\mpQP.f;
    u_offset = u_offset[1:n_control]
    # col major => row major
    Dth = Dth'[:,:]
    Uth_offset = Uth_offset'[:,:]

    return (M=Mext[n+1:end,:], Dth=Dth, du=du, dl=dl, 
            Uth_offset=Uth_offset, u_offset = u_offset, uscaling=uscaling)
end
