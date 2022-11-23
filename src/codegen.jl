using DAQP

function codegen(mpc;fname="mpc_workspace", dir="", opt_settings=nothing)
    # Generate QP workspace
    mpQP = mpc2mpqp(mpc)
    d = DAQP.Model()
    DAQP.setup(d,mpQP.H[:,:],mpQP.f[:],mpQP.A[:,:],mpQP.bu[:],mpQP.bl[:],mpQP.senses[:])
    if(!isnothing(opt_settings))
        DAQP.settings(d,opt_settings)
    end
    # TODO: move this to DAQP.jl
    @assert(d.has_model, "setup the model before code generation")
    exitflag = ccall((:render_daqp_workspace,DAQP.libdaqp),Cvoid,
                     (Ptr{DAQP.Workspace},Cstring,Cstring,), d.work,fname,dir);

    # Append MPC-specific data/functions
    mpLDP = dualize(mpQP,mpc.nu) 
    render_mpc_workspace(mpLDP,mpc.nu;fname,dir, fmode="a")
end

function render_mpc_workspace(mpLDP,n_control;fname="mpc_workspace",dir="",fmode="w")

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
    @printf(fh, "#define N_CONTROL %d\n\n",n_control);
    @printf(fh, "extern c_float Dth[%d];\n", nth*m);
    @printf(fh, "extern c_float du[%d];\n", m);
    @printf(fh, "extern c_float dl[%d];\n\n", m);

    @printf(fh, "extern c_float Xth[%d];\n\n", n_control*nth);


    # SRC 
    write_float_array(fsrc,mpLDP.Dth[:],"Dth");
    write_float_array(fsrc,mpLDP.du[:],"du");
    write_float_array(fsrc,mpLDP.dl[:],"dl");
    write_float_array(fsrc,mpLDP.Xth[:],"Xth");

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

