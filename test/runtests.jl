using LinearMPC
using Test
using LinearAlgebra
global templib

@testset verbose = true "LinearMPC.jl" begin
    @testset "Basic setup" begin
        A = randn(3,3)
        B = randn(3)
        C = [1.0 0 0; 0 1.0 0]
        Bd = randn(3)
        Dd = [1.0 0; 0 1]
        mpc = LinearMPC.MPC(A,B,0.1;C,Bd,Dd, Np = 10, Nc = 5)
        set_weights!(mpc,Q=[1.0;3.0], R = 2, Rr = [1.0;;])
        set_bounds!(mpc,umin = -[0.5], umax = [0.5])
        set_prestabilizing_feedback!(mpc)
        set_output_bounds!(mpc,ymin=[0.0;0.0], ymax=[5.0;1.0])
        setup!(mpc)
        set_horizon!(mpc,5)
        setup!(mpc)
    end
    @testset "MPC examples " begin
        mpc,range = LinearMPC.mpc_examples("invpend");
        LinearMPC.mpc2mpqp(mpc)
        mpc,range = LinearMPC.mpc_examples("dcmotor");
        LinearMPC.mpc2mpqp(mpc)
        mpc,range = LinearMPC.mpc_examples("aircraft");
        LinearMPC.mpc2mpqp(mpc)
        mpc,range = LinearMPC.mpc_examples("nonlin");
        LinearMPC.mpc2mpqp(mpc)
        mpc,range = LinearMPC.mpc_examples("mass", 10,10,params=Dict(:nx=>2));
        LinearMPC.mpc2mpqp(mpc)
        mpc,range = LinearMPC.mpc_examples("chained",10,10,params=Dict(:nx=>2));
        LinearMPC.mpc2mpqp(mpc)
        mpc,range = LinearMPC.mpc_examples("invpend_contact",6,6, params=Dict(:nwalls=>1));
        LinearMPC.mpc2mpqp(mpc)
    end
 
    @testset "Compute control" begin
        mpc,range = LinearMPC.mpc_examples("invpend")
        control = LinearMPC.compute_control(mpc,[5.0;5;0;0])
        @test norm(control.-1.7612519326) < 1e-6
    end


    @testset "Codegen IMPC" begin
        mpc,range = LinearMPC.mpc_examples("invpend")
        srcdir = tempname()
        LinearMPC.codegen(mpc;dir=srcdir)
        src = [f for f in readdir(srcdir) if last(f,1) == "c"]
        @test !isempty(src)
        if(!isnothing(Sys.which("gcc")))
            testlib = "mpctest."* Base.Libc.Libdl.dlext
            run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib $src`; dir=srcdir))
            u,x,r,d = zeros(1), [5.0;5;0;0], zeros(2), zeros(0)
            global templib = joinpath(srcdir,testlib)
            ccall(("mpc_compute_control", templib), Cint, (Ptr{Cdouble}, Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}), u,x,r,d)
            @test norm(u.-1.7612519326) < 1e-6
        end
    end

    @testset "Prestabilizing feedback" begin
        A,B = [0 1; 10 0], [0;1]
        mpc = LinearMPC.MPC(A,B,0.1;Np = 30)
        set_bounds!(mpc;umin = -1, umax = 1);
        unom = compute_control(mpc,zeros(2);r = [1,0])
        mpqp = LinearMPC.mpc2mpqp(mpc)

        cond_nom = cond(mpqp.H)

        set_prestabilizing_feedback!(mpc)

        uprestab = compute_control(mpc,zeros(2);r = [1,0])
        mpqp = LinearMPC.mpc2mpqp(mpc)
        cond_prestab = cond(mpqp.H)

        @test(norm(unom-uprestab) < 1e-10)
        @test(cond_prestab < cond_nom)
    end

    @testset "Move blocking" begin 
        using LinearMPC
        Np = 10
        mpc,_= LinearMPC.mpc_examples("aircraft",Np)
        mpc.settings.explicit_soft=false

        move_block!(mpc,Int[]) # Set empty
        mpqp = LinearMPC.mpc2mpqp(mpc)
        @test length(mpqp.f) == Np*mpc.model.nu

        move_block!(mpc,[1,1,2,3,3])
        mpqp = LinearMPC.mpc2mpqp(mpc)
        @test length(mpqp.f) == 5*mpc.model.nu

        # Pad
        move_block!(mpc,[1,1])
        mpqp = LinearMPC.mpc2mpqp(mpc)
        @test mpc.move_blocks == [1,9]

        # Clip
        move_block!(mpc,[2,3,3,6,8,9])
        mpqp = LinearMPC.mpc2mpqp(mpc)
        @test mpc.move_blocks == [2,3,3,2]

        move_block!(mpc,2)
        mpqp = LinearMPC.mpc2mpqp(mpc)
        @test mpc.move_blocks == [2,2,2,2,2]

        move_block!(mpc,3)
        mpqp = LinearMPC.mpc2mpqp(mpc)
        @test mpc.move_blocks == [3,3,3,1]
    end

    @testset "Explicit MPC" begin
        mpc,range = LinearMPC.mpc_examples("invpend")
        empc = LinearMPC.ExplicitMPC(mpc;range)
        LinearMPC.build_tree!(empc)
        control = LinearMPC.compute_control(empc,[5.0;5;0;0])
        @test norm(control.-1.7612519326) < 1e-6
        # Code generation
        srcdir = tempname()
        LinearMPC.codegen(empc;dir=srcdir,float_type="double")
        src = [f for f in readdir(srcdir) if last(f,1) == "c" && f != "example.c"]
        @test !isempty(src)
        if(!isnothing(Sys.which("gcc")))
            testlib = "mpctest."* Base.Libc.Libdl.dlext
            run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib $src`; dir=srcdir))
            u,x,r,d = zeros(1), [5.0;5;0;0], zeros(2), zeros(0)
            global templib = joinpath(srcdir,testlib)
            ccall(("mpc_compute_control", templib), Cint, (Ptr{Cdouble}, Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}), u,x,r,d)
            @test norm(u.-1.7612519326) < 1e-6
        end
    end

    @testset "Certification" begin
        # Load inverted pendulum example and certify the iteration complexity 
        mpc,range = LinearMPC.mpc_examples("invpend")
        mpc.settings.QP_double_sided = true
        result = LinearMPC.certify(mpc;range)
        @test length(result.partition) > 100
    end

    @testset "Invariant sets" begin
        # BBM17 - Example 10.12
        F = [0.5 0 ; 1 -0.5]
        xmin, xmax = -10*ones(2), 10*ones(2)
        wmin,wmax = -ones(2),ones(2)
        H,h = invariant_set(F,xmin,xmax;wmin,wmax,eps_shrink=0.0)
        @test norm(h-[10.0;10.0;10.0;10.0;8.05;8.05]) < 1e-1

        # BBM17 - Example 10.13
        F = [1.5 0 ; 1 -1.5]
        G = [1.0; 0;;]
        xmin, xmax = -10*ones(2), 10*ones(2)
        wmin,wmax = -0.1*ones(2),0.1*ones(2)
        umin,umax = -[5.0],[5.0]
        H,h = invariant_set(F,xmin,xmax;G,umin,umax,wmin,wmax,eps_shrink=0.0)
        @test norm(h-[3.72;3.72;2.008;2.008]) < 1e-2
    end


end
