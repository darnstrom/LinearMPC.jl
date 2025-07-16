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

    @testset "Reference Preview" begin
        # Test basic reference preview functionality
        A = [0 1; 10 0] 
        B = [0; 1]
        C = [1.0 0; 0 1.0]
        mpc = LinearMPC.MPC(A, B, 0.1; C, Np=5, Nc=3)
        set_bounds!(mpc; umin=[-20.0], umax=[20.0])
        set_objective!(mpc; Q=[1.0, 1.0], R=[0.1], Rr = [0.1])
        
        # Test with reference preview disabled (default)
        @test mpc.settings.reference_preview == false
        control_standard = compute_control(mpc, [1.0, 0.0]; r=[0.0, 0.0])
        @test length(control_standard) == 1
        
        # Enable reference preview
        mpc.settings.reference_preview = true
        setup!(mpc)  # Need to rebuild after changing settings
        
        # Test with single reference (should be broadcast)
        control_single = compute_control(mpc, [1.0, 0.0]; r=[0.0, 0.0])
        @test length(control_single) == 1
        
        # Test with reference trajectory matrix
        r_traj = [0.0 0.5 1.0 1.0 1.0;  # Output 1 trajectory
                  0.0 0.0 0.0 0.0 0.0]  # Output 2 trajectory
        control_preview = compute_control(mpc, [1.0, 0.0]; r=r_traj)
        @test length(control_preview) == 1
        
        # Test parameter dimensions
        nx, nr, nd, nuprev = LinearMPC.get_parameter_dims(mpc)
        @test nx == 2  # State dimension
        @test nr == 2 * 5  # 2 outputs × 5 prediction horizon
        @test nd == 0  # No disturbance
        @test nuprev == 1  # Control rate penalty enabled (Rr = [0.1])
        
        # Test that reference preview gives different result than standard
        # Use a more dynamic reference trajectory that shows clear benefit
        r_dynamic = [0.0 1.0 2.0 1.0 0.0;   # Changing reference for output 1
                     0.0 0.0 0.5 1.0 0.5]   # Changing reference for output 2
        
        mpc.settings.reference_preview = false
        setup!(mpc)
        control_no_preview = compute_control(mpc, [1.0, 0.0]; r=[0.0, 0.0])  # Use only current reference
        
        mpc.settings.reference_preview = true
        setup!(mpc)
        control_with_preview = compute_control(mpc, [1.0, 0.0]; r=r_dynamic)  # Use future knowledge
        
        # They should be different when reference trajectory is changing
        @test norm(control_no_preview - control_with_preview) > 1e-1
    end

    @testset "Reference Preview Simulation" begin
        # Test simulation with reference preview
        A = [1 1; 0 1] 
        B = [0; 1]
        C = [1.0 0; 0 1.0]
        mpc = LinearMPC.MPC(A, B; C, Np=5, Nc=3)
        set_bounds!(mpc; umin=[-2.0], umax=[2.0], ymin=[-1.0;-0.5],  ymax=[1.0;0.5])
        set_objective!(mpc; Q=[1.0, 1.0], R=[0.1])
        
        # Create reference trajectory
        N_sim = 20
        r_traj = zeros(2, N_sim)
        r_traj[1, :] = [zeros(10); ones(10)]  # Step reference for output 1
        r_traj[2, :] = zeros(N_sim)           # Zero reference for output 2
        
        # Test with reference preview enabled
        mpc.settings.reference_preview = true
        setup!(mpc)
        
        sim_preview = LinearMPC.Simulation(mpc; x0=[1.0, 0.0], N=N_sim, r=r_traj)
        @test size(sim_preview.xs) == (2, N_sim)
        @test size(sim_preview.us) == (1, N_sim)
        @test size(sim_preview.rs) == (2, N_sim)

        # Compute explicit solution
        empc = LinearMPC.ExplicitMPC(mpc;range=LinearMPC.ParameterRange(mpc),build_tree=true)
        sim_explicit_preview = LinearMPC.Simulation(empc; x0=[1.0, 0.0], N=N_sim, r=r_traj)
        @test size(sim_explicit_preview.xs) == (2, N_sim)
        @test size(sim_explicit_preview.us) == (1, N_sim)
        @test size(sim_explicit_preview.rs) == (2, N_sim)
        
        # Test with reference preview disabled for comparison
        mpc.settings.reference_preview = false
        setup!(mpc)
        
        sim_no_preview = LinearMPC.Simulation(mpc; x0=[1.0, 0.0], N=N_sim, r=r_traj)
        @test size(sim_no_preview.xs) == (2, N_sim)
        @test size(sim_no_preview.us) == (1, N_sim)
        
        # Control sequences should be different
        @test norm(sim_preview.us - sim_no_preview.us) > 1e-1

        # Implicit and explicit should be similar
        @test norm(sim_preview.ys - sim_explicit_preview.ys) < 1e-6
        
        # Error should be lower with reference preview
        e_preview = sim_preview.ys - sim_preview.rs
        e_no_preview = sim_no_preview.ys - sim_no_preview.rs
        @test norm(e_preview) / norm(e_no_preview) < 0.9
        @test norm(e_preview[:, end])    < 1e-3 # Test that the reference is approximately met in the end
        @test norm(e_no_preview[:, end]) < 1e-3
    end

    @testset "Reference Preview Error Handling" begin
        # Test error handling for reference preview
        A = [0 1; 10 0] 
        B = [0; 1]
        C = [1.0 0; 0 1.0]
        mpc = LinearMPC.MPC(A, B, 0.1; C, Np=5, Nc=3)
        set_bounds!(mpc; umin=[-2.0], umax=[2.0])
        set_objective!(mpc; Q=[1.0, 1.0], R=[0.1])
        mpc.settings.reference_preview = true
        setup!(mpc)
        
        # Test with wrong reference dimensions
        @test_throws ErrorException compute_control(mpc, [1.0, 0.0]; r=[0.0])  # Wrong vector length
        @test_throws ErrorException compute_control(mpc, [1.0, 0.0]; r=[0.0 0.0 0.0])  # Wrong matrix dimensions
        
        # Test with correct dimensions
        @test_nowarn compute_control(mpc, [1.0, 0.0]; r=[0.0, 0.0])  # Correct vector
        @test_nowarn compute_control(mpc, [1.0, 0.0]; r=[0.0 0.0; 0.0 0.0])  # Correct matrix
    end

    @testset "Codegen Reference Preview" begin
        # Test code generation with reference preview enabled
        A = [1 1; 0 1] 
        B = [0; 1]
        C = [1.0 0; 0 1.0]
        mpc = LinearMPC.MPC(A, B; C, Np=5, Nc=5)
        set_bounds!(mpc; umin=[-2.0], umax=[2.0])
        set_objective!(mpc; Q=[1.0, 1.0], R=[0.1])
        
        # Enable reference preview
        mpc.settings.reference_preview = true
        setup!(mpc)

        r_traj = [0.0 0.5 1.0 1.0 1.0;  # Output 1 trajectory
                  0.0 0.0 0.0 0.0 0.0]  # Output 2 trajectory

        x = [0.0, 0.0]
        
        u_julia = compute_control(mpc, x; r=r_traj)
        
        # Generate C code
        srcdir = tempname()
        LinearMPC.codegen(mpc; dir=srcdir)
        src = [f for f in readdir(srcdir) if last(f,1) == "c"]
        @test !isempty(src)
        
        if(!isnothing(Sys.which("gcc")))
            testlib = "mpctest."* Base.Libc.Libdl.dlext
            run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib $src`; dir=srcdir))
            
            # Test with reference preview (2 outputs × 5 prediction steps = 10 reference values)
            u = zeros(1)
            d = zeros(0)
            
            global templib = joinpath(srcdir, testlib)
            ccall(("mpc_compute_control", templib), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), u, x, r_traj, d)
                      
            # Test that Julia and C implementations give same result
            @test norm(u - u_julia) < 1e-10
        end
    end

    @testset "Evalute running cost" begin
        # Get inverted pendulum MPC
        mpc,_= LinearMPC.mpc_examples("invpend")
        # Run simulation
        x0 = [0.0,0.0,0.0,0.0];
        rs = [zeros(2,20) repeat([10;0],1,780) repeat([9;0],1,10)];
        dynamics = (x,u,d) -> mpc.model.F*x + mpc.model.G*u
        sim = LinearMPC.Simulation(dynamics, mpc; x0, N = 1000, r=rs)
        @test LinearMPC.evaluate_cost(mpc,sim) > 0
    end

    @testset "Hybrid MPC" begin
        mpc,_ = LinearMPC.mpc_examples("satellite",20)
        mpc.settings.reference_preview=true
        x0, N = zeros(3), 20;
        rs = [zeros(1,5) 0.5*ones(1,N-5);
              zeros(2,N)];
        dynamics = (x,u,d) -> mpc.model.F*x + mpc.model.G*u
        sim = LinearMPC.Simulation(dynamics, mpc; x0,N, r=rs)
        @test mpc.mpQP.has_binaries # generated mpQP has binaries
        @test abs(sim.ys[1,end]-0.5) < 1e-3;  # convergence to reference
        # Ensure binary controls are either umin or umax
        for bin_id in mpc.binary_controls
            @test all(isapprox.(sim.us[bin_id,:],mpc.umin[bin_id],atol=1e-5) .||
                      isapprox.(sim.us[bin_id,:],mpc.umax[bin_id],atol=1e-5))
        end
    end
end
