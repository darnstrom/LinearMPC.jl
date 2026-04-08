using LinearMPC
using Test
using LinearAlgebra
using Statistics
using Random
using RecipesBase
global templib
Random.seed!(1234)
# RecipesBase expects plotting backends to register supported keys; tests exercise recipes directly.
RecipesBase.is_key_supported(::Symbol) = true

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
        settings!(mpc,reference_tracking=false)
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
        mpc,range = LinearMPC.mpc_examples("satellite");
        LinearMPC.mpc2mpqp(mpc)
        mpc,range = LinearMPC.mpc_examples("ballplate");
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

    @testset "Codegen IMPC warm start" begin
        mpc,range = LinearMPC.mpc_examples("invpend")
        srcdir = tempname()
        LinearMPC.codegen(mpc;dir=srcdir)
        src = [f for f in readdir(srcdir) if last(f,1) == "c"]
        @test !isempty(src)
        if(!isnothing(Sys.which("gcc")))
            # Cold start
            testlib = "mpctest."* Base.Libc.Libdl.dlext
            run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib $src`; dir=srcdir))
            u,x,r,d = zeros(1), [5.0;5;0;0], zeros(2), zeros(0)
            global templib = joinpath(srcdir,testlib)
            Us_cold =  zeros(1,100)
            for i in 1:100 
                ccall(("mpc_compute_control", templib), Cint, (Ptr{Cdouble}, Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}), u,x,r,d)
                Us_cold[:,i] .= u
                x = mpc.model.true_dynamics(x,u,d)
            end
            # Cold start
            testlib = "mpctest_warm."* Base.Libc.Libdl.dlext
            run(Cmd(`gcc -lm -fPIC -O3 -DDAQP_WARMSTART -msse3 -xc -shared -o $testlib $src`; dir=srcdir))
            u,x,r,d = zeros(1), [5.0;5;0;0], zeros(2), zeros(0)
            global templib = joinpath(srcdir,testlib)
            Us_warm =  zeros(1,100)
            for i in 1:100 
                ccall(("mpc_compute_control", templib), Cint, (Ptr{Cdouble}, Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}), u,x,r,d)
                Us_warm[:,i] .= u
                x = mpc.model.true_dynamics(x,u,d)
            end

            @test all(abs.(Us_cold - Us_warm) .< 1e-9)
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

        move_block!(mpc,Int[]) # Set empty
        mpqp = LinearMPC.mpc2mpqp(mpc)
        @test length(mpqp.f) == Np*mpc.model.nu

        move_block!(mpc,[1,1,2,3,3])
        mpqp = LinearMPC.mpc2mpqp(mpc)
        @test length(mpqp.f) == 5*mpc.model.nu

        # Pad
        move_block!(mpc,[1,1])
        mpqp = LinearMPC.mpc2mpqp(mpc)
        @test mpc.move_blocks == [[1,9],[1,9]]

        # Clip
        move_block!(mpc,[2,3,3,6,8,9])
        mpqp = LinearMPC.mpc2mpqp(mpc)
        @test mpc.move_blocks == [[2,3,3,2],[2,3,3,2]]

        move_block!(mpc,2)
        mpqp = LinearMPC.mpc2mpqp(mpc)
        @test mpc.move_blocks == [[2,2,2,2,2],[2,2,2,2,2]]

        move_block!(mpc,3)
        mpqp = LinearMPC.mpc2mpqp(mpc)
        @test mpc.move_blocks == [[3,3,3,1],[3,3,3,1]]

        move_block!(mpc,[[1,2,3],[4,2]])
        mpqp = LinearMPC.mpc2mpqp(mpc)
        @test mpc.move_blocks == [[1,2,7],[4,6]]

        move_block!(mpc,[[1,2,3,15,20],[2]])
        mpqp = LinearMPC.mpc2mpqp(mpc)
        @test mpc.move_blocks == [[1,2,3,4],[10]]
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

    @testset "Reference Preview + Prestabilizing Feedback" begin
        Ac = diagm(1=>ones(2))
        Bc = [0.0; 0.0; 1.0]
        Cc = I(3)
        mpc = LinearMPC.MPC(Ac, Bc, 1.0; C=Cc, Np=10, Nc=10)

        set_objective!(mpc; Q=1e-9*[10000.0, 1, 1e-4], R=[1e-9], Qf=[1e6, 1e6, 1e6])
        set_input_bounds!(mpc; umin=[-1], umax=[1])
        add_constraint!(mpc;Ax=[0.0 1.0 0.0],lb=[-1], ub=[1],soft=true)
        add_constraint!(mpc;Ax=[0.0 0.0 1.0],lb=[-1], ub=[1],soft=false)

        LinearMPC.set_prestabilizing_feedback!(mpc)
        mpc.settings.reference_preview = true
        setup!(mpc)

        rs = zeros(3,100)
        rs[1,1:50] .= 1.0
        rs[1,51:end] .= 0.5
        sim = LinearMPC.Simulation(mpc; x0=zeros(3),r=rs, N=100)
        @test sim.ys[1,30] ≈ 1.0 atol=1e-5
        @test sim.ys[1,end] ≈ 0.5 atol=1e-5
    end

    @testset "Codegen Reference Preview - Full" begin
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

    @testset "Codegen Reference Preview - Condensed" begin
        # Test code generation with reference preview enabled
        A = [1 1; 0 1] 
        B = [0; 1]
        C = [1.0 0; 0 1.0]
        mpc = LinearMPC.MPC(A, B; C, Np=5, Nc=5)
        set_bounds!(mpc; umin=[-2.0], umax=[2.0])
        set_objective!(mpc; Q=[1.0, 1.0], R=[0.1])
        
        # Enable reference preview
        mpc.settings.reference_preview = true
        mpc.settings.reference_condensation = true
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
            
            u,d = zeros(1),zeros(0)
            
            global templib = joinpath(srcdir, testlib)
            ccall(("mpc_compute_control", templib), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), u, x, r_traj, d)
                      
            # Test that Julia and C implementations give same result
            @test norm(u - u_julia) < 1e-10
        end

        # Explicit
        empc = LinearMPC.ExplicitMPC(mpc;range=LinearMPC.ParameterRange(mpc),build_tree=true)
        @test norm(u_julia-compute_control(empc, x; r=r_traj))<1e-10

        # Generate code
        srcdir = tempname()
        LinearMPC.codegen(empc; dir=srcdir)
        src = [f for f in readdir(srcdir) if last(f,1) == "c"]
        @test !isempty(src)
        if(!isnothing(Sys.which("gcc")))
            testlib = "mpctest."* Base.Libc.Libdl.dlext
            run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib $src`; dir=srcdir))
            @test isfile(joinpath(srcdir,testlib))
            u,d = zeros(1),zeros(0)

            global templib = joinpath(srcdir, testlib)
            theta  = zeros(5)
            ccall(("mpc_update_parameter", templib), Cint, (Ptr{Cdouble},Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), theta, u, x, r_traj, d)
            @test norm([x;LinearMPC.condense_reference(mpc,vec(r_traj));u] - theta) < 1e-10
            
            ccall(("mpc_compute_control", templib), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), u, x, r_traj, d)
                      
            @test norm(u - u_julia) < 1e-10 || Sys.isapple() # TODO u=NaN on macOS-latest sometimes 
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

        # Generate C code
        srcdir = tempname()
        LinearMPC.codegen(mpc; dir=srcdir)
        src = [f for f in readdir(srcdir) if last(f,1) == "c"]
        @test "bnb.c" ∈ src

        if(!isnothing(Sys.which("gcc")))
            testlib = "mpctest."* Base.Libc.Libdl.dlext
            run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib $src`; dir=srcdir))

            # Test with reference preview (2 outputs × 5 prediction steps = 10 reference values)
            u,d = zeros(3),zeros(0)
            x = zeros(3)

            global templib = joinpath(srcdir, testlib)
            exitflag = ccall(("mpc_compute_control", templib), Cint,
                             (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), u, x, rs, d)
            # Test that Julia and C implementations give same result
            u_julia = compute_control(mpc, x; r=rs)
            @test norm(u - u_julia) < 1e-5
        end

    end

    @testset "Robust MPC" begin
        using LinearMPC
        F = [1.0 1 ;0 1]
        G = [1;0.5;;]
        mpc = LinearMPC.MPC(F,G;Np=10)
        set_prestabilizing_feedback!(mpc)
        set_bounds!(mpc;umin=-ones(1),umax=ones(1))
        set_output_bounds!(mpc;ymin=-0.15*ones(2),ymax=ones(2),soft=false)
        mpqp_nominal = LinearMPC.mpc2mpqp(mpc)
        dynamics = (x,u,d) -> mpc.model.F*x + mpc.model.G*u
        x0 = [0.9;0.5]
        sim_nominal = LinearMPC.Simulation(dynamics, mpc;x0,N=100,r=[0.0;0])
        @test minimum(sim_nominal.xs[2,:]) < -0.1 
        wmin,wmax = -[1e-2;1e-1],[1e-2;1e-1]
        set_disturbance!(mpc,wmin,wmax)
        mpqp_tightened = LinearMPC.mpc2mpqp(mpc)
        @test sum(mpqp_tightened.bu) < sum(mpqp_nominal.bu)
        @test sum(mpqp_tightened.bl) > sum(mpqp_nominal.bl)
        sim_tight = LinearMPC.Simulation(dynamics, mpc;x0,N=100,r=[0.0;0])
        @test minimum(sim_tight.xs[2,:]) > -0.1
    end

    @testset "Control trajectory" begin
        # Test code generation with reference preview enabled
        A = [1 1; 0 1]
        B = [0; 1]
        C = [1.0 0; 0 1.0]
        mpc = LinearMPC.MPC(A, B; C, Np=5, Nc=5)
        x = [0.5;1]
        u = compute_control(mpc,x)
        utraj = compute_control(mpc,x)
        @test u == utraj[1:1]
    end

    @testset "Observer" begin
        mpc0,th_range = LinearMPC.mpc_examples("invpend",100)
        move_block!(mpc0,[1,1,5,10,10])
        LinearMPC.set_state_observer!(mpc0;Q=1e2*[1e-3,1,1e-3,1],R=[1,0.1])
        empc = ExplicitMPC(mpc0;range=th_range,build_tree=true)

        for mpc in [mpc0,empc] 

            N=2000;
            rs = [zeros(2,20) repeat([10;0],1,N)]
            x = zeros(4)

            xs = zeros(4,N);
            ys = zeros(2,N);
            us = zeros(1,N);
            xhats = zeros(4,N);

            LinearMPC.set_state!(mpc.state_observer,x)
            for k = 1:N
                xs[:,k],ys[:,k] = x,mpc.model.C*x + [0.05;0.005].*randn(2)
                LinearMPC.correct!(mpc.state_observer,ys[:,k])
                xhats[:,k] = mpc.state_observer.x
                u = compute_control(mpc,mpc.state_observer.x;r=rs[:,k])
                us[:,k] .= u
                LinearMPC.predict!(mpc.state_observer,u)
                x = mpc.model.F*x + mpc.model.G*u + [0 0; 0.05 0; 0 0;0 0.005]*randn(2)
            end

            @test all(abs.(xs[1,end-50:end].-10) .< 1.0) # Check if x1 converged to 10

            # Test to generate C code
            srcdir = tempname()
            LinearMPC.codegen(mpc; dir=srcdir)
            src = [f for f in readdir(srcdir) if last(f,1) == "c"]

            if !isnothing(Sys.which("gcc"))
                testlib = "mpctest."* Base.Libc.Libdl.dlext
                run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib $src`; dir=srcdir))

                x = randn(4) 
                u = randn(1) 
                y = zeros(2)

                LinearMPC.set_state!(mpc,x)
                xref1 = copy(LinearMPC.predict!(mpc.state_observer,u))
                xref2 = copy(LinearMPC.correct!(mpc.state_observer,y))

                global templib = joinpath(srcdir, testlib)
                ccall(("mpc_predict_state", templib), Cint,
                      (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), x,u,C_NULL)
                @test norm(x-xref1) < 1e-9
                ccall(("mpc_correct_state", templib), Cint,
                      (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), x,y,C_NULL)
                @test norm(x-xref2) < 1e-9
            end
        end
    end
    @testset "Observer + disturbance" begin
        F,G = [1 1; 0 1],[0;1]
        Gd = [1 0;0 0]
        Dd = [0 1]

        mpc = LinearMPC.MPC(F,G;C=[1 0],Gd,Dd)
        LinearMPC.set_state_observer!(mpc;Q=[1.0, 1],R=[1e-2])
        get_measurement = (x,d) -> [x[1] + d[2] + 0.01*randn()]

        sim = LinearMPC.Simulation(mpc;x0=[1;0],
                                   d=[1;1],N=100,get_measurement)
        @test abs(mean(sim.ys[end-20:end])) < 1e-2

        srcdir = tempname()
        LinearMPC.codegen(mpc,dir=srcdir)
 
        src = [f for f in readdir(srcdir) if last(f,1) == "c"]

        if !isnothing(Sys.which("gcc"))
            testlib = "mpctest."* Base.Libc.Libdl.dlext
            run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib $src`; dir=srcdir))

            x = randn(2)
            u = randn(1)
            y = randn(1)
            d = randn(2)

            LinearMPC.set_state!(mpc,x)
            xref1 = copy(LinearMPC.predict!(mpc.state_observer,u,d))
            xref2 = copy(LinearMPC.correct!(mpc.state_observer,y,d))

            global templib = joinpath(srcdir,testlib)
            ccall(("mpc_predict_state", templib), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), x,u,d)
            @test norm(x-xref1) < 1e-9
            ccall(("mpc_correct_state", templib), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), x,y,d)
            @test norm(x-xref2) < 1e-9
        end
    end
    @testset "x0 uncertainty" begin
        F,G = [1 0.1; 0 1], [0.005;0.1;;] # double integrator with Ts=0.1 
        mpc= LinearMPC.MPC(F,G;Ts=0.1,Np=25,C=[1 0;])
        set_bounds!(mpc;umin=[-0.2],umax=[0.2],ymin=[-0.5],ymax=[0.5]) 
        set_x0_uncertainty!(mpc,0.1*ones(2))
        sim= Simulation(mpc;r=[0.5])
        @test abs(sim.xs[1,end] - 0.4) < 1e-6
    end
    @testset "Constant offset" begin
        F,G = [1 0.1; 0 1], [0.005;0.1;;] # double integrator with Ts=0.1
        mpc= LinearMPC.MPC(F,G;Ts=0.1,Np=25,C=[1 0;], f_offset = [0.1;0.1])
        set_objective!(mpc,R=0,Rr=1,Q=1)
        set_bounds!(mpc;umin=[-2],umax=[2],ymin=[-0.5],ymax=[0.5])
        LinearMPC.set_state_observer!(mpc;Q=1e-3)
        dynamics = (x,u,d) -> mpc.model.F*x + mpc.model.G*u + [0.1;0.1]
        get_measurement = (x,d) -> mpc.state_observer.C*x + 0.01*randn(1)
        sim= Simulation(dynamics,mpc;r=[0.5], get_measurement)
        @test all(abs.(sim.xs[1,end-50:end].-0.5) .< 0.1) # Check if x1 converged to 0.5 
    end
    @testset "Operating points" begin
        f = (x,u,d) -> [x[1] - x[2];x[2]+u[1]-1]
        xo,uo = [0.5;0.5],[0.5]
        Ts  = 0.1
        model = LinearMPC.Model(f,(x,u,d)->x,xo,uo,Ts)

        mpc = LinearMPC.MPC(model;Np=100)
        mpc.settings.reference_tracking = false
        set_objective!(mpc;Q=1,R=1,Rr=0)

        sim = LinearMPC.Simulation(mpc;x0 = [0.1;0], N = 100)
        @test norm(sim.xs[:,end]-xo) < 1e-4
        # Update opearting point
        xo,uo = [1;1],[0]
        set_operating_point!(mpc;xo=xo,uo=uo)
        sim = LinearMPC.Simulation(mpc;x0 = [0.1;0], N = 100)
        @test norm(sim.xs[:,end]-xo) < 1e-4
    end

    @testset "Linear Cost" begin
        # Test basic linear cost functionality
        A = [1 1; 0 1]
        B = [0; 1]
        C = [1.0 0; 0 1.0]
        mpc = LinearMPC.MPC(A, B; C, Np=5, Nc=3)
        set_bounds!(mpc; umin=[0.0], umax=[2.0])
        set_objective!(mpc; Q=[1.0, 1.0], R=[0.1])

        # Test with linear cost disabled (default)
        @test mpc.settings.linear_cost == false
        control_no_lincost = compute_control(mpc, [-1.0, 0.0]; r=[0.0, 0.0])
        @test length(control_no_lincost) == 1

        # Enable linear cost
        mpc.settings.linear_cost = true
        setup!(mpc)

        # Test parameter dimensions with linear cost
        nx, nr, nd, nuprev, nl = LinearMPC.get_parameter_dims(mpc)
        @test nx == 2  # State dimension
        @test nr == 2  # Reference dimension (no preview)
        @test nd == 0  # No disturbance
        @test nuprev == 0  # No rate penalty
        @test nl == 1 * 3  # nu × Nc = 1 × 3

        # Test that linear cost affects control
        control_zero_lincost = compute_control(mpc, [-1.0, 0.0]; r=[0.0, 0.0], l=zeros(1))
        control_pos_lincost = compute_control(mpc, [-1.0, 0.0]; r=[0.0, 0.0], l=[1.0])
        control_neg_lincost = compute_control(mpc, [-1.0, 0.0]; r=[0.0, 0.0], l=[-1.0])

        # Linear cost should affect the control differently for positive vs negative
        @test control_zero_lincost ≈ control_no_lincost
        # The linear cost should decrease on increase the control action depending on the sign
        @test control_pos_lincost[1] < control_zero_lincost[] < control_neg_lincost[1]

        # Test with linear cost trajectory matrix (nu × Nc)
        l_traj = [1.0 1.0 1.0]  # Different cost at each time step
        control_traj = compute_control(mpc, [-1.0, 0.0]; r=[0.0, 0.0], l=l_traj)
        @test control_traj[] ≈ control_pos_lincost[] # A trajectory of ones should give same result as a single 1

        # Test that omitting linear cost gives same as zero linear cost
        control_omit = compute_control(mpc, [-1.0, 0.0]; r=[0.0, 0.0])
        @test control_omit[] ≈ control_zero_lincost[]
    end

    @testset "Linear Cost Simulation" begin
        # Test simulation with linear cost trajectory
        A = [0 -0.37; 0.37 0.74]
        B = [0.37; 0.26]
        C = [1.0 0; 0 1.0]
        mpc = LinearMPC.MPC(A, B; C, Np=5, Nc=3)
        set_bounds!(mpc; umin=[-2.0], umax=[2.0])
        set_objective!(mpc; Q=[1.0, 1.0], R=[0.1])
        set_terminal_cost!(mpc)
        mpc.settings.linear_cost = true
        setup!(mpc)

        N_sim = 20
        r_traj = zeros(2, N_sim)

        l_traj = zeros(1, N_sim)
        l_traj[1, :] .= -0.5

        # Simulation with linear cost
        sim_with_l = LinearMPC.Simulation(mpc; x0=[1.0, 0.0], N=N_sim, r=r_traj, l=l_traj)
        @test size(sim_with_l.xs) == (2, N_sim)
        @test size(sim_with_l.us) == (1, N_sim)

        # Simulation without linear cost
        sim_no_l = LinearMPC.Simulation(mpc; x0=[1.0, 0.0], N=N_sim, r=r_traj)

        # Compute cost of trajectories and verify that the cost of the controller optimizing linear cost is smaller
        cost_l = sum(abs2, sim_with_l.xs) +
            0.1*sum(abs2, sim_with_l.us) +
            dot(sim_with_l.us, l_traj) +
            dot(sim_with_l.xs[:, end], mpc.weights.Qf, sim_with_l.xs[:, end])

        cost_no_l = sum(abs2, sim_no_l.xs) +
            0.1*sum(abs2, sim_no_l.us) +
            dot(sim_no_l.us, l_traj) +
            dot(sim_no_l.xs[:, end], mpc.weights.Qf, sim_no_l.xs[:, end])
        @test cost_l < cost_no_l

        # Since linear cost is negative, we should have a positive steady-state control
        @test sim_with_l.us[end] > 0.1
    end

    @testset "Linear Cost Codegen" begin
        # Test code generation with linear cost and move blocking
        A = [1 1; 0 1]
        B = [0; 1]
        C = [1.0 0; 0 1.0]
        mpc = LinearMPC.MPC(A, B; C, Np=10, Nc=10)
        set_bounds!(mpc; umin=[-2.0], umax=[2.0])
        set_objective!(mpc; Q=[1.0, 1.0], R=[0.1])
        mpc.settings.linear_cost = true

        # Add move blocking: [2,3,3,2] → 4 blocked moves covering 10 steps
        move_block!(mpc, [2, 3, 3, 2])
        setup!(mpc)

        x = [0.5, 0.1]
        r = [0.0, 0.0]
        # Full Np trajectory (nu × Np = 1 × 10) - mean will be computed per block
        l_traj = reshape(collect(range(0.1, 1.0, length=10)), 1, 10)

        u_julia = compute_control(mpc, x; r=r, l=l_traj)

        # If we preblock the cost trajectory by computing averages over move blocks and then repeating those averages over the blocks, we should get the same result
        l_preblocked = [fill(mean(l_traj[1:2]), 2); fill(mean(l_traj[3:5]), 3); fill(mean(l_traj[6:8]), 3); fill(mean(l_traj[9:10]), 2)]'
        u_preblocked = compute_control(mpc, x; r=r, l=l_preblocked)
        @test u_preblocked ≈ u_julia

        # Generate C code
        srcdir = tempname()
        LinearMPC.codegen(mpc; dir=srcdir)
        src = [f for f in readdir(srcdir) if last(f,1) == "c"]
        @test !isempty(src)

        if(!isnothing(Sys.which("gcc")))
            testlib = "mpctest."* Base.Libc.Libdl.dlext
            run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib $src`; dir=srcdir))

            u = zeros(1)
            d = zeros(0)

            global templib = joinpath(srcdir, testlib)
            # C code receives full l_traj and computes mean per block internally
            ccall(("mpc_compute_control", templib), Cint,
                  (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
                  u, x, r, d, l_traj)

            # Test that Julia and C implementations give same result
            @test u ≈ u_julia
        end
    end
    @testset "MPQP preprocessing" begin
        A = [0 1; 10 0]
        B = [0; 1]
        C = [1.0 0; 0 1.0]
        mpc = LinearMPC.MPC(A, B, 0.1)
        set_bounds!(mpc;umin=-1,umax=1)
        add_constraint!(mpc;Au=[-1.0;;],lb=[-0.9],ub=[1.5],ks=1:10)
        add_constraint!(mpc;Au=[1.0;;],lb=[-0.5],ub=[2.0],ks=1:10)
        setup!(mpc)
        @test isempty(mpc.mpQP.A)
        @test all(mpc.mpQP.bu .== 0.9*ones(10))
        @test all(mpc.mpQP.bl .== -0.5*ones(10))
    end

    @testset "Set offset" begin
        mpc = LinearMPC.MPC([0.778800783;;],[1.0;;];C=[0.44239843385;;])
        LinearMPC.set_objective!(mpc; Q=[1.0], R=[0.0], Rr=[0.1])
        LinearMPC.set_offset!(mpc; uo=[10.0],ho=[0.5])
        sim = LinearMPC.Simulation(mpc; x0=[0.0], r=[1.5], N=50)
        @test sim.us[end] ≈ 10.5
        @test sim.ys[end] ≈ 1.5
    end

    @testset "Unconstrained" begin
        mpc = LinearMPC.MPC([0.77880078307;;],[1.0];C=[2.211992169;;],Ts=100)
        LinearMPC.move_block!(mpc, [2, 2, 2, 24])
        LinearMPC.set_objective!(mpc; Q=[1], Rr=[0], R=[0])
        sim = LinearMPC.Simulation(mpc; x0=zeros(1), r=[5], N=20)
        @test sim.ys[end] ≈ 5.0
    end

    @testset "Game-theoretic MPC" begin
        F = [1 0.1; 0 1];
        G = [0 0;1 1];
        mpc = LinearMPC.MPC(F,G;C=[1 0;0 1],Np=10);
        set_objective!(mpc, [1]; Q=[1,0], Rr=1e3);
        set_objective!(mpc, [2]; Q=[0,1], Rr=1e3);
        set_bounds!(mpc;umin=-ones(2),umax=ones(2));
        move_block!(mpc,[1,1,8])
        empc = ExplicitMPC(mpc;build_tree=true)

        sim_imp = Simulation(mpc;x0=10*ones(2), r = [10,0], N=500);
        sim_exp = Simulation(empc;x0=10*ones(2), r = [10,0], N=500);

        @test !issymmetric(mpc.mpQP.H)

        @test all(abs.(sim_imp.us-sim_exp.us) .< 1e-4)

        @test sim_imp.ys[1,end] ≈ 10.0 atol=1e-4
        @test sim_imp.ys[2,end] ≈ 0.0 atol=1e-4

        @test sim_exp.ys[1,end] ≈ 10.0 atol=1e-4
        @test sim_exp.ys[2,end] ≈ 0.0 atol=1e-4
    end

    @testset "Model constructors and linearization helpers" begin
        F = [1.0 0.2; 0.0 1.0]
        G = [0.0; 1.0]
        C = [1.0 0.0]

        model = LinearMPC.Model(F, G; C, Gd=[0.5; -0.5])
        @test size(model.Gd) == (2, 1)
        @test size(model.Dd) == (1, 1)
        @test model.true_dynamics([1.0, 2.0], [3.0], [4.0]) ≈ F * [1.0, 2.0] + 3.0 .* G + model.Gd * [4.0]
        @test model.true_h([1.0, 2.0], [3.0], [4.0]) ≈ C * [1.0, 2.0]

        cmodel = LinearMPC.Model([0.0 1.0; 0.0 0.0], [0.0; 1.0], 0.1; C)
        @test cmodel.F ≈ [1.0 0.1; 0.0 1.0] atol=1e-8
        @test cmodel.G ≈ [0.005; 0.1] atol=1e-8

        f(x, u, d) = [x[1] + 2x[2] + u[1] + 3d[1]; x[2] + 2u[1] - d[1]]
        h(x, u, d) = [x[1] - d[1]]
        x0 = [1.0, -1.0]
        u0 = [0.5]
        d0 = [2.0]
        A, B, Bd, Cn, D, Dd, foff, hoff = LinearMPC.linearize(f, h, x0, u0; d=d0)
        @test A ≈ [1.0 2.0; 0.0 1.0]
        @test B ≈ [1.0; 2.0;;]
        @test Bd ≈ [3.0; -1.0;;]
        @test Cn ≈ [1.0 0.0]
        @test D ≈ zeros(1, 1)
        @test Dd ≈ [-1.0;;]
        @test foff ≈ zeros(2)
        @test hoff ≈ zeros(1)

        nmodel = LinearMPC.Model(f, h, x0, u0; d=d0)
        @test nmodel.F ≈ A
        @test nmodel.G ≈ B
        @test nmodel.Gd ≈ Bd
        @test nmodel.C ≈ Cn

        @test_throws ArgumentError LinearMPC.Model(F, G; C=ones(3, 3))
        @test_throws ArgumentError LinearMPC.Model(f, (x, u, d) -> [x[1] + u[1]], x0, u0; d=d0)
    end

    @testset "Reference and linear-cost formatting helpers" begin
        mpc = LinearMPC.MPC([1.0 1.0; 0.0 1.0], [0.0; 1.0]; C=[1.0 0.0; 0.0 1.0], Np=4, Nc=4)
        set_objective!(mpc; Q=[1.0, 1.0], R=[1.0])

        mpc.settings.reference_preview = true
        @test LinearMPC.format_reference(mpc, [1.0, 2.0]) ≈ repeat([1.0, 2.0], mpc.Np)
        @test LinearMPC.format_reference(mpc, [1.0 2.0 3.0 4.0 5.0; 10.0 20.0 30.0 40.0 50.0]) ≈
              [1.0, 10.0, 2.0, 20.0, 3.0, 30.0, 4.0, 40.0]
        @test LinearMPC.format_reference(mpc, [1.0 2.0; 10.0 20.0]) ≈
              [1.0, 10.0, 2.0, 20.0, 2.0, 20.0, 2.0, 20.0]
        @test_throws ErrorException LinearMPC.format_reference(mpc, [1.0])
        @test_throws ErrorException LinearMPC.format_reference(mpc, ones(1, 2))
        @test_throws ErrorException LinearMPC.format_reference(mpc, 1.0)

        mpc.settings.reference_preview = false
        @test LinearMPC.format_reference(mpc, [7.0 8.0 9.0; 1.0 2.0 3.0]) ≈ [7.0, 1.0]

        lmpc = LinearMPC.MPC([1.0 1.0; 0.0 1.0], [0.0; 1.0]; C=[1.0 0.0], Np=4, Nc=3)
        set_objective!(lmpc; Q=[1.0], R=[0.1])
        lmpc.settings.linear_cost = true
        setup!(lmpc)
        @test LinearMPC.format_linear_cost(lmpc, [2.0]) == [2.0, 2.0, 2.0]
        @test LinearMPC.format_linear_cost(lmpc, [1.0, 2.0, 3.0]) == [1.0, 2.0, 3.0]
        @test LinearMPC.format_linear_cost(lmpc, [4.0 5.0]) == [4.0, 5.0, 5.0]
        @test_throws ErrorException LinearMPC.format_linear_cost(lmpc, ones(2, 2))
        @test_throws ErrorException LinearMPC.format_linear_cost(lmpc, 1.0)

        blocked = LinearMPC.MPC([1.0 1.0; 0.0 1.0], [0.0; 1.0]; C=[1.0 0.0], Np=4, Nc=4)
        set_objective!(blocked; Q=[1.0], R=[0.1])
        blocked.settings.linear_cost = true
        move_block!(blocked, [2, 2])
        setup!(blocked)
        @test LinearMPC.format_linear_cost(blocked, [1.0 3.0 5.0 7.0]) ≈ [2.0, 6.0, 0.0]

        varying = LinearMPC.MPC([1.0 1.0; 0.0 1.0], [0.0 0.0; 1.0 1.0]; C=[1.0 0.0], Np=4, Nc=4)
        set_objective!(varying; Q=[1.0], R=[0.1, 0.1])
        varying.settings.linear_cost = true
        move_block!(varying, [[1, 3], [2, 2]])
        setup!(varying)
        @test_throws ArgumentError LinearMPC.format_linear_cost(varying, [1.0 2.0 3.0 4.0; 4.0 3.0 2.0 1.0])
    end

    @testset "Simulation preview and recipe helpers" begin
        rs = [1.0 2.0 3.0; 10.0 20.0 30.0]
        @test LinearMPC.get_preview(rs, 1, 4) == [2.0 3.0 3.0 3.0; 20.0 30.0 30.0 30.0]
        @test LinearMPC.get_reference_preview(rs, 1, 4) == LinearMPC.get_preview(rs, 1, 4)
        ls = [1.0 2.0 3.0 4.0]
        @test LinearMPC.get_linear_cost_preview(ls, 2, 3) == [2.0 3.0 4.0]

        mpc = LinearMPC.MPC([1.0;;], [1.0;;]; C=[1.0;;], Np=3)
        set_objective!(mpc; Q=[1.0], R=[1.0])
        set_bounds!(mpc; umin=[-1.0], umax=[1.0])
        setup!(mpc)
        sim = LinearMPC.Simulation(mpc; x0=[0.0], r=[0.0], N=3)
        attrs = RecipesBase.KW(:xids => [1])
        series = RecipesBase.apply_recipe(attrs, sim)
        @test length(series) == 6
        @test attrs[:layout] == reshape([(1, 1), (1, 1), (1, 1)], 1, :)
        @test series[1].args == (sim.ts, sim.rs[1, :])
        @test series[end].args == (sim.ts, sim.xs[1, :])
    end

    @testset "Explicit MPC parameter labeling and plotting helpers" begin
        mpc = LinearMPC.MPC([1.0 1.0; 0.0 1.0], [0.0; 1.0]; Gd=[1.0; 0.0], C=[1.0 0.0; 0.0 1.0], Dd=[0.0; 1.0], Np=3)
        set_objective!(mpc; Q=[1.0, 1.0], R=[1.0], Rr=[0.5])
        set_input_bounds!(mpc; umin=[-2.0], umax=[2.0])
        set_labels!(mpc; x=[:x1, :x2], u=[:u1], y=[:y1, :y2], d=[:d1])
        setup!(mpc)

        @test LinearMPC.label2id(mpc, :x2) == (2, "x2")
        @test LinearMPC.label2id(mpc, :y1r) == (3, "y1^r")
        @test LinearMPC.label2id(mpc, :d1) == (5, "d1")
        @test LinearMPC.label2id(mpc, :u1p) == (6, "u1^-")
        @test LinearMPC.label2id(mpc, :unknown) == (nothing, "unknown")
        @test LinearMPC.make_subscript("x12") == "x_12"
        @test LinearMPC.make_subscript("state") == "state"

        empc = ExplicitMPC(mpc; range=LinearMPC.ParameterRange(mpc), build_tree=false)
        free_ids, fix_ids, lx, ly = LinearMPC.get_parameter_plot(empc, :y1r, :x1)
        @test free_ids == [1, 3]
        @test all(i ∉ free_ids for i in fix_ids)
        @test lx == "x_1"
        @test ly == "y_1^r"
        @test_throws ArgumentError LinearMPC.get_parameter_plot(empc, :unknown, :x1)
        @test_throws ArgumentError LinearMPC.get_parameter_plot(empc, :x1, :unknown)

        hidden_zero = string(LinearMPC.fixed_parameters_string(empc, [1, 2], [0.0, 2.0]))
        shown_zero = string(LinearMPC.fixed_parameters_string(empc, [1, 2], [0.0, 2.0]; show_zero=true))
        @test !occursin("0.0", hidden_zero)
        @test occursin("2.0", hidden_zero)
        @test occursin("0.0", shown_zero)

        attrs = RecipesBase.KW(:parameters => [:x1, :y1r], :control => :u1)
        series = RecipesBase.apply_recipe(attrs, empc)
        @test length(series) == 1
        @test attrs[:CR_attr][1] == 1
        @test attrs[:CR_attr][2] == [1, 3]
        @test_throws ArgumentError RecipesBase.apply_recipe(RecipesBase.KW(:parameters => [:x1], :control => :u1), empc)

        cert = LinearMPC.certify(mpc; range=LinearMPC.ParameterRange(mpc))
        cert_attrs = RecipesBase.KW(:parameters => [:x1, :y1r])
        cert_series = RecipesBase.apply_recipe(cert_attrs, cert)
        @test length(cert_series) == 1
        @test cert_attrs[:CR_attr][1] == 0
    end

    @testset "MPQP preprocessing helpers" begin
        duplicate_groups = LinearMPC.find_duplicate_rows([1.0 2.0; 1.0000001 2.0; 3.0 4.0]; digits=5)
        @test duplicate_groups == [[1, 2], [3]]

        cdup = LinearMPC.DenseConstraints(
            [1.0 0.0; 1.0 0.0],
            [10.0, 3.0, 5.0],
            [-10.0, -5.0, -1.0],
            zeros(3, 1),
            BitVector([false, false, false]),
            BitVector([false, false, false]),
            # One simple bound followed by two general constraints with duplicate rows.
            Cint[0, 1, 1],
        )
        dedup = LinearMPC.remove_duplicate(cdup)
        @test size(dedup.A, 1) == 1
        @test dedup.bu == [10.0, 3.0]
        @test dedup.bl == [-10.0, -1.0]

        cfold = LinearMPC.DenseConstraints(
            [2.0 0.0],
            [10.0, 4.0],
            [-10.0, 0.0],
            zeros(2, 0),
            BitVector([false, false]),
            BitVector([false, false]),
            Cint[0, 0],
        )
        folded = LinearMPC.remove_redundant(cfold)
        @test isempty(folded.A)
        @test folded.bu == [2.0]
        @test folded.bl == [0.0]

        cpruned = LinearMPC.DenseConstraints(
            [-2.0 0.0; 0.0 0.0],
            [4.0, 1.0],
            [2.0, -1.0],
            zeros(2, 1),
            BitVector([false, false]),
            BitVector([false, false]),
            Cint[0, 0],
        )
        pruned = LinearMPC.remove_redundant(cpruned)
        @test size(pruned.A) == (1, 2)
        @test pruned.A[1, :] == [1.0, 0.0]
        @test pruned.bu == [-1.0]
        @test pruned.bl == [-2.0]

        mpqp_soft = LinearMPC.MPQP(
            [1.0;;], [0.0], zeros(1, 0), zeros(1, 0),
            [1.0;;], [1.0], [-1.0], zeros(1, 0),
            Cint[LinearMPC.DAQP.SOFT], Cint[0], Cint[],
            false, true, zeros(1), ones(1), -ones(1), zeros(1),
        )
        singlesided_single = LinearMPC.make_singlesided(mpqp_soft; single_soft=true)
        singlesided_multi = LinearMPC.make_singlesided(mpqp_soft; single_soft=false)
        @test size(singlesided_single.H, 1) == 2
        @test size(singlesided_single.A, 2) == 2
        @test size(singlesided_multi.H, 1) == 2
        @test size(singlesided_multi.A, 2) == 2

        mpqp_inf = LinearMPC.MPQP(
            [1.0;;], [0.0], zeros(1, 0), zeros(1, 0),
            [1.0;;], [1e21], [-1.0], zeros(1, 0),
            Cint[0], Cint[0], Cint[],
            false, true, zeros(1), ones(1), -ones(1), zeros(1),
        )
        pruned_inf = LinearMPC.make_singlesided(mpqp_inf)
        @test length(pruned_inf.b) == 1
        @test length(pruned_inf.bounds_table) == 1
    end

    @testset "Range, cost, and constraint helpers" begin
        base = LinearMPC.MPC([1.0;;], [1.0;;]; C=[1.0;;], Np=2)
        set_objective!(base; Q=[1.0], R=[1.0])
        setup!(base)
        pr = LinearMPC.ParameterRange(base)
        @test isempty(pr.umin)
        @test isempty(pr.lmin)

        parametric = LinearMPC.MPC([1.0;;], [1.0;;]; C=[1.0;;], Np=3, Nc=2)
        set_input_bounds!(parametric; umin=[-2.0], umax=[2.0])
        set_objective!(parametric; Q=[1.0], R=[1.0], Rr=[0.5])
        parametric.settings.linear_cost = true
        setup!(parametric)
        pr_param = LinearMPC.ParameterRange(parametric)
        region = LinearMPC.range2region(pr_param)
        @test pr_param.umin == [-2.0]
        @test length(pr_param.lmin) == parametric.nl
        @test region.lb == [pr_param.xmin; pr_param.rmin; pr_param.dmin; pr_param.umin; pr_param.lmin]
        @test region.ub == [pr_param.xmax; pr_param.rmax; pr_param.dmax; pr_param.umax; pr_param.lmax]

        mpc = LinearMPC.MPC([1.0;;], [1.0;;]; C=[1.0;;], Np=2)
        xs = [1.0 2.0]
        us = [0.0 1.0]
        rs = [0.0 1.0]
        @test LinearMPC.evaluate_cost(mpc, xs, us, rs; Q=[2.0;;], R=[3.0;;], Rr=[4.0;;], S=[5.0;;]) ≈ 10.5

        c = LinearMPC.Constraint([1.0;;], [1.0 0.0], zeros(0, 0), zeros(0, 0), zeros(0, 0), zeros(0, 0), [1.0], [-1.0], 1:1, false, false, 0)
        @test LinearMPC.constraint_violation(c, [0.8, 0.0], [0.5]) ≈ 0.3
        @test LinearMPC.constraint_violation(c, [0.8 0.2; 0.0 0.0], [0.5 0.0]) ≈ [0.3, 0.0]
        @test_throws AssertionError LinearMPC.constraint_violation(c, [0.8 0.0], [0.5 0.0 0.1])
    end

    @testset "Setup warnings and early-return branches" begin
        mpc = LinearMPC.MPC([1.0;;], [1.0;;]; C=[1.0;;], Np=2)
        @test_logs (:error, r"# of controls are 1") set_input_bounds!(mpc; umin=[-1.0, -2.0], umax=[1.0, 2.0])

        n_constraints = length(mpc.constraints)
        add_constraint!(mpc)
        add_constraint!(mpc; Ax=[1.0;;], ub=zeros(0), lb=zeros(0))
        @test length(mpc.constraints) == n_constraints

        tracked = LinearMPC.MPC([1.0;;], [1.0;;]; C=[1.0;;], Np=2)
        set_objective!(tracked; Q=[1.0], R=[1.0])
        @test false == (@test_logs (:warn, r"LQR cost not valid for reference tracking problems") LinearMPC.set_terminal_cost!(tracked))

        @test_logs (:warn, r"The setting \"does_not_exist\" does not exist") (:warn, r"The setting \"does_not_exist\" does not exist") settings!(tracked; does_not_exist=true)
        @test_logs (:warn, r"The setting \"still_missing\" does not exist") settings!(tracked, Dict(:still_missing => true))
    end
end
