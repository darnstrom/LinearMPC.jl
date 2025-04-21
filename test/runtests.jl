using LinearMPC
using Test
using LinearAlgebra
global templib

@testset verbose = true "LinearMPC.jl" begin
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
        src = [f for f in readdir(srcdir) if last(f,1) == "c" && f != "example.c"]
        testlib = "mpctest."* Base.Libc.Libdl.dlext
        if(!isnothing(Sys.which("gcc")))
            run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib $src`; dir=srcdir))
        end
        theta = [5;5;zeros(5)]
        control = ones(1)
        global templib = joinpath(srcdir,testlib)
        ccall(("mpc_compute_control", templib), Cint, (Ptr{Cdouble}, Ptr{Cdouble}), theta, control)
        @test norm(control.-1.7612519326) < 1e-6
    end

    # TODO: add codegen EMPC test

end
