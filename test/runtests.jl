using LinearMPC
using Test

@testset "LinearMPC.jl" begin
  mpc,range = LinearMPC.mpc_examples("invpend"); 
  mpc,range = LinearMPC.mpc_examples("dcmotor"); 
  mpc,range = LinearMPC.mpc_examples("aircraft"); 
end
