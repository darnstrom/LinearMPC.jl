using LinearMPC
using Test

@testset "LinearMPC.jl" begin
  mpQP,P_theta = LinearMPC.mpc_examples("invpend"); 
  mpQP,P_theta = LinearMPC.mpc_examples("dcmotor"); 
  mpQP,P_theta = LinearMPC.mpc_examples("aircraft"); 
end
