using LinearMPC
using Test

@testset "LinearMPC.jl" begin
  mpc,range = LinearMPC.mpc_examples("invpend");
  LinearMPC.mpc2mpqp(mpc)
  mpc,range = LinearMPC.mpc_examples("dcmotor");
  LinearMPC.mpc2mpqp(mpc)
  mpc,range = LinearMPC.mpc_examples("aircraft");
  LinearMPC.mpc2mpqp(mpc)
  mpc,range = LinearMPC.mpc_examples("nonlin");
  LinearMPC.mpc2mpqp(mpc)
  mpc,range = LinearMPC.mpc_examples("mass", 10,10,nx=2);
  LinearMPC.mpc2mpqp(mpc)
  mpc,range = LinearMPC.mpc_examples("chained",10,10,nx=2);
  LinearMPC.mpc2mpqp(mpc)
  mpc,range = LinearMPC.mpc_examples("invpend_contact",10,5);
  LinearMPC.mpc2mpqp(mpc)
end
