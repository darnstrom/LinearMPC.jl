module LinearMPC
using LinearAlgebra

include("types.jl");
include("utils.jl");
include("setup.jl");

include("mpc2mpqp.jl");
include("mpc_examples.jl");

using Printf
include("codegen.jl");
end
