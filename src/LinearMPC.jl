module LinearMPC
using LinearAlgebra
using JuMP
using Gurobi
using MatrixEquations

include("types.jl");
include("utils.jl");

include("mpc2mpqp.jl");
include("mpc_examples.jl");

using Printf
include("codegen.jl");
end
