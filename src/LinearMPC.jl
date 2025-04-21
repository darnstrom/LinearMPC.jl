module LinearMPC
using LinearAlgebra
using DAQPBase
const DAQP = DAQPBase

include("types.jl");

using ParametricDAQP
include("explicit.jl");

include("utils.jl");
export compute_control
include("setup.jl");
export set_bounds!,add_constraints!,set_output_bounds!,set_weights!
export set_terminal_cost!,set_prestabilizing_feedback!

include("mpc2mpqp.jl");
include("mpc_examples.jl");


using Printf
include("codegen.jl");

end
