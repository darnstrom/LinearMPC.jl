module LinearMPC
using LinearAlgebra
using DAQPBase
const DAQP = DAQPBase

include("model.jl")
include("types.jl");

using ParametricDAQP
include("explicit.jl");

include("utils.jl");
export compute_control
include("setup.jl");
export setup!
export set_bounds!,add_constraint!,set_input_bounds!, set_output_bounds!
export set_objective!, set_weights!, set_horizon!
export set_terminal_cost!,set_prestabilizing_feedback!
export move_block!,set_labels!

include("mpc2mpqp.jl");
include("mpc_examples.jl");

include("simulation.jl");

using Printf
include("codegen.jl");

using ASCertain
include("certify.jl");

using PolyDAQP
include("invariant.jl");
export invariant_set

end
