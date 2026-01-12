module LinearMPC
using LinearAlgebra
using Statistics
using DAQPBase
const DAQP = DAQPBase
export DAQP

include("model.jl")
include("types.jl");

using ParametricDAQP
include("explicit.jl");
export ExplicitMPC
export build_tree!

include("utils.jl");
export compute_control
export compute_control_trajectory
include("setup.jl");
export setup!
export set_bounds!,add_constraint!,set_input_bounds!, set_output_bounds!
export set_objective!, set_weights!, set_horizon!
export set_binary_controls!
export set_disturbance!, set_x0_uncertainty!
export set_terminal_cost!,set_prestabilizing_feedback!
export set_state_observer!
export set_operating_point!, set_offset!
export move_block!,set_labels!
export settings!

include("mpc2mpqp.jl");
include("mpc_examples.jl");

include("simulation.jl");
export Simulation

using Printf
include("codegen.jl");

using ASCertain
include("certify.jl");

using PolyDAQP
include("invariant.jl");
export invariant_set

include("robust.jl");

include("observer.jl");
export predict_state!,correct_state!
export set_state!,get_state,update_state!

end
