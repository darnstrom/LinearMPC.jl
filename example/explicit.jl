using LinearMPC
mpc,range = LinearMPC.mpc_examples("dcmotor") # Load DCMotor example
# Adjust range 
range.xmin[1:2] = -[0.3;2.0] 
range.xmax[1:2] = [0.3;2.0]
# Compute explicit controller over range
empc = LinearMPC.ExplicitMPC(mpc;range)

LinearMPC.plot_regions(empc,:x1,:x2,r=[0.5,0.0])
LinearMPC.plot_feedback(empc,:u1,:x1,:x2,r=[0.5,0.0])
