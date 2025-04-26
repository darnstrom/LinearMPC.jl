using LinearMPC
mpc,range = LinearMPC.mpc_examples("dcmotor") # Load DCMotor example
# Adjust range 
range.xmin[1:2] = -[0.3;2.0] 
range.xmax[1:2] = [0.3;2.0]
# Compute explicit controller over range
empc = LinearMPC.ExplicitMPC(mpc;range)

LinearMPC.plot_regions(empc,:α,:β)

# Plot slice of critical regions + control law
fix_vals = [0.0,0.0,0.5,0,0]
fig = LinearMPC.plot_regions(empc;fix_vals)
fig = LinearMPC.plot_feedback(empc;fix_vals, u_id = 1)
