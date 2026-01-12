using LinearMPC
mpc,range = LinearMPC.mpc_examples("dcmotor") # Load DCMotor example
# Adjust range 
range.xmin[1:2] = -[0.3;2.0] 
range.xmax[1:2] = [0.3;2.0]
# Compute explicit controller over range
empc = LinearMPC.ExplicitMPC(mpc;range)

using Plots
p1 = plot(empc;parameters=[:x1,:x2],r=[0.5,0.0],title="Critical regions");
p2 = plot(empc;parameters=[:x1,:x2],control=:u1, r=[0.5,0.0], title="Feedback law");
plot(p1,p2, size = (800, 400))
