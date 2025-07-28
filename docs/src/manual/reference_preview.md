# Reference Preview

Reference preview allows the MPC controller to use knowledge of future reference values, enabling better tracking performance for time-varying setpoints.

## Setup

Enable reference preview in the MPC settings:
```julia
mpc.settings.reference_preview = true
setup!(mpc)  # Rebuild the controller with new settings
```

## Usage

### Reference Format
Provide references as a matrix of size `(ny, Np)` where:
- `ny` is the number of outputs
- `Np` is the prediction horizon

```julia
# For a system with 2 outputs and prediction horizon of 5
r_trajectory = [1.0 1.5 2.0 2.0 2.0;   # Reference for output 1
                0.0 0.0 0.5 1.0 1.0]   # Reference for output 2

u = compute_control(mpc, x; r=r_trajectory)
```

### Single Reference
Vector inputs are automatically broadcast across the prediction horizon:
```julia
# This reference will be used for all time steps
u = compute_control(mpc, x; r=[1.0, 0.0])
```

## Example

```@example reference_preview
using LinearMPC

# Create MPC with reference preview
A = [1 1; 0 1]
B = [0; 1] 
C = [1.0 0; 0 1.0]
mpc = LinearMPC.MPC(A, B; C, Np=10, Nc=5)
set_objective!(mpc; Q=[1.0, 1.0], R=[0.1])

# Enable reference preview
mpc.settings.reference_preview = true
setup!(mpc)

# Define reference trajectory (2 outputs × 10 time steps)
r_traj = [zeros(1,5) ones(1,5);      # Step change for output 1
          zeros(1,10)]               # Zero reference for output 2

# Compute control with preview
x = randn(2)
u = compute_control(mpc, x; r=r_traj)
nothing # hide
```

## Simulation

The `Simulation` function automatically handles reference preview when enabled:
```@example reference_preview
# Reference trajectory for simulation
N_sim = 30
r_sim = zeros(2, N_sim)
r_sim[1, 15:end] .= 1.0  # Step reference at time 20

sim = Simulation(mpc; x0=zeros(2), N=N_sim, r=r_sim)
nothing # hide
```

```@example reference_preview
using Plots
plot(sim)
```
Notice how the controller initiates the move from 0 to 1 already before the reference actually changes. If you look really closely, you can see that the output actually initially goes down rather than up before the reference change occurs, this is a phenomenon arising due to the quadratic cost function which makes it worthwhile to allow a small ("small squared" is really small) deviation in order to build up momentum in advance of the large reference change in order to reduce the large error ("large squared" is really large) arising from the step much faster.

## Code Generation

Reference preview is supported in generated C code. The generated `mpc_compute_control` function expects a reference array of size `ny * Np` when preview is enabled.
