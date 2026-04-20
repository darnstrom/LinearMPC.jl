# Disturbance Preview

Disturbance preview lets the controller use a known disturbance trajectory over the prediction horizon instead of assuming a constant measured disturbance.

## Setup

Enable disturbance preview in the MPC settings:

```@tab
# julia
mpc.settings.disturbance_preview = true
setup!(mpc)
# python
mpc.settings({"disturbance_preview": True})
mpc.setup()
```

## Usage

Provide disturbances as a matrix of size `(nd, Np)` where:
- `nd` is the number of disturbance channels
- `Np` is the prediction horizon

```@tab
# julia
d_traj = [0.0 0.0 0.5 1.0 1.0]  # nd × Np
u = compute_control(mpc, x; d=d_traj)
# python
import numpy as np

d_traj = np.array([[0.0, 0.0, 0.5, 1.0, 1.0]])
u = mpc.compute_control(x, d=d_traj)
```

Vector inputs are automatically broadcast across the horizon:

```@tab
# julia
u = compute_control(mpc, x; d=[0.5])
# python
u = mpc.compute_control(x, d=[0.5])
```

## Example

```@tabsetup disturbance_preview
# julia
using LinearMPC

A = [1.0 1.0; 0.0 1.0]
B = [0.0; 1.0]
Gd = [0.0; 1.0]
C = [1.0 0.0]

mpc = LinearMPC.MPC(A, B; Gd, C, Np=5, Nc=5)
set_bounds!(mpc; umin=[-0.5], umax=[0.5])
set_objective!(mpc; Q=[10.0], R=[0.1])
mpc.settings.disturbance_preview = true
setup!(mpc)

x = [0.0, 0.0]
d_traj = [0.0 1.0 1.0 1.0 1.0]
u = compute_control(mpc, x; d=d_traj)
# python
import numpy as np
from lmpc import MPC

A = np.array([[1.0, 1.0], [0.0, 1.0]])
B = np.array([[0.0], [1.0]])
Gd = np.array([[0.0], [1.0]])
C = np.array([[1.0, 0.0]])

mpc = MPC(A, B, Gd=Gd, C=C, Np=5, Nc=5)
mpc.set_bounds(umin=[-0.5], umax=[0.5])
mpc.set_objective(Q=[10.0], R=[0.1])
mpc.settings({"disturbance_preview": True})
mpc.setup()

x = np.array([0.0, 0.0])
d_traj = np.array([[0.0, 1.0, 1.0, 1.0, 1.0]])
u = mpc.compute_control(x, d=d_traj)
```

## Simulation

`Simulation` automatically passes a disturbance preview window to the controller when disturbance preview is enabled:

```@tabsetup disturbance_preview
# julia
N_sim = 20
d_sim = [zeros(1, 8) ones(1, 12)]
sim = Simulation(mpc; x0=[0.0, 0.0], N=N_sim, d=d_sim)
# python
from lmpc import Simulation

N_sim = 20
d_sim = np.block([np.zeros((1, 8)), np.ones((1, 12))])
sim = Simulation(mpc, x0=np.array([0.0, 0.0]), N=N_sim, d=d_sim)
```

## Code Generation

Disturbance preview is supported in generated C code. The generated `mpc_compute_control` function expects a disturbance array of size `nd * Np` when preview is enabled.
