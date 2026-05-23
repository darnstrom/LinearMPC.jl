# [Hybrid MPC](@id man_hybridmpc)

**LinearMPC.jl** can be used to control hybrid systems, where controls and states might take both continuous and binary values.

Binary controls are either equal to its upper or lower bound. For example, if the lower and upper bounds for the control $i$ is $\underline{u}_i$ and $\overline{u}_i$, making $u_i$ binary means that $u_i \in \{\underline{u}_i, \overline{u}_i\}$ rather than $\underline{u}_i \leq u_i \leq \overline{u}_i$. Controls can be made binary with the function `set_binary_controls!`, which is exemplified below.

More generally, **LinearMPC.jl** allows any constraint to be defined as `binary`. A constraint can be enforced to be binary by settings the setting the optional argument `binary` to  `true` in the functions `add_constraint!`.

## Illustrative example
As an illustrative example, we consider the control of the attitude of a satellite[^Axehill04]. The actuators consist of one reaction wheel and two thrusters. The thruster takes on binary values (they are either 'on' or 'off')

[^Axehill04]: Axehill, Daniel, and Hansson, Anders. "A preprocessing algorithm for MIQP solvers with applications to MPC." _43rd IEEE Conference on Decision and Control (CDC)_ (2004) 

The dynamics of the system is given by

```math
\dot{x} = \begin{bmatrix}
    0 & 1 & 0 \\ 0 & 0 & 0\\ 0 & 0 & 0 
\end{bmatrix} x  
+ \begin{bmatrix}
    0 & 0 & 0 \\
    2.5 & 1 & 1 \\
    -10 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
 u_1 \\
 u_2 \\
 u_3 
\end{bmatrix}
```
where $x_1$ i the attitude of the satellite. 

The thrusters give rise to the binary controls $u_2 \in \{0,1\}$ and $u_3 \in \{-1,0\}$.

An MPC controller that controls the attitude of the satellite can be set up with **LinearMPC.jl** as follows: 

```@tabsetup hybrid_mpc
# julia
using LinearMPC
# Setup dynamics + horizon
A = [0.0 1 0; 0 0 0; 0 0 0]
B = [0 0 0; 2.5 1 1; -10 0 0]
mpc = LinearMPC.MPC(A, B, 0.1; Np=20)

# Setup the binary controls u_2 ∈ {0,1} and u_3 ∈ {-1,0}
set_binary_controls!(mpc, [2, 3])
set_bounds!(mpc; umin=[-Inf; 0; -1], umax=[Inf; 1; 0])

# Setup objective to prioritize tracking of the attitude x1
set_objective!(mpc; Q=[0.5e4, 1e-2, 1e-1], R=[10, 10, 10], Rr=0)

# Enable reference preview
mpc.settings.reference_preview = true
# python
import numpy as np
from lmpc import MPC

# Setup dynamics + horizon
A = np.array([[0.0, 1, 0], [0, 0, 0], [0, 0, 0]])
B = np.array([[0, 0, 0], [2.5, 1, 1], [-10, 0, 0]])
mpc = MPC(A, B, Ts=0.1, Np=20)

# Setup the binary controls u_2 in {0,1} and u_3 in {-1,0}
mpc.set_binary_controls([2, 3])
mpc.set_bounds(umin=[-np.inf, 0, -1], umax=[np.inf, 1, 0])

# Setup objective to prioritize tracking of the attitude x1
mpc.set_objective(Q=[0.5e4, 1e-2, 1e-1], R=[10, 10, 10], Rr=0)

# Enable reference preview
mpc.settings({"reference_preview": True})
```
!!! note "Binary control horizon"
    `set_binary_controls!` takes in a third argument which specifies for how many time step the control should be binary (by default, this is equal to the control horizon.) After the binary control horizon, the control is allowed to take continuous values, which can reduce the computational time significantly, with minor effect on the solution.

We simulate the controller with an attitude reference change to 0.5 after 5 time steps with the following code

```@tabsetup hybrid_mpc
# julia
x0, N = zeros(3), 20
rs = [zeros(1, 5) 0.5*ones(1, N-5);
      zeros(2, N)]
dynamics = (x, u, d) -> mpc.model.F*x + mpc.model.G*u
sim = LinearMPC.Simulation(dynamics, mpc; x0, N, r=rs)
# python
from lmpc import Simulation

x0, N = np.zeros(3), 20
rs = np.block([[np.zeros((1, 5)), 0.5*np.ones((1, N-5))],
               [np.zeros((2, N))]])
sim = Simulation(mpc, x0=x0, N=N, r=rs)
```

The result of the simulation can be plotted with 

```@tabexample hybrid_mpc
# julia
using Plots
plot(sim)
# python
import matplotlib.pyplot as plt
plt.plot(sim.ts, sim.ys.T)
plt.xlabel("Time step")
plt.show()
```

We can see that the attitude is able to reach the setpoint of 0.5. Moreover, we also we that $u_2$ and $u_3$ only take values in $\{0,1\}$ and $\{-1,0\}$, respectively.

## Mixed logical dynamical systems

For MLD systems of the form

```math
\begin{aligned}
x_{k+1} &= A x_k + B_u u_k + B_\delta delta_k + B_z z_k + bb, \\
y_k &= C x_k + D_u u_k + D_\delta \delta_k + D_z z_k + bd, \\
E_\delta \delta_k + E_z z_k &\le E_u u_k + E_x x_k + be,
\end{aligned}
```

**LinearMPC.jl** provides an `MLDModel` constructor. Internally, the auxiliary binary variables `δ` and auxiliary continuous variables `z` are appended to the stage decision vector, so the existing mixed-integer mpQP/DAQP pipeline is reused. The constructor and constraint helper use the block names `Bu`, `Bdelta`, `Bz`, `Du`, `Ddelta`, `Dz`, `Eu`, `Edelta`, `Ez`, `Ex`, and `be`.

```@tab
# julia
using LinearMPC

mld = LinearMPC.MLDModel(
    [-0.8;;],            # A
    [1.0;;],             # Bu
    zeros(1, 1),         # Bdelta
    [1.6;;];             # Bz
    C=[1.0;;],
    zmin=[-10.0],
    zmax=[10.0],
    delta_labels=[:delta1],
    z_labels=[:z1],
)

mpc = LinearMPC.MPC(mld; Np=1, Nc=1)
set_input_bounds!(mpc; umin=[-1.0, 0.0, -10.0], umax=[1.0, 1.0, 10.0])
set_objective!(mpc; Q=[1.0], R=[0.1, 0.0, 0.0])

# [delta1 = 1] ↔ [x >= 0]
add_indicator_constraint!(mpc, 1; Ax=[-1.0;;], m=-10.0, M=10.0)

# z1 = delta1 * x
add_product_constraint!(mpc, 1, 1; Ax=[1.0;;], m=-10.0, M=10.0)

decision = compute_control(mpc, [2.0]; r=[0.0])
# decision = [u1, delta1, z1]
```

The helper functions below are available in addition to the generic `add_constraint!` interface:

1. `add_mld_constraint!` adds inequalities with separate `u`, `δ`, and `z` blocks, or directly from the matrices `Eu`, `Edelta`, `Ez`, `Ex`, and `be`.
2. `add_logic_constraint!` adds linear inequalities over the auxiliary binary variables.
3. `add_indicator_constraint!` encodes relations of the form `[δ = 1] ↔ [h(x,u) ≤ 0]`.
4. `add_product_constraint!` encodes products `z = δ h(x,u)`.
5. `add_ifthenelse_constraint!` encodes affine `if-then-else` relations.

!!! note
    For `MLDModel` controllers, `compute_control` returns the full stage decision vector `[u; \delta; z]`. This makes the returned vector directly compatible with the MLD dynamics used inside the model and simulation routines.
