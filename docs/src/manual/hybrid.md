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

In **LinearMPC.jl**, MLD models can be built with an augmented
input vector

```math
\tilde{u}_k = \begin{bmatrix} u_k \\ \delta_k \\ z_k \end{bmatrix},
```

where the physical inputs `u`, binary auxiliary variables `δ`, and continuous
auxiliary variables `z` are all treated as ordinary controller inputs. The
binary entries are then marked with `set_binary_controls!`, and the helper
functions below add the mixed-integer relations on top of the existing
framework.

In this setup, penalties on auxiliary variables are normally added directly
through the input cost (`R`, `Eu`, `eu`) instead of by extending the output
objective.

The helper functions available in `setup.jl` are:

1. `add_logic_constraint!` for relations involving only binary auxiliary variables.
2. `add_indicator_constraint!` for big-M indicator relations.
3. `add_product_constraint!` for products `z = δ h(x,u)`.
4. `add_ifthenelse_constraint!` and `add_ifthenelse_input_constraint!` for affine branch relations.

### Indicator + product

```julia
using LinearMPC

mpc = LinearMPC.MPC([-0.8;;], [1.0 0.0 1.6];
                    C=[1.0;;], Np=1, Nc=1)
set_input_bounds!(mpc; umin=[-1.0, 0.0, -10.0], umax=[1.0, 1.0, 10.0])
set_binary_controls!(mpc, [2])  # δ
set_objective!(mpc; Q=[0.0], Qf=[0.0], R=[1e-6, 1e-6, 1e-6], eu=[0.0, 0.0, -1.0])

# [δ = 1] ↔ [x >= 0]
add_indicator_constraint!(mpc, 2; Ax=[-1.0;;], Au=zeros(1, 3), m=-10.0, M=10.0)

# z = δ*x
add_product_constraint!(mpc, 3, 2; Ax=[1.0;;], Au=zeros(1, 3), m=-10.0, M=10.0)
```

### If-then-else branches

```julia
using LinearMPC, LinearAlgebra

alpha = pi / 3
c, s = cos(alpha), sin(alpha)

mpc = LinearMPC.MPC(zeros(2, 2), [0.0 1.0 0.0; 0.0 0.0 1.0];
                    C=Matrix{Float64}(I, 2, 2), Np=1, Nc=1)
set_input_bounds!(mpc; umin=[0.0, -6.0, -6.0], umax=[1.0, 6.0, 6.0])
set_binary_controls!(mpc, [1])  # sign selector
target = [1.0, 0.0] # choose an affine objective that favors the desired branch output
set_objective!(mpc; Q=[0.0, 0.0], Qf=[0.0, 0.0], R=[1e-6, 1e-6, 1e-6],
               eu=[0.0, -target[1], -target[2]])

# [sign = 1] ↔ [x1 <= 0]
add_indicator_constraint!(mpc, 1; Ax=[1.0 0.0], Au=zeros(1, 3), m=-5.0, M=5.0, ϵ=0.0)

add_ifthenelse_constraint!(mpc, 2, 1;
    Ax_then=0.8 .* [c s],   Au_then=zeros(1, 3),
    Ax_else=0.8 .* [c -s],  Au_else=zeros(1, 3),
    c_then=[0.0], c_else=[0.0],
    m_then=-6.0, M_then=6.0, m_else=-6.0, M_else=6.0,
)
add_ifthenelse_constraint!(mpc, 3, 1;
    Ax_then=0.8 .* [-s c],  Au_then=zeros(1, 3),
    Ax_else=0.8 .* [s c],   Au_else=zeros(1, 3),
    c_then=[0.0], c_else=[0.0],
    m_then=-6.0, M_then=6.0, m_else=-6.0, M_else=6.0,
)
```

