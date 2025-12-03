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

```@example hybrid_mpc
using LinearMPC
# Setup dynamics + horizon
A = [0.0 1 0; 0 0 0; 0 0 0]
B = [0 0 0; 2.5 1 1; -10 0 0]
mpc = LinearMPC.MPC(A,B,0.1;Np=20)

# Setup the binary controls u_2 {0,1} and u3 in {-1,0}
set_binary_controls!(mpc,[2,3])
set_bounds!(mpc;umin=[-Inf;0;-1],umax=[Inf;1;0])

# Setup objetive to prioritize tracking of the attitude x1
set_objective!(mpc;Q=[0.5e4, 1e-2, 1e-1], R = [10,10,10], Rr = 0)

# Enable reference preview
mpc.settings.reference_preview = true
nothing # hide
```
!!! note "Binary control horizon"
    `set_binary_controls!` takes in a third argument which specifies for how many time step the control should be binary (by default, this is equal to the control horizon.) After the binary control horizon, the control is allowed to take continuous values, which can reduce the computational time significantly, with minor effect on the solution.

We simulate the controller with an attitude reference change to 0.5 after 5 time steps with the following code

```@example hybrid_mpc
x0, N = zeros(3), 20; 
rs = [zeros(1,5) 0.5*ones(1,N-5);
      zeros(2,N)];
dynamics = (x,u,d) -> mpc.model.F*x + mpc.model.G*u
sim = LinearMPC.Simulation(dynamics, mpc; x0,N, r=rs)
nothing # hide
```

The result of the simulation can be plotted with 

```@example hybrid_mpc
using Plots
plot(sim)
```

We can see that the attitude is able to reach the setpoint of 0.5. Moreover, we also we that $u_2$ and $u_3$ only take values in $\{0,1\}$ and $\{-1,0\}$, respectively.
