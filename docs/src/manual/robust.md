# [Robust MPC](@id man_robust)

**LinearMPC.jl** can handle uncertaintiy in the system dynamics and in the current state. 

## Uncertainty in the dynamics
Model predictive control requires a model of the system dynamics. In practice, there is often a mismatch between this model and the true dynamics. This mismatch can be captured with an additive disturbance $w_k$ in the dynamics
```math
x_{k+1} = F x_k + G u_k + w_k.
```
If the mismatch is no accounted for, constraints that are predicted to be satisfied might now be so in practice. We illustrate this with the following example of the control of a double integrator.

```@tabexample robust_mpc
# julia
using LinearMPC, Plots

F, G = [1 0.1; 0 1], [0.005; 0.1;;]  # double integrator with Ts=0.1
true_dynamics = (x, u, d) -> F*x + G*u + 0.01*(rand(2) .- 0.5)

mpc_nominal = LinearMPC.MPC(F, G; Ts=0.1, Np=25, C=[1 0;])
set_bounds!(mpc_nominal; umin=[-0.2], umax=[0.2], ymin=[-0.5], ymax=[0.5])

sim_nominal = Simulation(true_dynamics, mpc_nominal; r=[0.5])
hline([0.5], label="Constraint bound", linestyle=:dash)
plot!(sim_nominal.ys[1,:], xlabel="Time step", ylabel="Position [m]", label="Nominal MPC")
# exec-only
ylims!(0.4, 0.6)
# python
import numpy as np
import matplotlib.pyplot as plt
from lmpc import MPC, Simulation

F = np.array([[1, 0.1], [0, 1]])
G = np.array([[0.005], [0.1]])
true_dynamics = lambda x, u, d: F @ x + G @ u + 0.01 * (np.random.rand(2) - 0.5)

mpc_nominal = MPC(F, G, C=np.array([[1, 0]]))
mpc_nominal.set_bounds(umin=[-0.2], umax=[0.2], ymin=[-0.5], ymax=[0.5])

sim_nominal = Simulation(mpc_nominal, f=true_dynamics, r=[0.5])
plt.axhline(0.5, linestyle="--", label="Constraint bound")
plt.plot(sim_nominal.ys[0, :], label="Nominal MPC")
plt.xlabel("Time step"); plt.ylabel("Position [m]")
plt.legend(); plt.show()
```
As the plot shows, the constraint $x_1 \leq 0.5$ is violated due to the model mismatch. 

To account for model mismatch, the function `set_disturbance!(mpc,wmin,wmax)` defines upper and lower bounds on the process noise, and then tightens the constraints to ensure that the original constraint is satisfied (as long as the actual disturbance $w$ is between `wmin` and `wmax`.
!!! note "Tube MPC"
    Doing this a priori constraint tightening is superior to many commonly used tube MPC approaches.[^Zanon21].

[^Zanon21]: Zanon, Mario, and Gros, Sébastien. "On the similarity between two popular tube MPC formulations." _European Control Conference (ECC)_ (2021) 

```@tabexample robust_mpc
# julia
mpc_robust = LinearMPC.MPC(F, G; Ts=0.1, Np=25, C=[1 0;])
set_prestabilizing_feedback!(mpc_robust)
set_disturbance!(mpc_robust, [-0.005; -0.005], [0.005; 0.005])
set_bounds!(mpc_robust; umin=[-0.2], umax=[0.2], ymin=[-0.5], ymax=[0.5])

sim_robust = Simulation(true_dynamics, mpc_robust; r=[0.5])
plot!(sim_robust.ys[1,:], label="Robust MPC")
# exec-only
ylims!(0.4, 0.6)
# python
mpc_robust = MPC(F, G, C=np.array([[1, 0]]))
mpc_robust.set_prestabilizing_feedback()
mpc_robust.set_disturbance([-0.005, -0.005], [0.005, 0.005])
mpc_robust.set_bounds(umin=[-0.2], umax=[0.2], ymin=[-0.5], ymax=[0.5])

sim_robust = Simulation(mpc_robust, f=true_dynamics, r=[0.5])
plt.plot(sim_robust.ys[0, :], label="Robust MPC")
plt.legend(); plt.show()
```
By constructing an MPC `mpc_robust`, we see that the constraint are satisfied despite the disturbance $w$. One can see that the resulting response is quite conservative. This is because the constraint are tightened to handle the worst-case noise realization. For the example and scenario in question, the worst-case disturbance is for ${w = \left(\begin{smallmatrix}0.005 \\ 0.005 \end{smallmatrix}\right)}$. If we rerun the simulations for the worst-case disturbance, we get a large violation for `mpc_nominal`, while `mpc_robust` still manages to fulfill the constraints.

```@tabexample robust_mpc
# julia
worst_case_dynamics = (x, u, d) -> F*x + G*u + 0.005*ones(2)
sim_nominal_wc = Simulation(worst_case_dynamics, mpc_nominal; r=[0.5])
sim_robust_wc  = Simulation(worst_case_dynamics, mpc_robust;  r=[0.5])
hline([0.5], label="Constraint bound", linestyle=:dash)
plot!(sim_nominal_wc.ys[1,:], xlabel="Time step", ylabel="Position [m]", label="Nominal MPC")
plot!(sim_robust_wc.ys[1,:],  label="Robust MPC")
# exec-only
ylims!(0.4, 0.6)
# python
worst_case_dynamics = lambda x, u, d: F @ x + G @ u + 0.005 * np.ones(2)
sim_nominal_wc = Simulation(mpc_nominal, f=worst_case_dynamics, r=[0.5])
sim_robust_wc  = Simulation(mpc_robust,  f=worst_case_dynamics, r=[0.5])
plt.axhline(0.5, linestyle="--", label="Constraint bound")
plt.plot(sim_nominal_wc.ys[0, :], label="Nominal MPC")
plt.plot(sim_robust_wc.ys[0, :],  label="Robust MPC")
plt.xlabel("Time step"); plt.ylabel("Position [m]")
plt.legend(); plt.show()
```
!!! note "Measurable disturbance"
    Note that if the disturbance $w$ is measurable/known, it can be accounted for without having to tighten the constraints by providing it as a measurable disturbance $d$ (see `Model` for details.)

### Disturbances through an input matrix
A disturbance that does not act on every state variable independently is described by an input matrix $B_w$,
```math
x_{k+1} = F x_k + G u_k + B_w w_k, \qquad w_{\min} \leq w_k \leq w_{\max},
```
given with the keyword argument `Bw` (`nx × nw`, the identity by default). A disturbance on the control input, for example, enters through $B_w = G$:
```julia
set_disturbance!(mpc_robust, [-0.05], [0.05]; Bw = G)
```
A box on $G w$ also contains this disturbance, but it treats the state variables as independent and therefore tightens the constraints more than the input matrix does.

## Uncertainty in the current state
Uncertainty in the current state can also be handled robustly. If the current state $\hat{x}$ is assumed to be in the box 
```math
-\delta \leq x - \hat{x} \leq \delta,
```
constraints can be tightened with 

```@tab
# julia
set_x0_uncertainty!(mpc,delta)
# python
mpc.set_x0_uncertainty(delta)
```
`delta` is a vector with one entry per state variable, or a scalar applied to every state variable.

## The tightening
The control applied at time step $k$ of the prediction is $u_k = v_k - K x_k$, where $v_k$ are the decision variables and $K$ is the prestabilizing feedback (`set_prestabilizing_feedback!`, zero otherwise). The control $u_0$ is computed from the estimate $\hat{x}_0$ and is therefore exact. The deviation $e_k$ of the state from its nominal prediction then evolves as
```math
e_1 = F e_0 + B_w w_0, \qquad e_{k+1} = \Phi e_k + B_w w_k, \qquad \Phi = F - G K,
```
and a constraint row $a_x x_k + a_u u_k + a_{up} u_{k-1} \leq b$ deviates from its nominal prediction by a linear function of $e_0, w_0, \ldots, w_{k-1}$. Its upper bound is decreased by the worst case of this deviation over the boxes,
```math
|M_{k-1} F|\,\delta + \sum_{s=0}^{k-1} \max_{w_{\min} \leq w \leq w_{\max}} M_s B_w w,
\qquad M_s = (a_x - a_u K)\Phi^s - a_{up} K \Phi^{s-1},
```
where the last term is present for $s \geq 1$, and the lower bound is increased correspondingly. Constraints on $x_0$ are not tightened, since $x_0$ cannot be influenced. Without prestabilizing feedback, $\Phi = F$, and for an unstable or integrating system the tightening grows with the prediction horizon. Prestabilizing feedback keeps it bounded when $\Phi$ is stable.
