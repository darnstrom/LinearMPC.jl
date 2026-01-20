# [Game-Theoretic MPC](@id man_robust)

**LinearMPC.jl** can also solve general Nash equilibria to game-theoretic linear MPC problems.
Specifically it can handle objective of the form 

```math
\forall i \in \mathcal{I}: \min_{u^i} J^i(u^i,x), 
```
where $J^i$ is a linear-quadratic objective.
Moreover, $\{u^i\}_{i\in\mathcal{I}}$ completely partitions the set of nominal controls $u$.

## Adding objectives 
For adding an objective, one can use the function `set_objective!` as usual, but with an additional argument `uids` which specify which controls corresponding to the agents objective. 

For example, if there are 5 controls and 2 players, where player 1 governs control 1,4, and player 2 governs control 2,3,5, one can set up the objectives as
```julia
set_objective!(mpc,[1,4];Q=Q1, R=R1, Rr=Rr1)
set_objective!(mpc,[2,3,5];Q=Q2, R=R2, Rr=Rr2)
```
where `Qi`, `Ri`, and `Rri`, corresponds to the cost weights for player `i` on states, controls, and control change, respectively.

## Example  
The following example illustrate how generalized Nash equilibria can be computed for a double integrator with two inputs, where each input is control by its own player. 

The objective for Player 1 is: 

```math
J^1(u^1,x) = \sum_{k=0}^{9} {\left((x_{k}-r)^T Q^1 (x_{k}-r) + 10^3\|\Delta u^1_{k}\|^2_2 \right)},
```

with $Q^1 = \begin{bmatrix}
    1& 0  \\
    0& 0
\end{bmatrix}$ and $u^1$ being the first control.

The objective for Player 2 is: 
```math
J^2(u^2,x) = \sum_{k=0}^{9} {\left((x_{k}-r)^T Q^2 (x_{k}-r) + 10^3\|\Delta u^2_{k}\|^2_2 \right)},
```

with $Q^2 = \begin{bmatrix}
    0& 0  \\
    0& 1 
\end{bmatrix}$ and $u^2$ being the second control.

In additions to this we have the constraints $\|u\|_{\infty} \leq 1$.

This problem can be setup with 

```@example game_mpc
using LinearMPC
F,G = [1 0.1; 0 1], [0 0;1 1];
mpc = LinearMPC.MPC(F,G;C=[1 0;0 1],Np=10);

set_objective!(mpc, [1]; Q=[1,0], Rr=1e3);
set_objective!(mpc, [2]; Q=[0,1], Rr=1e3);

set_bounds!(mpc;umin=-ones(2),umax=ones(2));
nothing # hide

```

As usual, the closed-loop behaviour for a scenario can be simulated with

```@example game_mpc
sim_game = Simulation(mpc;x0=10*ones(2), r = [10,0], N=100);

using Plots
plt_game = plot(sim_game, label="Game-theoretic MPC")
```

The first player tries to drive the first state (the position) to 10, while the second player tries to drive the second state (the velocity) to 0.

We can also compare the game-theoretic closed-loop behaviour with a cooperative centralized MPC controller: 

```@example game_mpc
F,G = [1 0.1; 0 1], [0 0;1 1];
mpc = LinearMPC.MPC(F,G;C=[1 0;0 1],Np=10);
set_objective!(mpc; Q=[1,1], Rr=1e3);
set_bounds!(mpc;umin=-ones(2),umax=ones(2));

sim_centralized = Simulation(mpc;x0=10*ones(2), r = [10,0], N=100);
plot!(plt_game, sim_centralized, label="Centralized MPC", color=:red)
plot!(subplot=1,legend=true) # hide
```

As is expected, the centralized MPC leads to a better reference tracking since it is coordinating both control signals at the same time, while the players in the game-theoretic MPC are greedily trying to fulfill their own objectives (which are somewhat conflicting, since player 1 tries to drive the position to 10, while player 2 is try to bring the system to rest with zero velocity.) 
