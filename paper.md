---
title: 'LinearMPC.jl: A Julia Package for Embedded Linear Model Predictive Control'
tags:
  - Julia
  - model predictive control
  - embedded systems
  - code generation
  - quadratic programming
  - control systems
authors:
  - name: Daniel Arnström
    orcid: 0000-0003-0970-0620
    corresponding: true
    affiliation: 1
affiliations:
  - name: Independent Researcher, Sweden 
    index: 1
date: 24 March 2026
bibliography: paper.bib
---

# Summary

`LinearMPC.jl` is a Julia package for Model Predictive Control (MPC) of linear systems. It provides a user-friendly development environment for designing, simulating, and deploying MPC controllers while targeting high-performance embedded implementations. The package bridges the gap between rapid prototyping in a high-level language and deployment on resource-constrained hardware by generating allocation-free, library-free C code that can run on any microcontroller. `LinearMPC.jl` supports both online optimization via the dual active-set Quadratic Programming (QP) solver DAQP [@Arnstrom:2022] and explicit MPC solutions computed by the multi-parametric QP solver ParametricDAQP [@Arnstrom:2024]. Additional features include robust MPC with constraint tightening, state estimation via Kalman filters, hybrid MPC with binary controls, game-theoretic MPC for Nash equilibria, and complexity certification for real-time guarantees.

# Statement of need

Model Predictive Control is the dominant advanced control strategy in industry, used in applications ranging from process control to autonomous vehicles and robotics. At each sampling instant, MPC solves an optimization problem that accounts for a dynamical model, constraints, and a performance objective. For linear systems, this optimization problem is a QP. Deploying MPC on embedded systems, such as microcontrollers in automotive, aerospace, and robotic applications, poses several challenges: the controller must execute within strict timing deadlines, the code must be lightweight and allocation-free, and the implementation must be verifiable.

Existing tools address parts of this workflow. `ModelPredictiveControl.jl` [@ModelPredictiveControl] provides a high-level Julia interface for MPC design but does not target embedded deployment. `acados` [@Verschueren:2022] provides fast embedded solvers for nonlinear optimal control. MATLAB's Model Predictive Control Toolbox offers a mature commercial solution. However, there is a gap for a tool that combines (i) an expressive, high-level Julia interface for rapid prototyping, (ii) generation of lightweight, allocation-free C code for any embedded target, (iii) state-of-the-art explicit MPC computation, and (iv) complexity certification for hard real-time guarantees. `LinearMPC.jl` fills this gap by providing a unified workflow from design to deployment for linear MPC.

# Software design

The architecture of `LinearMPC.jl` follows a layered design (\autoref{fig:architecture}):

1. **Model layer.** Users define system dynamics as discrete-time or continuous-time state-space models, transfer functions (via `ControlSystems.jl` integration), or nonlinear models that are automatically linearized using `ForwardDiff.jl`. Measured disturbances are also supported.

2. **MPC formulation layer.** An expressive API lets users configure the controller: setting objectives with tracking, rate-of-change, and linear cost terms; defining input, output, and general polyhedral constraints (hard or soft); specifying prediction and control horizons with move blocking; and adding prestabilizing feedback for numerical conditioning.

3. **Optimization layer.** The MPC problem is condensed into a dense QP and solved using the DAQP solver [@Arnstrom:2022], a dual active-set method specialized for embedded applications. For hybrid systems with binary controls, DAQP's branch-and-bound capabilities are used. Game-theoretic MPC problems are solved by computing Nash equilibria of coupled QPs.

4. **Explicit MPC layer.** The parametric solution of the QP, a piecewise affine function over polyhedral regions [@Bemporad:2002], can be computed using ParametricDAQP [@Arnstrom:2024], which is approximately 100 times faster than competing solvers. A binary search tree [@Tondel:2003] is constructed for efficient online point location.

5. **Deployment layer.** The `codegen` function generates standalone, allocation-free C code for both online and explicit MPC controllers, including state observers when configured. The generated code depends on no external libraries and can be compiled for any target with a C compiler.

6. **Verification layer.** Complexity certification via the `ASCertain` package [@Arnstrom:2022b] provides worst-case iteration bounds for the DAQP solver, enabling hard real-time guarantees before deployment.

![Simplified workflow of `LinearMPC.jl`: from model definition and MPC design in Julia to embedded C code deployment.\label{fig:architecture}](architecture.pdf)

The central optimization problem solved at each sampling instant is
\begin{equation}\label{eq:mpc}
\begin{aligned}
    &\underset{u_0,\dots,u_{N-1}}{\text{minimize}}&& \frac{1}{2}\sum_{k=0}^{N-1} {\left((Cx_{k}-r)^T Q (C x_{k}-r) + u_{k}^T R u_{k} + \Delta u_{k}^T R_r \Delta u_k + l_k^T u_k\right)}\\
    &\text{subject to} && x_{k+1} = F x_k + G u_k, \quad k=0,\dots, N-1\\
    &&& x_0 = \hat{x} \\
    &&& \underline{b} \leq A_x x_k + A_u u_k  \leq \overline{b}, \quad k=0, \dots, N-1
\end{aligned}
\end{equation}
where $\hat{x}$ is the current state estimate and $r$ is the desired reference. The objective weights $Q$, $R$, $R_r$ and the constraint matrices $A_x$, $A_u$ are configured by the user. Optional terminal costs, cross terms, reference preview, and time-varying linear costs extend this basic formulation.

# Key features

- **Rapid prototyping.** Controllers are designed and tested in Julia with a concise, expressive API. Closed-loop simulations and plotting are built in.
- **Embedded code generation.** Allocation-free, library-free C code is generated for both online and explicit MPC, including Kalman filter observers.
- **Explicit MPC.** State-of-the-art computation of piecewise affine control laws, with visualization of critical regions and feedback surfaces.
- **Robust MPC.** A priori constraint tightening handles bounded process noise and state estimation uncertainty, with guarantees that constraints are satisfied under worst-case disturbances.
- **Hybrid MPC.** Binary control variables are supported for systems with on/off actuators, solved via branch-and-bound within DAQP.
- **Game-theoretic MPC.** Generalized Nash equilibria can be computed for multi-agent linear MPC problems where each agent has its own objective.
- **Complexity certification.** Worst-case solver iteration counts can be certified offline, enabling hard real-time guarantees.
- **Numerical conditioning.** Prestabilizing feedback reparametrization mitigates ill-conditioning for unstable systems with long prediction horizons.

# Example usage

The following example demonstrates the core workflow: defining a system, designing an MPC controller, simulating, and generating embedded C code:

```julia
using LinearMPC

# Continuous-time inverted pendulum on a cart
A = [0 1 0 0; 0 -10 9.81 0; 0 0 0 1; 0 -20 39.24 0]
B = 100*[0; 1.0; 0; 2.0;;]
C = [1.0 0 0 0; 0 0 1.0 0]

# Create MPC with sample time 0.01, horizon 50/5
mpc = LinearMPC.MPC(A, B, 0.01; C, Np=50, Nc=5)
set_objective!(mpc; Q=[1.2^2, 1], R=[0.0], Rr=[1.0])
set_bounds!(mpc; umin=[-2.0], umax=[2.0])

# Compute a control action
u = compute_control(mpc, [0,0,0,0], r=[1, 0])

# Simulate closed-loop and generate C code
sim = Simulation(mpc; x0=zeros(4), r=[1,0], N=100)
LinearMPC.codegen(mpc; dir="mpc_codegen")
```

# AI usage disclosure

Generative AI (GitHub Copilot, powered by Claude) was used to assist in drafting this manuscript. The core software design decisions, algorithmic implementations, and all technical content were made by the human author. All AI-assisted output was reviewed and edited by the author.

# References
