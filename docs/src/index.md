# LinearMPC.jl

!!! note "Documentation under development"
    The documentation for LineaerMPC.jl is currently under development

The focus of **LinearMPC.jl** is Model Predictive Control (MPC) of linear systems. The aim of the package is to produce _high-performant_ and _lightweight_ C-code that can easily be used on embedded systems, while at the same time give a user-friendly and expressive development environment for MPC. The package supports code generation for the Quadratic Programming solver [DAQP](https://github.com/darnstrom/daqp), and for explicit solutions computed by [ParametricDAQP.jl](https://github.com/darnstrom/ParametricDAQP.jl).

## Model Predictive Control 
In MPC, an optimal control decision is computed at every sampling instance by solving an _optimization probelm_. On a high-level, the optimization problems solved are of the form

```math
\begin{aligned}
        &\underset{u_0,\dots,u_{N-1}}{\text{minimize}}&& \textcolor{purple}{\frac{1}{2}\sum_{k=0}^{N-1} {\left((Cx_{k}-r)^T Q (C x_{k}-r) + u_{k}^T R u_{k} + \Delta u_{k}^T R_r \Delta u_k\right)}}\\
        &\text{subject to} &&\textcolor{blue}{{x_{k+1} = F x_k + G u_k}}, \quad k=0,\dots, N-1\\
        &&& \textcolor{red}{x_0 = \hat{x}} \\
        &&& \textcolor{green}{\underline{b} \leq A_x x_k + A_u u_k  \leq \overline{b}}, \quad k=0, \dots, N-1
\end{aligned}
```
where an $\textcolor{purple}{\text{objective}}$ is minimized , subject to a $\textcolor{blue}{\text{dynamical system}}$ that is simulated over a horizon $N$, $\textcolor{red}{\text{starting from}}$ an estimate $\hat{x}$ of the current state. Additionally, $\textcolor{green}{\text{constraints}}$ like actuator limits and state constraints are accounted for. The objective is comprised by the deviation of an output $y= Cx" from a reference value $r$, the control effort ($u^T R u$) , and the change of the control action $\Delta u^T R_r \Delta u$. 

LinearMPC.jl generates a _condensed_ problem by eliminating the equality constraint from the dynamics. The resuliting optimization problem is a dense Quadratic Program (QP). *LinearMPC.jl* uses the QP sovler [DAQP](https://github.com/darnstrom/daqp), a dual active-set solver that has been specialized to solved such problems.

One can also explicitly express the solution map problem as a piecewise affine function over polyhedral regions. *LinearMPC.jl* supports the computation of such explicit solutions by interfacing the multi-parameteric QP solver [ParametricDAQP.jl](https://github.com/darnstrom/ParametricDAQP.jl.)

## Why LinearMPC.jl? 
* Code generation of **high-performant**, **allocation-free**, **library-free**, and **lightweight** C-code that can be embedded on _any_ micro controller. 
* State-of-the computation of explicit solutions (**~100x** faster than other software packages)  
* Tools to determine **real-time certificates** of the complexity of the solver, allowing for MPC in with guarantees on the memory and computationala requirements before deploying the solver.

## Why not LinearMPC.jl? 
As its name suggests, the package is specialized for MPC for _linear_ systems. If a linear (or linearized) model does not suffices for your use case, consider the following packages.
* [ModelPredictiveControl.jl](https://github.com/JuliaControl/ModelPredictiveControl.jl) - A high-level MPC packages in Julia. While it does not generate embeddable C-code, it is an excellent package during development of MPC controllers. 
* [acados](https://github.com/acados/acados) - provides fast and embedded solvers for nonlinear optimal control, specifically designed for real-time applications and embedded systems. This is the current state-of-the art if you are interested in real-time nonlinear MPC on embedded systems. 
* Are you more of a Python person? You can still use LinearMPC.jl through its sister package [lmpc](https://github.com/darnstrom/lmpc).
