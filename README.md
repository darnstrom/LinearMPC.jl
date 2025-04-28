# **LinearMPC.jl**
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![](https://img.shields.io/badge/docs-online-brightgreen)](https://darnstrom.github.io/LinearMPC.jl)

**LinearMPC.jl** is a Julia package for Model Predictive Control (MPC) of linear systems. It aims to produce _high-performant_ and _lightweight_ C-code that can easily be used on embedded systems, while at the same time give a user-friendly and expressive development environment for MPC. The package supports code generation for the Quadratic Programming solver [DAQP](https://github.com/darnstrom/daqp), and for explicit solutions computed by [ParametricDAQP.jl](https://github.com/darnstrom/ParametricDAQP.jl). 

A simplified version (see the documentation for a more complete formulation) of the solved problem is

$$
\begin{align}
        &\underset{u_0,\dots,u_{N-1}}{\text{minimize}}&& \frac{1}{2}\sum_{k=0}^{N-1} {\left((Cx_{k}-r)^T Q (C x_{k}-r) + u_{k}^T R u_{k} + \Delta u_{k}^T R_r \Delta u_k\right)}\\
        &\text{subject to} &&{x_{k+1} = F x_k + G u_k}, \quad k=0,\dots, N-1\\
        &&& x_0 = \hat{x} \\
        &&& {\underline{b} \leq A_x x_k + A_u u_k  \leq \overline{b}}, \quad k=0, \dots, N-1
        \end{align}
$$

where $\hat{x}$ is the current state and $r$ is the desired reference value of $Cx$.


## Example 
The following code show a simple MPC example of controlling an inverted pendulum on a cart, inspired by [this](https://se.mathworks.com/help/mpc/ug/mpc-control-of-an-inverted-pendulum-on-a-cart.html) example in the Model Predictive Toolbox in MATLAB.
```julia
using LinearMPC
# Continuous time system dx = A x + B u, y = C x
A = [0 1 0 0; 0 -10 9.81 0; 0 0 0 1; 0 -20 39.24 0]; 
B = 100*[0;1.0;0;2.0;;];
C = [1.0 0 0 0; 0 0 1.0 0];


# create an MPC control with sample time 0.01, prediction/control horizon 50/5
Ts = 0.01
mpc = LinearMPC.MPC(A,B,Ts;C,Np=50,Nc=5);

# set the objective functions weights
set_objective!(mpc;Q=[1.2^2,1], R=[0.0], Rr=[1.0])

# set actuator limits
set_bounds!(mpc; umin=[-2.0],umax=[2.0])
# additional functions for adding constraints: set_output_bounds!, add_constraint!
```

A control, given the state `x` and reference value `r`, is computed with
```julia
# compute control u at state x = [0,0,0,0] and reference r. 
u = compute_control(mpc,[0,0,0,0], r = [1, 0])
```

Embeddable C-code for the MPC controller is generated with the command
```julia
LinearMPC.codegen(mpc;dir="codgen_dir")
```
which produces _allocation-free_ C-code in the directory `codegen_dir` for setting up optimization problems and solving them with the Quadratic Programming solver [DAQP](https://github.com/darnstrom/daqp).

The C-function `mpc_compute_control(control, state, reference, disturbance)` computes the optimal control, given the current `state`, `reference`, and measured disturbances `disturbance`, which are all floating-point arrays. The optimal control is stored in the floating-point array `control`. 


## Citation
If you find the package useful, consider citing one of the following papers, which are the backbones of the package:
```
@article{arnstrom2022daqp,
  author={Arnström, Daniel and Bemporad, Alberto and Axehill, Daniel},
  journal={IEEE Transactions on Automatic Control},
  title={A Dual Active-Set Solver for Embedded Quadratic Programming Using Recursive {LDL}$^{T}$ Updates},
  year={2022},
  volume={67},
  number={8},
  pages={4362-4369},
  doi={10.1109/TAC.2022.3176430}
}
```

```
@inproceedings{arnstrom2024pdaqp,
  author={Arnström, Daniel and Axehill, Daniel},
  booktitle={2024 IEEE 63rd Conference on Decision and Control (CDC)}, 
  title={A High-Performant Multi-Parametric Quadratic Programming Solver}, 
  year={2024},
  volume={},
  number={},
  pages={303-308},
}
```
