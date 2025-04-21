# **LinearMPC.jl**
**LinearMPC.jl** is a julia package for Model Predictive Control (MPC) of linear systems. Its aim is to produce _high-performant_ and _lightweight_ C-code that can easily be used on embedded systems. The package both support code generation for the Quadratic Programming solver [DAQP](https://github.com/darnstrom/daqp), and for explicit solution computed by [ParametricDAQP.jl](https://github.com/darnstrom/ParametricDAQP.jl). 

A simplified version (see the documentation for a more complete formulation) of the problem solved is

$$
\begin{align}
        &\underset{\mathbf{u}}{\text{minimize}}&& \frac{1}{2}\sum_{k=0}^{N_p-1} {\left(x_{k}^T Q x_{k} + u_{k}^T R u_{k} + \Delta u_{k}^T R_r \Delta u_k\right)}\\
        &\text{subject to} &&\textcolor{red}{x_{k+1} = F x_k + G u_k}, \quad k=0,\dots, N_p-1\\
        &&& x_0 = x \\
        &&& {\underline{b} \leq A_x x_k + A_u u_k  \overline{b}}, \quad k=0, \dots, N-1
        \end{align}
$$


## Example 
The following code show a simple MPC example of controlling an inverted pendulum on a cart, inspired by [this](https://se.mathworks.com/help/mpc/ug/mpc-control-of-an-inverted-pendulum-on-a-cart.html) example in the Model Predictive Toolbox in MATLAB.
```julia
using LinearMPC
# Continuous time system dx = A x + B u
A = [0 1 0 0; 0 -10 9.81 0; 0 0 0 1; 0 -20 39.24 0]; 
B = 100*[0;1.0;0;2.0;;];
C = [1.0 0 0 0; 0 0 1.0 0];


# create an MPC control with sample time 0.01, prediction horizon 10 and control horizon 5 
Np,Nc = 10,5
Ts = 0.01;
mpc = LinearMPC.MPC(A,B,Ts;C,Nc,Np);

# set the objective functions weights
Q,R,Rr= [1.2^2,1], [0.0], [1.0]
set_weights!(mpc;Q,R,Rr)

# set actuator limits
umin,umax = [-2.0], [2.0]
set_bounds!(mpc; umin,umax)
```

A control state `x` and setpoint `r` can then be computed with 
```julia
r = [1;0];
x = [0;0;0;0]
u = compute_control(mpc,x;r)
```

The following code then generates embeddable C-code for the MPC controller
```julia
LinearMPC.codegen(mpc;dir="codgen_dir")
```
The directory `codegen_dir` then contain C-code for setting up optimization problems and for solving them with the Quadratic Programming solver [DAQP](https://github.com/darnstrom/daqp).

The C-function `mpc_compute_control(control, state, reference, disturbance)` populates the floating-point array `control` with the optimal control given the current `state`, `reference`, and measured `disturbance` (all of which are also floating-point arrays).


## Citation
If you find the package useful consider citing one of the following papers, which are the backbones of the package: 
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
