# **LinearMPC.jl**
**LinearMPC.jl** is a Julia package for Model Predictive Control (MPC) of linear systems. Its aim is to produce _high-performant_ and _lightweight_ C-code that can easily be used on embedded systems. The package supports code generation for the Quadratic Programming solver [DAQP](https://github.com/darnstrom/daqp), and for explicit solutions computed by [ParametricDAQP.jl](https://github.com/darnstrom/ParametricDAQP.jl). 

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
# Continuous time system dx = A x + B u
A = [0 1 0 0; 0 -10 9.81 0; 0 0 0 1; 0 -20 39.24 0]; 
B = 100*[0;1.0;0;2.0;;];
C = [1.0 0 0 0; 0 0 1.0 0];


# create an MPC control with sample time 0.01, prediction horizon 10 and control horizon 5 
Np,Nc = 10,5
Ts = 0.01
mpc = LinearMPC.MPC(A,B,Ts;C,Nc,Np);

# set the objective functions weights
Q,R,Rr= [1.2^2,1], [0.0], [1.0]
set_weights!(mpc;Q,R,Rr)

# set actuator limits
umin,umax = [-2.0], [2.0]
set_bounds!(mpc; umin,umax)
```

A control, given the state `x` and reference value `r`, is computed with
```julia
r = [1;0]
x = [0;0;0;0]
u = compute_control(mpc,x;r)
```

Embeddable C-code for the MPC controller is generated with the command
```julia
LinearMPC.codegen(mpc;dir="codgen_dir")
```
which produce _allocation-free_ C-code in the directory `codegen_dir` for setting up optimization problems and solving them with the Quadratic Programming solver [DAQP](https://github.com/darnstrom/daqp).

The C-function `mpc_compute_control(control, state, reference, in_dist, out_dist)` computes the optimal control, given the current `state`, `reference`, and measured intput/output disturbances `in_dist`/`out_dist`, which are all floating-point arrays. The optimal control is stored in the floating-point array `control`. 


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
