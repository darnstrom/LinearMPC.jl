# [Simple Example](@id man_simple)
A simplified form of the MPC problem that **LinearMPC.jl** solves is

```math
\begin{aligned}
        &\underset{u_0,\dots,u_{N-1}}{\text{minimize}}&& \textcolor{black}{\frac{1}{2}\sum_{k=0}^{N-1} {\left((Cx_{k}-r)^T Q (C x_{k}-r) + u_{k}^T R u_{k} + \Delta u_{k}^T R_r \Delta u_k\right)}}\\
        &\text{subject to} &&\textcolor{black}{{x_{k+1} = F x_k + G u_k}}, \quad k=0,\dots, N-1\\
        &&& \textcolor{black}{x_0 = \hat{x}} \\
        &&& \textcolor{black}{\underline{b} \leq A_x x_k + A_u u_k  \leq \overline{b}}, \quad k=0, \dots, N-1
\end{aligned}
```

Consider now the specific case when we have a given dynamical system $x_{k+1} = F x_k + G u_k$ with

```math
F \triangleq \begin{bmatrix}
1 & 0.5  \\
0 & 1
\end{bmatrix},\quad
G \triangleq \begin{bmatrix}
0   \\
1
\end{bmatrix}.
```

The quantities that we want to control are $y_1 = x_1$, and $y_2 = x_1+x_2$, which can be put on the form $y = C x$ with

```math
C \triangleq  \begin{bmatrix}
1 & 0    \\
1 & 1
\end{bmatrix},
```

and the weights of the objective are

```math
Q \triangleq  \begin{bmatrix}
1 & 0    \\
0 & 1
\end{bmatrix},\quad
R = 0, \quad
Rr= 1, \quad
```

Moreove, we assume that we have the input constraint $-3 \leq u \leq 3$, and the output constraints  $0 \leq y_1 \leq 1$, and $0 \leq y_2 \leq 2$.
Finally we further assume that we have a constraint that $-1 \leq 2 x_1 - x_2 \leq 2$.


```julia
using LinearMPC
# Dynamics
F = [1 0.5; 0 1]
G = [0;1]

# create mpc struct
mpc  = LinearMPC.MPC(F,G; C=[1 0; 1  1])

# objective
set_objective!(mpc, Q=[1,1], R=0, Rr = [1])

# add constraints
set_bounds!(mpc,umin=[-3],umax=[3], ymin= [0, 0],ymax = [1,2])
add_constraints!(mpc,Ax = [2 -1], lb = [-1], ub = [2])

# the prediction/control horizons
set_horizon!(mpc, Np=10, Nc=10)

```
That is it! No we are ready to try the controller out.

### Testing the MPC controller in simulation
The function `compute_control(mpc,x;r)` computes a control action given the state `x` and reference value `r`. For example,
```julia
u = compute_control(mpc,[1,2];r=[0,0])
```
gives that the optimal control action at `x=[1,2]` with `r=[0,0]` is `u`.

Still, it is hard based on solving just one problem to get a feel for if the controller is doing what we desire. LinearMPC.jl can be used to create a `Simulation` for the given controller.

```julia
# simulate the system
sim = LinearMPC.Simulation(mpc;x0=[0,0],r=[1,0])

```
To visualize the result we can use `Plots`:

```julia
using Plots
plot(sim)
```

We can make the value of $y_1$ reach its reference value quicker by weighting it more:

```julia
set_objective!(mpc, Q=[1e3 1], Rr = [1])
sim = LinearMPC.Simulation(mpc;x0=[0,0],r=[1,0])
plot(sim)
```

This is closer toof the desired behaviour that we are after! Now we are read to to deploy the code on our embedded system.

### Code generation
To generate corresponding C-code for the MPC controller, we use the function `codegen`

```julia
codegen
```
